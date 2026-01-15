// resolve_symmetry_and_cell(json& rtdbjson)
//
// Goal:
//   - Decide the *working* (effective) simulation cell + symmetry to be used by physics.
//   - Keep user intent in geometries[geomname]["symmetry"] (and coords/unita) untouched.
//   - Write derived, self-consistent results into rtdbjson["effective_symmetry"] (+ optionally effective_geometry).
//
// Key precedence rule (as you stated):
//   - If rtdbjson["nwpw"]["simulation_cell"]["unita"] is present -> it wins.
//   - Otherwise use rtdbjson["geometries"][geomname]["unita"].
//
// Symmetry policy:
//   - If geometry symmetry specified => start from that.
//   - Else if autosym is true => detect symmetry from lattice+coords.
//   - Else => trivial symmetry (identity only).
//
// Primitive conversion policy:
//   - If primitive requested (from geometry.symmetry.primitive or a global switch),
//     then conventional->primitive conversion overwrites *effective* lattice, coords, and ops.
//
// Restart policy (recommended):
//   - If rtdbjson already contains effective_symmetry and a restart flag indicates reuse,
//     keep it unless user forced a new symmetry_specified or changed primitive/cell.
//
// NOTE:
//   This function calls helper hooks you’ll implement (or wrap spglib):
//     - build_symmetry_from_spec(...)
//     - autosym_detect(...)
//     - make_primitive_cell(...)
//   I’ve left them as clear placeholders so you can wire your Symmetry constructor later.

#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <array>
#include <iomanip>
#include <iostream>

#include "json.hpp"
using json = nlohmann::json;

#include "Symmetry.hpp"
#include "autospace.hpp"

#include "apply_kpoint_symmetry_pruning.hpp"

namespace pwdft {


/**************************************************
 *                                                *
 *                 read_bool_flag                 *
 *                                                *
 **************************************************/
/**
 * @brief Read a boolean flag from JSON.
 *
 * Accepts:
 *  - boolean values (true / false)
 *  - integer values (0 = false, nonzero = true)
 *
 * If the key does not exist or has an unsupported type,
 * the default value is returned.
 *
 * @param j        JSON object
 * @param key      Key to read
 * @param defval   Default value if key is missing or invalid
 * @return         Parsed boolean value
 */
static bool read_bool_flag(const json& j,
                           const std::string& key,
                           bool defval = false)
{
    if (!j.contains(key))
        return defval;

    const auto& v = j.at(key);

    if (v.is_boolean())
        return v.get<bool>();

    if (v.is_number_integer())
        return (v.get<int>() != 0);

    return defval;
}



/**************************************************
 *                                                *
 *                 invert_3x3                     *
 *                                                *
 **************************************************/
/**
 * @brief Invert a 3×3 matrix.
 *
 * Computes the inverse of a 3×3 matrix using the analytic adjugate /
 * determinant formula. The matrix is assumed to be stored in
 * row-major order:
 *
 *     M = | a b c |
 *         | d e f |
 *         | g h i |
 *
 * This routine is primarily intended for lattice and reciprocal-lattice
 * transformations, where numerical stability and determinism are required.
 *
 * @param M     Input 3×3 matrix (row-major).
 * @param Minv  Output 3×3 inverse matrix (row-major).
 *
 * @throws std::runtime_error if the matrix is singular or nearly singular
 *         (|det| < 1e-14).
 *
 * @note This function performs no pivoting and assumes a well-conditioned
 *       lattice matrix. Callers should ensure that `M` represents a valid
 *       simulation cell or metric tensor.
 */
static void invert_3x3(const double M[9], double Minv[9])
{
    const double a = M[0], b = M[1], c = M[2];
    const double d = M[3], e = M[4], f = M[5];
    const double g = M[6], h = M[7], i = M[8];

    const double det =
          a*(e*i - f*h)
        - b*(d*i - f*g)
        + c*(d*h - e*g);

    if (std::abs(det) < 1e-14)
        throw std::runtime_error("invert_3x3: singular lattice");

    const double invdet = 1.0 / det;

    Minv[0] =  (e*i - f*h) * invdet;
    Minv[1] = -(b*i - c*h) * invdet;
    Minv[2] =  (b*f - c*e) * invdet;

    Minv[3] = -(d*i - f*g) * invdet;
    Minv[4] =  (a*i - c*g) * invdet;
    Minv[5] = -(a*f - c*d) * invdet;

    Minv[6] =  (d*h - e*g) * invdet;
    Minv[7] = -(a*h - b*g) * invdet;
    Minv[8] =  (a*e - b*d) * invdet;
}


/**************************************************
 *                                                *
 *              make_primtive_cell                *
 *                                                *
 **************************************************/
/**
 * @brief Convert a centered conventional cell into a primitive cell.
 *
 * Constructs a primitive lattice and corresponding atomic coordinates
 * from a centered space-group description. This routine is only applied
 * when a space group with centering (I, F, etc.) is detected and the user
 * has explicitly requested a primitive cell.
 *
 * The procedure:
 *  1. Identifies centering translations from identity-rotation symmetry
 *     operators.
 *  2. Constructs primitive lattice vectors using standard International
 *     Tables conventions (currently I- and F-centering).
 *  3. Rewrites the lattice vectors in-place.
 *  4. Transforms Cartesian coordinates into the primitive basis and
 *     wraps them into the unit cell.
 *
 * This transformation is *structural*, not merely cosmetic: it reduces
 * the simulation cell volume and implicitly changes the Brillouin zone.
 * Downstream k-point generation and symmetry pruning must therefore be
 * performed **after** this routine is applied.
 *
 * @param sym        Effective space-group symmetry.
 * @param unita      On input: 3×3 lattice matrix of the conventional cell
 *                   (row-major).
 *                   On output: 3×3 lattice matrix of the primitive cell.
 * @param coords_xyz Atomic coordinates (Cartesian, flat 3N array).
 *                   On output, coordinates are expressed in the primitive
 *                   lattice basis and wrapped into [0,1).
 *
 * @note
 *  - This routine currently supports I- and F-centered lattices.
 *    Other centerings (A, B, C, R) may be added later.
 *  - If the symmetry has no centering (P-lattice) or is not a space group,
 *    the function returns without modification.
 *  - No attempt is made to re-identify symmetry after reduction; the
 *    resulting primitive cell is assumed consistent with the original
 *    space-group intent.
 *
 * @warning
 *  Primitive-cell construction changes reciprocal-space topology.
 *  Any cached k-points or Brillouin-zone data must be regenerated.
 */
static void make_primitive_cell(const Symmetry& sym, double unita[9], std::vector<double>& coords_xyz)
{
    // Only meaningful for space groups with centering
    if (!sym.is_space_group() || sym.num_centering() == 1)
        return;

    // --------------------------------------------------
    // Step 1: determine centering type from translations
    // --------------------------------------------------
    // We look at symmetry ops with identity rotation
    // and extract their fractional translations.

    std::vector<std::array<double,3>> centering;

    for (const auto& op : sym.operators())
    {
        bool is_identity = true;
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                if (op.R[i][j] != (i==j))
                    is_identity = false;

        if (!is_identity) continue;

        centering.push_back({op.t[0], op.t[1], op.t[2]});
    }

    if (centering.size() <= 1)
        return;

    // --------------------------------------------------
    // Step 2: build primitive lattice vectors
    // --------------------------------------------------
    // Minimal hard-coded centering reductions
    // (this matches International Tables)

    double A[3] = { unita[0], unita[1], unita[2] };
    double B[3] = { unita[3], unita[4], unita[5] };
    double C[3] = { unita[6], unita[7], unita[8] };

    double Ap[3], Bp[3], Cp[3];

    if (sym.num_centering() == 2)
    {
        // I-centering (body-centered)
        for (int i=0;i<3;i++)
        {
            Ap[i] = 0.5*(B[i] + C[i]);
            Bp[i] = 0.5*(A[i] + C[i]);
            Cp[i] = 0.5*(A[i] + B[i]);
        }
    }
    else if (sym.num_centering() == 4)
    {
        // F-centering
        for (int i=0;i<3;i++)
        {
            Ap[i] = 0.5*(B[i] + C[i]);
            Bp[i] = 0.5*(A[i] + C[i]);
            Cp[i] = 0.5*(A[i] + B[i]);
        }
    }
    else
    {
        // Other centerings (A/B/C/R) can be added later
        return;
    }

    // Write primitive unita
    for (int i=0;i<3;i++) unita[i]   = Ap[i];
    for (int i=0;i<3;i++) unita[3+i] = Bp[i];
    for (int i=0;i<3;i++) unita[6+i] = Cp[i];

    // --------------------------------------------------
    // Step 3: reduce coordinates into primitive cell
    // --------------------------------------------------
    // Convert cartesian → fractional (old cell)
    // then map to primitive basis and wrap into [0,1)

    // Build primitive lattice matrix
    double M[9] = {
        Ap[0], Ap[1], Ap[2],
        Bp[0], Bp[1], Bp[2],
        Cp[0], Cp[1], Cp[2]
    };

    double Minv[9];
    invert_3x3(M, Minv);

    std::vector<double> new_coords;

    for (size_t i=0;i<coords_xyz.size();i+=3)
    {
       const double x = coords_xyz[i];
       const double y = coords_xyz[i+1];
       const double z = coords_xyz[i+2];
 
       double xf = Minv[0]*x + Minv[1]*y + Minv[2]*z;
       double yf = Minv[3]*x + Minv[4]*y + Minv[5]*z;
       double zf = Minv[6]*x + Minv[7]*y + Minv[8]*z;
 
       xf -= std::floor(xf);
       yf -= std::floor(yf);
       zf -= std::floor(zf);
 
       new_coords.push_back(xf);
       new_coords.push_back(yf);
       new_coords.push_back(zf);
    }

    coords_xyz.swap(new_coords);

}



/**************************************************
 *                                                *
 *              symmetry_fingerprint              *
 *                                                *
 **************************************************/
/**
 * @brief Construct a deterministic fingerprint for an effective symmetry state.
 *
 * Builds a compact, deterministic string representation capturing the
 * *structural identity* of the resolved symmetry as seen by downstream physics.
 * The fingerprint is intended for:
 *
 *   - Restart consistency checks
 *   - Cache validation
 *   - Debugging and provenance tracking
 *
 * The fingerprint encodes:
 *   - Symmetry name (e.g. space-group symbol or "identity")
 *   - Group order
 *   - Whether a primitive cell was requested
 *   - Translation interpretation (e.g. "fractional")
 *   - The effective lattice matrix (`unita`)
 *
 * No symmetry operators are serialized; the fingerprint assumes that
 * symmetry operators are deterministically reconstructible from the
 * symmetry name and settings.
 *
 * @param sym                  Resolved symmetry object.
 * @param unita                Effective 3×3 lattice matrix (row-major).
 * @param primitive_requested  Whether primitive-cell reduction was requested.
 * @param translation_type     Interpretation of symmetry translations
 *                             (e.g. "fractional" or "cartesian").
 *
 * @return A deterministic string fingerprint uniquely identifying the
 *         effective symmetry + lattice state.
 *
 * @note This function is intentionally inexpensive and human-readable.
 *       Callers may hash the returned string if a compact identifier is needed.
 *       This is not intended to be a cryptographic hash.
 */
static std::string symmetry_fingerprint(const pwdft::Symmetry& sym,
                                        const double unita[9],
                                        bool primitive_requested,
                                        const std::string& translation_type)
{
   // Keep this deterministic and cheap.
   // If you already have a hash util, use it.
   std::ostringstream oss;
   oss.setf(std::ios::fixed); oss<<std::setprecision(12);
   oss << sym.name() << "|" << sym.order()
       << "|prim=" << (primitive_requested?1:0)
       << "|tt=" << translation_type << "|unita=";
   for (int i=0;i<9;++i) oss << unita[i] << ",";
   return oss.str(); // later: hash this string
}


/**************************************************
 *                                                *
 *                  get_geomname                  *
 *                                                *
 **************************************************/
// ------------------------------
// helper: safe get geometry name
// ------------------------------
static std::string get_geomname(const json& rtdbjson)
{
   std::string geomname = "geometry";
   if (rtdbjson.contains("geometry") && rtdbjson["geometry"].is_string())
      geomname = rtdbjson["geometry"].get<std::string>();
   return geomname;
}

/**************************************************
 *                                                *
 *               read_unita_3x3                   *
 *                                                *
 **************************************************/
/**
 * @brief Read a 3×3 lattice matrix from a flat JSON array.
 *
 * Deserializes a 3×3 lattice matrix (`unita`) from a JSON array of at least
 * nine numeric elements stored in row-major order:
 *
 *   [ a11, a12, a13,
 *     a21, a22, a23,
 *     a31, a32, a33 ]
 *
 * This function performs only structural validation and numeric conversion.
 * No assumptions are made about units, orthogonality, centering, or symmetry.
 *
 * @param[in]  j      JSON array containing the lattice matrix.
 * @param[out] unita  Output 3×3 lattice matrix (row-major).
 *
 * @return `true` if the JSON array contained at least nine numeric values
 *         and the matrix was successfully read; `false` otherwise.
 *
 * @note This helper is intentionally low-level and side-effect free.
 *       Higher-level logic (primitive conversion, symmetry resolution,
 *       restart policy) must be handled by the caller.
 */
static bool read_unita_3x3(const json& j, double unita[9])
{
   if (!j.is_array() || j.size() < 9) return false;
   for (int i=0;i<9;++i)
   {
      if (!(j[i].is_number_float() || j[i].is_number_integer())) return false;
      unita[i] = j[i].get<double>();
   }
   return true;
}

/**************************************************
 *                                                *
 *              write_unita_3x3                   *
 *                                                *
 **************************************************/
/**
 * @brief Write a 3×3 lattice matrix into JSON as a flat 9-element array.
 *
 * Serializes a 3×3 lattice matrix (unita) into a JSON array of length 9
 * using row-major ordering:
 *
 *   [ a11, a12, a13,
 *     a21, a22, a23,
 *     a31, a32, a33 ]
 *
 * This format is used consistently in the RTDB for simulation cells,
 * effective symmetry metadata, and restart consistency.
 *
 * @param[out] j     JSON array to receive the lattice matrix (overwritten)
 * @param[in]  unita Pointer to 9 doubles storing the lattice matrix
 *
 * @note No validation or unit conversion is performed.
 * @note Interpretation (Cartesian vs fractional basis) is handled elsewhere.
 */
static void write_unita_3x3(json& j, const double unita[9])
{
   j = json::array();
   for (int i=0;i<9;++i) j.push_back(unita[i]);
}

/**************************************************
 *                                                *
 *            read_symbols_coords                 *
 *                                                *
 **************************************************/
/**
 * @brief Read atomic symbols and Cartesian coordinates from a geometry JSON block.
 *
 * Extracts atomic symbols and Cartesian coordinates from a geometry object.
 * Coordinates may be stored either as a flat 3N array or as an Nx3 array.
 *
 * Supported coordinate formats:
 *   - Flat array: [x1,y1,z1,x2,y2,z2,...]
 *   - Matrix form: [[x1,y1,z1],[x2,y2,z2],...]
 *
 * Symbols are optional; if absent, only coordinates are read.
 * If coordinates are missing or malformed, coords_xyz will be empty.
 *
 * @param[in]  geomjson   Geometry JSON object
 * @param[out] symbols    Atomic symbols (cleared and refilled)
 * @param[out] coords_xyz Flat Cartesian coordinate array (3N)
 *
 * @throws std::runtime_error if flat coordinate array length is not divisible by 3
 *
 * @note No unit conversion is performed.
 * @note Fractional vs Cartesian interpretation is handled elsewhere.
 */
static void read_symbols_and_coords(const json& geomjson,
                                    std::vector<std::string>& symbols,
                                    std::vector<double>& coords_xyz)
{
   symbols.clear();
   coords_xyz.clear();

   if (geomjson.contains("symbols") && geomjson["symbols"].is_array())
      symbols = geomjson["symbols"].get<std::vector<std::string>>();

   if (!geomjson.contains("coords") || !geomjson["coords"].is_array())
      return;

   const json& c = geomjson["coords"];

   // Case 1: flat 3N array  ← YOUR CURRENT INPUT
   if (!c.empty() && c[0].is_number())
   {
      if (c.size() % 3 != 0)
         throw std::runtime_error("coords array length not divisible by 3");

      for (size_t i = 0; i < c.size(); ++i)
         coords_xyz.push_back(c[i].get<double>());

      return;
   }

   // Case 2: Nx3 array
   if (!c.empty() && c[0].is_array())
   {
      for (const auto& row : c)
      {
         if (row.size() < 3) continue;
         coords_xyz.push_back(row[0].get<double>());
         coords_xyz.push_back(row[1].get<double>());
         coords_xyz.push_back(row[2].get<double>());
      }
      return;
   }
}



/**************************************************
 *                                                *
 *              write_coords_xyz                  *
 *                                                *
 **************************************************/
/**
 * @brief Write Cartesian coordinates into RTDB JSON as an Nx3 array.
 *
 * Converts a flat 3N vector of Cartesian coordinates into a JSON array
 * of N rows, each containing three values `[x, y, z]`.
 *
 * This helper is used when recording a snapshot of the coordinates that
 * were used during symmetry resolution. The stored coordinates are
 * informational only (for restart consistency and debugging) and are
 * not interpreted as the active runtime geometry.
 *
 * @param[out] out_coords   JSON array to receive the Nx3 coordinates.
 * @param[in]  coords_xyz  Flat vector of length 3N containing Cartesian
 *                          coordinates ordered as (x₀,y₀,z₀,x₁,y₁,z₁,...).
 */
static void write_coords_xyz(json& out_coords, const std::vector<double>& coords_xyz)
{
   out_coords = json::array();
   const size_t n = coords_xyz.size()/3;
   for (size_t i=0;i<n;++i)
   {
      json row = json::array();
      row.push_back(coords_xyz[3*i+0]);
      row.push_back(coords_xyz[3*i+1]);
      row.push_back(coords_xyz[3*i+2]);
      out_coords.push_back(row);
   }
}


// ======================================================
// Main function
// ======================================================
/**************************************************
 *                                                *
 *          resolve_symmetry_and_cell             *
 *                                                *
 **************************************************/
/**
 * @brief Resolve symmetry intent and simulation cell into an authoritative
 *        effective symmetry description.
 *
 * This routine is the single entry point for converting user intent
 * (geometry symmetry directives, autosym flags, primitive requests)
 * into an immutable, self-consistent symmetry description stored in
 * `rtdb["effective_symmetry"]`.
 *
 * Responsibilities:
 *  - Determine the authoritative lattice (`unita`) with clear precedence:
 *      1. nwpw.simulation_cell.unita
 *      2. geometry.unita
 *  - Interpret symmetry intent from the geometry block (specified vs autosym).
 *  - Construct a concrete Symmetry object (identity, point group, or space group).
 *  - Optionally convert the lattice and coordinates to a primitive cell.
 *  - Freeze the resolved symmetry into `effective_symmetry` for downstream use.
 *  - Ensure the simulation cell seen by physics matches the effective symmetry.
 *  - Apply symmetry-based pruning of Monkhorst–Pack k-points.
 *
 * Design principles:
 *  - Symmetry is resolved exactly once and never re-interpreted downstream.
 *  - Restart safety: previously resolved symmetry may be reused if intent
 *    has not changed.
 *  - Operators are not stored in RTDB; symmetry is reconstructed deterministically
 *    from identifiers (name / number / fingerprint).
 *  - Downstream physics must treat `effective_symmetry` as authoritative.
 *
 * @param rtdbstring  Serialized RTDB JSON string.
 * @return Updated RTDB JSON string with resolved effective symmetry and
 *         symmetry-pruned Brillouin zone.
 *
 * @throws std::runtime_error if required geometry or lattice information
 *         is missing or inconsistent.
 */
std::string resolve_symmetry_and_cell(std::string rtdbstring)
{
   auto rtdbjson = json::parse(rtdbstring);

   const std::string geomname = get_geomname(rtdbjson);

   // ---- sanity: geometry block must exist
   if (!rtdbjson.contains("geometries") || !rtdbjson["geometries"].contains(geomname))
      throw std::runtime_error("resolve_symmetry_and_cell: missing geometries[geomname]");

   json& geomjson = rtdbjson["geometries"][geomname];

   // ---- precedence: unita comes from nwpw.simulation_cell first (if present), else from geometry
   double unita_in[9];
   bool have_unita = false;

   if (rtdbjson.contains("nwpw") &&
       rtdbjson["nwpw"].contains("simulation_cell") &&
       rtdbjson["nwpw"]["simulation_cell"].contains("unita") &&
       read_unita_3x3(rtdbjson["nwpw"]["simulation_cell"]["unita"], unita_in))
   {
      have_unita = true;
   }
   else if (geomjson.contains("unita") && read_unita_3x3(geomjson["unita"], unita_in))
   {
      have_unita = true;
   }

   if (!have_unita)
      throw std::runtime_error("resolve_symmetry_and_cell: could not find a valid unita");

   // ---- read coords + symbols from geometry (intent)
   std::vector<std::string> symbols;
   std::vector<double> coords_xyz;
   read_symbols_and_coords(geomjson, symbols, coords_xyz);

   // If coords missing, you may want to treat as error for crystals
   // (or allow empty for cell-only operations).
   // Here we allow empty but symmetry will likely degrade to identity.
   const bool have_coords = (coords_xyz.size() >= 3);

   // ---- pull symmetry intent from geometry
   bool symmetry_specified = false;
   std::string symmetry_group_name;
   int symmetry_group_number = 0;
   int symmetry_group_setting = 0;
   double symmetry_tolerance = 1.0e-6;
   bool symmetry_primitive_requested = false;
   bool primitive_suggested = false;

   if (geomjson.contains("symmetry") && geomjson["symmetry"].is_object())
   {
      const json& sj = geomjson["symmetry"];
      if (sj.contains("specified") && sj["specified"].is_boolean())
         symmetry_specified = sj["specified"].get<bool>();

      if (sj.contains("group_name") && sj["group_name"].is_string())
         symmetry_group_name = sj["group_name"].get<std::string>();

      if (sj.contains("group_number") && sj["group_number"].is_number_integer())
         symmetry_group_number = sj["group_number"].get<int>();

      if (sj.contains("setting") && sj["setting"].is_number_integer())
         symmetry_group_setting = sj["setting"].get<int>();

      if (sj.contains("tolerance") && (sj["tolerance"].is_number_float() || sj["tolerance"].is_number_integer()))
         symmetry_tolerance = sj["tolerance"].get<double>();

      if (sj.contains("primitive") && sj["primitive"].is_boolean())
         symmetry_primitive_requested = sj["primitive"].get<bool>();
   }

   // ---- autosym flag (you said you already store geomjson["autosym"])
   bool autosym   = read_bool_flag(geomjson, "autosym");
   bool autospace = read_bool_flag(geomjson, "autospace");
   if (autospace && !have_coords)
     autospace = false;  // degrade safely


   // ---- restart policy knobs (keep simple; you can refine)
   const bool have_effective = (rtdbjson.contains("effective_symmetry") && rtdbjson["effective_symmetry"].is_object());
   bool reuse_effective = false;

   // If you have a formal restart flag, wire it here.
   // For now: reuse only if effective exists AND user did not explicitly specify symmetry this run.
   // (i.e., specified symmetry should win over old effective_symmetry)
   if (have_effective && !symmetry_specified && rtdbjson["effective_symmetry"].value("primitive", false) == symmetry_primitive_requested)
      reuse_effective = true;


   // ---- if reusing, still ensure nwpw.simulation_cell.unita is consistent with effective
   if (reuse_effective)
   {
      json& es = rtdbjson["effective_symmetry"];

      if (es.contains("unita") && es["unita"].is_array())
      {
         // Force simulation_cell.unita to match effective
         rtdbjson["nwpw"]["simulation_cell"]["unita"] = es["unita"];
      }
      // NOTE: if you also store effective coords, you can push them into a working geometry here.

      return rtdbjson.dump();
   }

   // ======================================================
   // Build a new effective symmetry from intent
   // ======================================================
   pwdft::Symmetry sym;

   std::string sym_source = "identity";
   std::string sym_backend = "";

   if (symmetry_specified)
   {
       if (!symmetry_group_name.empty())
       {
           sym = pwdft::Symmetry(symmetry_group_name);
           sym_source = "specified";
       }
       else if (symmetry_group_number > 0)
       {
           sym = pwdft::Symmetry(symmetry_group_number);
           sym_source = "specified";
       }
       else
       {
        sym = pwdft::Symmetry(); // C1
       }
   }
   else if (autospace && have_coords)
   {
      AutoSpaceResult as = detect_space_group(unita_in, symbols, coords_xyz, symmetry_tolerance);
 
      if (as.success)
      {
          if (as.group_number > 0)
          {
              sym = pwdft::Symmetry(as.group_number);
              sym_source = "autospace";
              sym_backend = "spglib";
          }
          else if (!as.group_name.empty())
          {
              sym = pwdft::Symmetry(as.group_name);
              sym_source = "autospace";
              sym_backend = "spglib";
          }
          else
          {
              sym = pwdft::Symmetry(); // fallback safety
              sym_source = "autospace_failed";
          }
 
          // Primitive suggestion is advisory only
          if (as.primitive_suggested && !symmetry_primitive_requested)
          {
              // record suggestion, do NOT apply
              rtdbjson["effective_symmetry"]["primitive_suggested"] = true;
          }
      }
   }
   else if (autosym && have_coords)
   {
      // NOT IMPLEMENTED YET — placeholder
      // later: sym = detect_symmetry(...)
      sym = pwdft::Symmetry();
      sym_source = "autosym_failed";
      //sym_backend = "internal";
   }
   else
   {
      sym = pwdft::Symmetry(); // identity
      sym_source = "identity";
   }


   if (symmetry_primitive_requested && sym.is_space_group())
   {
      make_primitive_cell(sym, unita_in, coords_xyz);
   }



   // ======================================================
   // Write outputs into RTDB JSON (effective world + sim cell)
   // ======================================================
   json es;
   es["source"]    = sym_source;
   es["primitive"] = symmetry_primitive_requested;
   es["name"]      = sym.name();
   es["order"]     = sym.order();
   es["tolerance"] = symmetry_tolerance;
   es["coords_type"] = symmetry_primitive_requested ? "fractional" : "cartesian";
   es["translation_type"] = "fractional";    // how SymOp.t must be interpreted
   es["primitive_lattice_only"] = symmetry_primitive_requested;
   if (!sym_backend.empty()) 
      es["backend"] = sym_backend;


   //es["type"] = sym.is_space_group() ? "space_group" : "point_group";
   //es["type"] = sym.is_trivial() ? "trivial" : (sym.is_space_group() ? "space_group" : "point_group");
   if (sym.is_trivial())
      es["type"] = "trivial";
   else if (sym.is_space_group())
      es["type"] = "space_group";
   else
      es["type"] = "point_group";



   if (sym.is_space_group())
   {
      es["num_centering"] = sym.num_centering();

      if (symmetry_primitive_requested)
         es["num_centering"] = 1;
   }

   write_unita_3x3(es["unita"], unita_in);

   // write symmetry operators
   //es["ops"] = json::array();
   //for (const auto& op : sym.operators())
   //{
   //    json jop;
   //    jop["R"] = json::array();
   //    for (int i = 0; i < 3; ++i)
   //    {
   //       json row = json::array();
   //       for (int j = 0; j < 3; ++j)
   //          row.push_back(op.R[i][j]);
   //       jop["R"].push_back(row);
   //    }
//
   //    jop["t"] = { op.t[0], op.t[1], op.t[2] };
//
   //    es["ops"].push_back(jop);
   //}
   // Do NOT store full ops list in RTDB.
   // We can reconstruct from (name/group_number/setting) deterministically.
   es["store_ops"] = false;
   es["sym_fingerprint"] = symmetry_fingerprint(sym, unita_in, symmetry_primitive_requested, "fractional");





   // Snapshot of coordinates used to resolve symmetry (NOT runtime geometry)
   // These are cartesian coordinates taken from geometries[geomname] at symmetry resolution time.
   // They are stored ONLY for restart consistency and debugging.
   es["coords_xyz"] = json::array();
   write_coords_xyz(es["coords_xyz"], coords_xyz);


   // Commit effective symmetry
   rtdbjson["effective_symmetry"] = es;

   // Physics must see the effective lattice
   write_unita_3x3( rtdbjson["nwpw"]["simulation_cell"]["unita"], unita_in);



   // Symmetry pruning of kpoints
   bool prune_kpoints = true; // future option?

   // NOTE:
   // K-point symmetry pruning is applied after effective symmetry and
   // simulation_cell are finalized. Downstream physics must not re-interpret symmetry.

   if (prune_kpoints)
      return apply_kpoint_symmetry_pruning(rtdbjson.dump());
   else
      return rtdbjson.dump();

   //return rtdbjson.dump();
}


}
