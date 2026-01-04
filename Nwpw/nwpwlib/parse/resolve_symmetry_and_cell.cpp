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

#include "json.hpp"
using json = nlohmann::json;

#include "Symmetry.hpp"


namespace pwdft {


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
// ------------------------------
// helper: read 3x3 unita
// ------------------------------
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
// --------------------------------------
// helper: read coords (Nx3) and symbols
// --------------------------------------
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
// --------------------------------------
// helper: write coords (Nx3) flat 3N
// --------------------------------------
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
// PLACEHOLDER hooks (you implement / wire to Symmetry)
// ======================================================







// ======================================================
// Main function
// ======================================================
/**************************************************
 *                                                *
 *          resolve_symmetry_and_cell             *
 *                                                *
 **************************************************/
//void resolve_symmetry_and_cell(json& rtdbjson)
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
   bool autosym = false;
   if (geomjson.contains("autosym") && geomjson["autosym"].is_boolean())
      autosym = geomjson["autosym"].get<bool>();

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
   else if (autosym && have_coords)
   {
       // NOT IMPLEMENTED YET — placeholder
       // later: sym = detect_symmetry(...)
       sym = pwdft::Symmetry();
       sym_source = "autosym";
   }
   else
   {
       sym = pwdft::Symmetry(); // identity
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



   return rtdbjson.dump();
}


}
