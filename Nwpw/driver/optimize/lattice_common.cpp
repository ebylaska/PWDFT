// lattice_common.cpp
//
// Shared helpers used by the driver and by every per-system lattice
// minimizer. Nothing here is driver-specific; nothing here is
// crystal-system-specific except scale_cubic_cell, which will eventually
// be joined by scale_tetragonal_cell, scale_orthorhombic_cell, etc., or
// replaced by a generic set_lattice(...).
//

#include "lattice_common.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

namespace pwdft {

using json = nlohmann::json;

// ---------------------------------------------------------------------------
// print_lattice_state
//
// Diagnostics: dump the current cell, the frozen cell (if present), and the
// numerical grid descriptor from the RTDB. Called by compute_egs_values on
// every evaluation; lattice minimizers may also call it per step.
// ---------------------------------------------------------------------------

void print_lattice_state(const json& rtdbjson,
                         std::ostream& coutput,
                         const std::string& tag)
{
    const std::string geomname =
        (rtdbjson.contains("geometry") &&
         rtdbjson["geometry"].is_string())
            ? rtdbjson["geometry"].get<std::string>()
            : "geometry";

    if (!rtdbjson.contains("geometries") ||
        !rtdbjson["geometries"].is_object() ||
        !rtdbjson["geometries"].contains(geomname))
    {
        coutput << tag << "Geometry '" << geomname << "' is missing\n";
        return;
    }

    const json& geometry = rtdbjson["geometries"].at(geomname);

    const json* current_unita =
        (geometry.contains("unita") && geometry["unita"].is_array())
            ? &geometry["unita"]
            : nullptr;

    const json* frozen_unita = nullptr;

    if (rtdbjson.contains("nwpw") &&
        rtdbjson["nwpw"].is_object() &&
        rtdbjson["nwpw"].contains("simulation_cell") &&
        rtdbjson["nwpw"]["simulation_cell"].is_object())
    {
        const json& simulation_cell = rtdbjson["nwpw"]["simulation_cell"];

        if (simulation_cell.contains("unita_frozen") &&
            simulation_cell["unita_frozen"].is_array())
        {
            frozen_unita = &simulation_cell["unita_frozen"];
        }
    }

    if (current_unita != nullptr && current_unita->size() == 9)
    {
        coutput << tag << "Current unita:\n";
        for (int j = 0; j < 3; ++j)
        {
            coutput << std::setprecision(10)
                    << tag << "  "
                    << (*current_unita)[3*j + 0].get<double>() << " "
                    << (*current_unita)[3*j + 1].get<double>() << " "
                    << (*current_unita)[3*j + 2].get<double>() << '\n';
        }
    }

    if (frozen_unita != nullptr && frozen_unita->size() == 9)
    {
        coutput << tag << "Frozen unita_frozen:\n";
        for (int j = 0; j < 3; ++j)
        {
            coutput << std::setprecision(10)
                    << tag << "  "
                    << (*frozen_unita)[3*j + 0].get<double>() << " "
                    << (*frozen_unita)[3*j + 1].get<double>() << " "
                    << (*frozen_unita)[3*j + 2].get<double>() << '\n';
        }
    }

    if (rtdbjson.contains("driver") &&
        rtdbjson["driver"].is_object() &&
        rtdbjson["driver"].contains("numerical_grid") &&
        rtdbjson["driver"]["numerical_grid"].is_object())
    {
        const json& grid = rtdbjson["driver"]["numerical_grid"];

        coutput << tag
                << "FFT grid = " << grid.value("nx", -1)
                << " x "         << grid.value("ny", -1)
                << " x "         << grid.value("nz", -1)
                << " waves0 = "  << grid.value("nwave0", -1)
                << " waves1 = "  << grid.value("nwave1", -1)
                << " npack0 = "  << grid.value("npack0", -1)
                << " npack1 = "  << grid.value("npack1", -1)
                << '\n';
    }
    else
    {
        coutput << tag << "Numerical grid: unavailable\n";
    }
}

// ---------------------------------------------------------------------------
// read_unita / write_unita
//
// Convert between a JSON array of 9 doubles (row-major 3x3 direct lattice)
// and std::array<double,9>. read_unita returns false on any malformed input;
// callers decide whether that is fatal.
// ---------------------------------------------------------------------------

bool read_unita(const json& value, std::array<double, 9>& unita)
{
    if (!value.is_array() || value.size() != 9)
        return false;

    for (int i = 0; i < 9; ++i)
    {
        if (!value[i].is_number())
            return false;
        unita[i] = value[i].get<double>();
    }

    return true;
}

void write_unita(json& value, const std::array<double, 9>& unita)
{
    value = json::array();
    for (double x : unita)
        value.push_back(x);
}

// ---------------------------------------------------------------------------
// unita_relative_difference
//
// Frobenius-norm relative difference between two 3x3 lattices. Denominator
// is floored at 1 to avoid a blow-up for a near-zero reference cell.
// ---------------------------------------------------------------------------

double unita_relative_difference(const std::array<double, 9>& current,
                                 const std::array<double, 9>& frozen)
{
    double difference_squared = 0.0;
    double frozen_squared     = 0.0;

    for (int i = 0; i < 9; ++i)
    {
        const double difference = current[i] - frozen[i];
        difference_squared += difference * difference;
        frozen_squared     += frozen[i] * frozen[i];
    }

    const double denominator = std::max(1.0, std::sqrt(frozen_squared));
    return std::sqrt(difference_squared) / denominator;
}

// ---------------------------------------------------------------------------
// scale_cubic_cell
//
// Isotropic scaling of a cubic cell by `scale`. Updates geometry.unita,
// geometry.coords, and (if present) nwpw.simulation_cell.unita, then writes
// the modified RTDB back into the caller's string.
//
// NOTE: this name is a lie for non-cubic cells — it scales all nine lattice
// components uniformly, which is only physically correct when a = b = c and
// all angles are 90 degrees. When you add tetragonal/orthorhombic/etc.
// minimizers, introduce a generic set_lattice(...) alongside this and have
// scale_cubic_cell forward to it. Do not call it from a non-cubic minimizer.
// ---------------------------------------------------------------------------

void scale_cubic_cell(std::string& rtdbstring, const double scale)
{
    json rtdbjson = json::parse(rtdbstring);

    const std::string geomname =
        (rtdbjson.contains("geometry") &&
         rtdbjson["geometry"].is_string())
            ? rtdbjson["geometry"].get<std::string>()
            : "geometry";

    json& geometry = rtdbjson["geometries"][geomname];

    for (int i = 0; i < 9; ++i)
    {
        geometry["unita"][i] =
            geometry["unita"][i].get<double>() * scale;
    }

    if (geometry.contains("coords") && geometry["coords"].is_array())
    {
        for (std::size_t i = 0; i < geometry["coords"].size(); ++i)
        {
            geometry["coords"][i] =
                geometry["coords"][i].get<double>() * scale;
        }
    }

    if (rtdbjson.contains("nwpw") &&
        rtdbjson["nwpw"].is_object() &&
        rtdbjson["nwpw"].contains("simulation_cell") &&
        rtdbjson["nwpw"]["simulation_cell"].is_object())
    {
        json& cell = rtdbjson["nwpw"]["simulation_cell"];

        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            for (int i = 0; i < 9; ++i)
            {
                cell["unita"][i] =
                    cell["unita"][i].get<double>() * scale;
            }
        }
    }

    rtdbstring = rtdbjson.dump();
}

// ---------------------------------------------------------------------------
// compute_egs_values
//
// Prepare an RTDB request for the selected calculation, invoke the electronic
// minimizer, extract the requested quantities from the backend section, and
// return them as a JSON object with the keys:
//
//     energy, gradient, stress, stress_sym, lstress
//
// Missing backend fields come back as null JSON values.
//
// option:
//     1 -> energy    (current_task = "task pspw energy")
//     2 -> gradient  (current_task = "task pspw gradient")
//     3 -> stress    (current_task = "task pspw stress", includestress = true)
//
// rtdbstring is replaced with the RTDB returned by the minimizer, so callers
// that want a fixed reference lattice across several evaluations should keep
// their own copy of the input string and pass fresh copies per call.
//
// Throws std::invalid_argument for a null minimizer or an out-of-range option,
// std::runtime_error if the minimizer reports failure or no backend section is
// present, and nlohmann::json::exception if the returned RTDB is malformed.
// ---------------------------------------------------------------------------

json compute_egs_values(const int option,
                        MPI_Comm comm_world0,
                        electronic_minimizer minimizer,
                        std::string& rtdbstring,
                        std::ostream& coutput)
{
    if (minimizer == nullptr)
        throw std::invalid_argument("compute_egs_values: null minimizer");
    if (option < 1 || option > 3)
        throw std::invalid_argument("compute_egs_values: invalid option");

    json request = json::parse(rtdbstring);

    if (option == 1)
    {
        request["current_task"] = "task pspw energy";
        request["nwpw"]["includestress"] = false;
    }
    else if (option == 2)
    {
        request["current_task"] = "task pspw gradient";
        request["nwpw"]["includestress"] = false;
    }
    else // option == 3
    {
        request["current_task"] = "task pspw stress";
        request["nwpw"]["includestress"] = true;
    }

    request["driver"]["cell_optimization"] = true;
    request["driver"]["use_frozen_lattice"] = true;

    const std::string tag = "@";
    print_lattice_state(request, coutput, tag);

    std::string rtdbstring1 = request.dump();

    const int ierr = minimizer(comm_world0, rtdbstring1, coutput);

    if (ierr != 0)
        throw std::runtime_error(
            "compute_egs_values: minimizer failed with error " +
            std::to_string(ierr));

    const json output = json::parse(rtdbstring1);

    const json* backend = nullptr;
    if (output.contains("pspw") && output["pspw"].is_object())
        backend = &output["pspw"];
    else if (output.contains("band") && output["band"].is_object())
        backend = &output["band"];
    else
        throw std::runtime_error(
            "compute_egs_values: neither pspw nor band results found");

    const json& backend_result = *backend;

    json result = json::object();
    result["energy"]     = backend_result.value("energy",     json{});
    result["gradient"]   = backend_result.value("gradient",   json{});
    result["stress"]     = backend_result.value("stress",     json{});
    result["stress_sym"] = backend_result.value("stress_sym", json{});
    result["lstress"]    = backend_result.value("lstress",    json{});

    rtdbstring = std::move(rtdbstring1);

    return result;
}


// lattice_common.cpp (additions)

/*************************************
 *                                   *
 *        read_tetragonal_lattice    *
 *                                   *
 *************************************/

std::pair<double, double> read_tetragonal_lattice(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    const auto& unita = rtdb.at("geometries").at(geomname).at("unita");

    const double a = unita.at(0).get<double>();
    const double c = unita.at(8).get<double>();

    return { a, c };
}


/*************************************
 *                                   *
 *        set_tetragonal_cell        *
 *                                   *
 *************************************/
void set_tetragonal_cell(std::string& rtdbstring, double a_new, double c_new)
{
    if (!(a_new > 0.0) || !(c_new > 0.0))
        throw std::runtime_error("set_tetragonal_cell: non-positive lattice parameter");

    json rtdb = json::parse(rtdbstring);

    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    json& geometry = rtdb["geometries"][geomname];

    const double a_old = geometry["unita"].at(0).get<double>();
    const double c_old = geometry["unita"].at(8).get<double>();

    if (!(a_old > 0.0) || !(c_old > 0.0))
        throw std::runtime_error("set_tetragonal_cell: invalid current lattice");

    const double scale_a = a_new / a_old;
    const double scale_c = c_new / c_old;

    // Diagonal-only update. Assumes an axis-aligned tetragonal cell:
    //   unita = [a, 0, 0, 0, a, 0, 0, 0, c]
    geometry["unita"][0] = a_new;
    geometry["unita"][4] = a_new;
    geometry["unita"][8] = c_new;

    // Cartesian coords: x,y scale with a; z scales with c
    if (geometry.contains("coords") && geometry["coords"].is_array())
    {
        auto& coords = geometry["coords"];
        for (std::size_t i = 0; i + 2 < coords.size(); i += 3)
        {
            coords[i + 0] = coords[i + 0].get<double>() * scale_a;
            coords[i + 1] = coords[i + 1].get<double>() * scale_a;
            coords[i + 2] = coords[i + 2].get<double>() * scale_c;
        }
    }

    // Keep nwpw.simulation_cell.unita in sync if present
    if (rtdb.contains("nwpw") &&
        rtdb["nwpw"].is_object() &&
        rtdb["nwpw"].contains("simulation_cell") &&
        rtdb["nwpw"]["simulation_cell"].is_object())
    {
        json& cell = rtdb["nwpw"]["simulation_cell"];

        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            cell["unita"][0] = cell["unita"][0].get<double>() * scale_a;
            cell["unita"][4] = cell["unita"][4].get<double>() * scale_a;
            cell["unita"][8] = cell["unita"][8].get<double>() * scale_c;
        }
    }

    rtdbstring = rtdb.dump();
}


/*************************************
 *                                   *
 *       read_hexagonal_lattice      *
 *                                   *
 *************************************/

std::pair<double, double> read_hexagonal_lattice(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    const auto& unita = rtdb.at("geometries").at(geomname).at("unita");

    const double u0 = unita.at(0).get<double>();
    const double u1 = unita.at(1).get<double>();
    const double u2 = unita.at(2).get<double>();
    const double a  = std::sqrt(u0*u0 + u1*u1 + u2*u2);

    const double u6 = unita.at(6).get<double>();
    const double u7 = unita.at(7).get<double>();
    const double u8 = unita.at(8).get<double>();
    const double c  = std::sqrt(u6*u6 + u7*u7 + u8*u8);

    return { a, c };
}

/*************************************
 *                                   *
 *        set_hexagonal_cell         *
 *                                   *
 *************************************/
void set_hexagonal_cell(std::string& rtdbstring, double a_new, double c_new)
{
    if (!(a_new > 0.0) || !(c_new > 0.0))
        throw std::runtime_error(
            "set_hexagonal_cell: non-positive lattice parameter");

    json rtdb = json::parse(rtdbstring);

    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    json& geometry = rtdb["geometries"][geomname];

    const double u0 = geometry["unita"].at(0).get<double>();
    const double u1 = geometry["unita"].at(1).get<double>();
    const double u2 = geometry["unita"].at(2).get<double>();
    const double u6 = geometry["unita"].at(6).get<double>();
    const double u7 = geometry["unita"].at(7).get<double>();
    const double u8 = geometry["unita"].at(8).get<double>();

    const double a_old = std::sqrt(u0*u0 + u1*u1 + u2*u2);
    const double c_old = std::sqrt(u6*u6 + u7*u7 + u8*u8);

    if (!(a_old > 0.0) || !(c_old > 0.0))
        throw std::runtime_error(
            "set_hexagonal_cell: invalid current lattice");

    const double scale_a = a_new / a_old;
    const double scale_c = c_new / c_old;

    // Scale a1 and a2 by scale_a; a3 by scale_c.
    // Assumes the standard hexagonal orientation: a1 along x, a2 in the
    // xy-plane at 120 degrees, a3 along z.
    for (int i = 0; i < 3; ++i)
        geometry["unita"][i] =
            geometry["unita"][i].get<double>() * scale_a;
    for (int i = 3; i < 6; ++i)
        geometry["unita"][i] =
            geometry["unita"][i].get<double>() * scale_a;
    for (int i = 6; i < 9; ++i)
        geometry["unita"][i] =
            geometry["unita"][i].get<double>() * scale_c;

    // Cartesian coords: in-plane (x, y) by scale_a; z by scale_c.
    if (geometry.contains("coords") && geometry["coords"].is_array())
    {
        auto& coords = geometry["coords"];
        for (std::size_t i = 0; i + 2 < coords.size(); i += 3)
        {
            coords[i + 0] = coords[i + 0].get<double>() * scale_a;
            coords[i + 1] = coords[i + 1].get<double>() * scale_a;
            coords[i + 2] = coords[i + 2].get<double>() * scale_c;
        }
    }

    // Keep nwpw.simulation_cell.unita in sync if present.
    if (rtdb.contains("nwpw") &&
        rtdb["nwpw"].is_object() &&
        rtdb["nwpw"].contains("simulation_cell") &&
        rtdb["nwpw"]["simulation_cell"].is_object())
    {
        json& cell = rtdb["nwpw"]["simulation_cell"];

        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            for (int i = 0; i < 3; ++i)
                cell["unita"][i] =
                    cell["unita"][i].get<double>() * scale_a;
            for (int i = 3; i < 6; ++i)
                cell["unita"][i] =
                    cell["unita"][i].get<double>() * scale_a;
            for (int i = 6; i < 9; ++i)
                cell["unita"][i] =
                    cell["unita"][i].get<double>() * scale_c;
        }
    }

    rtdbstring = rtdb.dump();
}

std::array<double, 3> read_orthorhombic_lattice(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    const auto& unita = rtdb.at("geometries").at(geomname).at("unita");

    // Diagonal orthorhombic cell: lengths on the diagonal.
    const double a = unita.at(0).get<double>();
    const double b = unita.at(4).get<double>();
    const double c = unita.at(8).get<double>();

    return { a, b, c };
}

void set_orthorhombic_cell(std::string& rtdbstring, double a_new, double b_new, double c_new)
{
    if (!(a_new > 0.0) || !(b_new > 0.0) || !(c_new > 0.0))
        throw std::runtime_error(
            "set_orthorhombic_cell: non-positive lattice parameter");

    json rtdb = json::parse(rtdbstring);

    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    json& geometry = rtdb["geometries"][geomname];

    const double a_old = geometry["unita"].at(0).get<double>();
    const double b_old = geometry["unita"].at(4).get<double>();
    const double c_old = geometry["unita"].at(8).get<double>();

    if (!(a_old > 0.0) || !(b_old > 0.0) || !(c_old > 0.0))
        throw std::runtime_error(
            "set_orthorhombic_cell: invalid current lattice");

    const double scale_a = a_new / a_old;
    const double scale_b = b_new / b_old;
    const double scale_c = c_new / c_old;

    // Row-major 3x3. Orthorhombic cells are diagonal:
    //   a1 = (a, 0, 0), a2 = (0, b, 0), a3 = (0, 0, c)
    for (int i = 0; i < 3; ++i)
        geometry["unita"][i] =
            geometry["unita"][i].get<double>() * scale_a;
    for (int i = 3; i < 6; ++i)
        geometry["unita"][i] =
            geometry["unita"][i].get<double>() * scale_b;
    for (int i = 6; i < 9; ++i)
        geometry["unita"][i] =
            geometry["unita"][i].get<double>() * scale_c;

    // Cartesian coords: x by scale_a, y by scale_b, z by scale_c.
    if (geometry.contains("coords") && geometry["coords"].is_array())
    {
        auto& coords = geometry["coords"];
        for (std::size_t i = 0; i + 2 < coords.size(); i += 3)
        {
            coords[i + 0] = coords[i + 0].get<double>() * scale_a;
            coords[i + 1] = coords[i + 1].get<double>() * scale_b;
            coords[i + 2] = coords[i + 2].get<double>() * scale_c;
        }
    }

    // Keep nwpw.simulation_cell.unita in sync if present.
    if (rtdb.contains("nwpw") &&
        rtdb["nwpw"].is_object() &&
        rtdb["nwpw"].contains("simulation_cell") &&
        rtdb["nwpw"]["simulation_cell"].is_object())
    {
        json& cell = rtdb["nwpw"]["simulation_cell"];

        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            for (int i = 0; i < 3; ++i)
                cell["unita"][i] =
                    cell["unita"][i].get<double>() * scale_a;
            for (int i = 3; i < 6; ++i)
                cell["unita"][i] =
                    cell["unita"][i].get<double>() * scale_b;
            for (int i = 6; i < 9; ++i)
                cell["unita"][i] =
                    cell["unita"][i].get<double>() * scale_c;
        }
    }

    rtdbstring = rtdb.dump();
}


MonoclinicLattice read_monoclinic_lattice(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    const auto& unita = rtdb.at("geometries").at(geomname).at("unita");

    // Canonical b-unique monoclinic setting:
    //   a1 = (a,  0,  0)
    //   a2 = (0,  b,  0)
    //   a3 = (c cosβ, 0, c sinβ)
    const double a = unita.at(0).get<double>();
    const double b = unita.at(4).get<double>();
    const double u6 = unita.at(6).get<double>();
    const double u8 = unita.at(8).get<double>();
    const double c  = std::hypot(u6, u8);
    const double beta = std::atan2(u8, u6);

    return { a, b, c, beta };
}

void set_monoclinic_cell(std::string& rtdbstring,
                    double a_new, double b_new, double c_new,
                    double beta_new_rad)
{
    if (!(a_new > 0.0) || !(b_new > 0.0) || !(c_new > 0.0))
        throw std::runtime_error("set_monoclinic_cell: non-positive length");
    if (!(beta_new_rad > 0.0) || !(beta_new_rad < M_PI))
        throw std::runtime_error("set_monoclinic_cell: beta out of range");

    json rtdb = json::parse(rtdbstring);

    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    json& geometry = rtdb["geometries"][geomname];

    // --- read old cell ---
    const double a_old = geometry["unita"].at(0).get<double>();
    const double b_old = geometry["unita"].at(4).get<double>();
    const double u6_old = geometry["unita"].at(6).get<double>();
    const double u8_old = geometry["unita"].at(8).get<double>();
    const double c_old  = std::hypot(u6_old, u8_old);
    const double beta_old = std::atan2(u8_old, u6_old);

    if (!(a_old > 0.0) || !(b_old > 0.0) || !(c_old > 0.0))
        throw std::runtime_error("set_monoclinic_cell: invalid current cell");

    const double cb_old = std::cos(beta_old);
    const double sb_old = std::sin(beta_old);
    const double cb_new = std::cos(beta_new_rad);
    const double sb_new = std::sin(beta_new_rad);

    // --- update unita to the new cell ---
    geometry["unita"][0] = a_new;
    geometry["unita"][1] = 0.0;
    geometry["unita"][2] = 0.0;
    geometry["unita"][3] = 0.0;
    geometry["unita"][4] = b_new;
    geometry["unita"][5] = 0.0;
    geometry["unita"][6] = c_new * cb_new;
    geometry["unita"][7] = 0.0;
    geometry["unita"][8] = c_new * sb_new;

    // --- update coords, preserving fractional positions ---
    //
    // Old cartesian (x, y, z) -> fractional (f1, f2, f3):
    //   f2 = y / b_old
    //   f3 = z / (c_old * sb_old)
    //   f1 = (x - f3 * c_old * cb_old) / a_old
    //
    // Fractional -> new cartesian:
    //   x' = f1 * a_new + f3 * c_new * cb_new
    //   y' = f2 * b_new
    //   z' = f3 * c_new * sb_new

    if (geometry.contains("coords") && geometry["coords"].is_array())
    {
        auto& coords = geometry["coords"];
        for (std::size_t i = 0; i + 2 < coords.size(); i += 3)
        {
            const double x = coords[i + 0].get<double>();
            const double y = coords[i + 1].get<double>();
            const double z = coords[i + 2].get<double>();

            const double f2 = y / b_old;
            const double f3 = z / (c_old * sb_old);
            const double f1 = (x - f3 * c_old * cb_old) / a_old;

            coords[i + 0] = f1 * a_new + f3 * c_new * cb_new;
            coords[i + 1] = f2 * b_new;
            coords[i + 2] = f3 * c_new * sb_new;
        }
    }

    // --- keep nwpw.simulation_cell.unita in sync ---
    if (rtdb.contains("nwpw") &&
        rtdb["nwpw"].is_object() &&
        rtdb["nwpw"].contains("simulation_cell") &&
        rtdb["nwpw"]["simulation_cell"].is_object())
    {
        json& cell = rtdb["nwpw"]["simulation_cell"];

        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            // Read the old values from here too, then apply the same update.
            const double ao = cell["unita"].at(0).get<double>();
            const double bo = cell["unita"].at(4).get<double>();
            const double u6o = cell["unita"].at(6).get<double>();
            const double u8o = cell["unita"].at(8).get<double>();
            const double co = std::hypot(u6o, u8o);
            const double bto = std::atan2(u8o, u6o);

            // For simplicity: if the simulation_cell had the same cell as
            // geometry, just overwrite. If not, this is a corner case we
            // don't handle; caller should not mix cells.
            (void)ao; (void)bo; (void)co; (void)bto;

            cell["unita"][0] = a_new;
            cell["unita"][1] = 0.0;
            cell["unita"][2] = 0.0;
            cell["unita"][3] = 0.0;
            cell["unita"][4] = b_new;
            cell["unita"][5] = 0.0;
            cell["unita"][6] = c_new * cb_new;
            cell["unita"][7] = 0.0;
            cell["unita"][8] = c_new * sb_new;
        }
    }

    rtdbstring = rtdb.dump();
}

RhombohedralLattice read_rhombohedral_lattice(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    const auto& unita = rtdb.at("geometries").at(geomname).at("unita");

    auto row = [&](int i) {
        return std::array<double, 3>{
            unita.at(3*i + 0).get<double>(),
            unita.at(3*i + 1).get<double>(),
            unita.at(3*i + 2).get<double>()
        };
    };

    const auto a1 = row(0);
    const auto a2 = row(1);
    const auto a3 = row(2);

    auto norm = [](const std::array<double, 3>& v) {
        return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    };
    auto dot = [](const std::array<double, 3>& u,
                  const std::array<double, 3>& v) {
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    };

    const double l1 = norm(a1);
    const double l2 = norm(a2);
    const double l3 = norm(a3);
    const double a  = (l1 + l2 + l3) / 3.0;

    // Average all three pairwise angles; in a well-formed R cell they're equal.
    const double cos12 = dot(a1, a2) / (l1 * l2);
    const double cos13 = dot(a1, a3) / (l1 * l3);
    const double cos23 = dot(a2, a3) / (l2 * l3);
    const double cos_alpha = (cos12 + cos13 + cos23) / 3.0;

    double alpha = std::acos(std::clamp(cos_alpha, -1.0, 1.0));

    return { a, alpha };
}

void set_rhombohedral_cell(std::string& rtdbstring, double a_new, double alpha_new_rad)
{
    if (!(a_new > 0.0))
        throw std::runtime_error("set_rhombohedral_cell: non-positive a");
    if (!(alpha_new_rad > 0.0) || !(alpha_new_rad < M_PI))
        throw std::runtime_error("set_rhombohedral_cell: alpha out of range");

    json rtdb = json::parse(rtdbstring);

    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    json& geometry = rtdb["geometries"][geomname];

    // --- read old cell ---
    auto row_old = [&](int i) {
        return std::array<double, 3>{
            geometry["unita"][3*i + 0].get<double>(),
            geometry["unita"][3*i + 1].get<double>(),
            geometry["unita"][3*i + 2].get<double>()
        };
    };
    const auto a1_old = row_old(0);
    const auto a2_old = row_old(1);
    const auto a3_old = row_old(2);

    // --- build new cell ---
    const double r = 2.0 * std::sin(alpha_new_rad / 2.0) / std::sqrt(3.0);
    const double z2 = 1.0 - r*r;
    if (z2 <= 0.0)
        throw std::runtime_error(
            "set_rhombohedral_cell: alpha too small, z^2 negative");
    const double z = std::sqrt(z2);

    const double rt3_2 = r * std::sqrt(3.0) / 2.0;

    auto row_new = [&](int i) {
        switch (i) {
            case 0: return std::array<double, 3>{  a_new * r,       0.0,        a_new * z };
            case 1: return std::array<double, 3>{ -a_new * r / 2.0, a_new * rt3_2, a_new * z };
            default:return std::array<double, 3>{ -a_new * r / 2.0,-a_new * rt3_2, a_new * z };
        }
    };
    const auto a1_new = row_new(0);
    const auto a2_new = row_new(1);
    const auto a3_new = row_new(2);

    // --- invert old cell (rows a1, a2, a3) ---
    const double m00 = a1_old[0], m01 = a1_old[1], m02 = a1_old[2];
    const double m10 = a2_old[0], m11 = a2_old[1], m12 = a2_old[2];
    const double m20 = a3_old[0], m21 = a3_old[1], m22 = a3_old[2];

    const double det = m00*(m11*m22 - m12*m21)
                     - m01*(m10*m22 - m12*m20)
                     + m02*(m10*m21 - m11*m20);
    if (std::abs(det) < 1.0e-14)
        throw std::runtime_error("set_rhombohedral_cell: singular old cell");

    const double inv00 =  (m11*m22 - m12*m21) / det;
    const double inv01 = -(m01*m22 - m02*m21) / det;
    const double inv02 =  (m01*m12 - m02*m11) / det;
    const double inv10 = -(m10*m22 - m12*m20) / det;
    const double inv11 =  (m00*m22 - m02*m20) / det;
    const double inv12 = -(m00*m12 - m02*m10) / det;
    const double inv20 =  (m10*m21 - m11*m20) / det;
    const double inv21 = -(m00*m21 - m01*m20) / det;
    const double inv22 =  (m00*m11 - m01*m10) / det;

    // --- update unita ---
    geometry["unita"][0] = a1_new[0]; geometry["unita"][1] = a1_new[1]; geometry["unita"][2] = a1_new[2];
    geometry["unita"][3] = a2_new[0]; geometry["unita"][4] = a2_new[1]; geometry["unita"][5] = a2_new[2];
    geometry["unita"][6] = a3_new[0]; geometry["unita"][7] = a3_new[1]; geometry["unita"][8] = a3_new[2];

    // --- rescale coords, preserving fractional positions ---
    if (geometry.contains("coords") && geometry["coords"].is_array())
    {
        auto& coords = geometry["coords"];
        for (std::size_t i = 0; i + 2 < coords.size(); i += 3)
        {
            const double x = coords[i + 0].get<double>();
            const double y = coords[i + 1].get<double>();
            const double zc = coords[i + 2].get<double>();

            // frac = cart . A^{-1}
            const double f1 = x*inv00 + y*inv10 + zc*inv20;
            const double f2 = x*inv01 + y*inv11 + zc*inv21;
            const double f3 = x*inv02 + y*inv12 + zc*inv22;

            // cart_new = frac . A_new
            coords[i + 0] = f1*a1_new[0] + f2*a2_new[0] + f3*a3_new[0];
            coords[i + 1] = f1*a1_new[1] + f2*a2_new[1] + f3*a3_new[1];
            coords[i + 2] = f1*a1_new[2] + f2*a2_new[2] + f3*a3_new[2];
        }
    }

    // --- sync nwpw.simulation_cell.unita ---
    if (rtdb.contains("nwpw") &&
        rtdb["nwpw"].is_object() &&
        rtdb["nwpw"].contains("simulation_cell") &&
        rtdb["nwpw"]["simulation_cell"].is_object())
    {
        json& cell = rtdb["nwpw"]["simulation_cell"];
        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            cell["unita"][0] = a1_new[0]; cell["unita"][1] = a1_new[1]; cell["unita"][2] = a1_new[2];
            cell["unita"][3] = a2_new[0]; cell["unita"][4] = a2_new[1]; cell["unita"][5] = a2_new[2];
            cell["unita"][6] = a3_new[0]; cell["unita"][7] = a3_new[1]; cell["unita"][8] = a3_new[2];
        }
    }

    rtdbstring = rtdb.dump();
}

bool
invert_3x3(const std::array<double, 9>& M, std::array<double, 9>& Minv)
{
    const double m00 = M[0], m01 = M[1], m02 = M[2];
    const double m10 = M[3], m11 = M[4], m12 = M[5];
    const double m20 = M[6], m21 = M[7], m22 = M[8];

    const double det = m00*(m11*m22 - m12*m21)
                     - m01*(m10*m22 - m12*m20)
                     + m02*(m10*m21 - m11*m20);
    if (std::abs(det) < 1.0e-14)
        return false;

    const double id = 1.0 / det;
    Minv[0] =  (m11*m22 - m12*m21) * id;
    Minv[1] = -(m01*m22 - m02*m21) * id;
    Minv[2] =  (m01*m12 - m02*m11) * id;
    Minv[3] = -(m10*m22 - m12*m20) * id;
    Minv[4] =  (m00*m22 - m02*m20) * id;
    Minv[5] = -(m00*m12 - m02*m10) * id;
    Minv[6] =  (m10*m21 - m11*m20) * id;
    Minv[7] = -(m00*m21 - m01*m20) * id;
    Minv[8] =  (m00*m11 - m01*m10) * id;
    return true;
}

TriclinicLattice
read_triclinic_lattice(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    const auto& unita = rtdb.at("geometries").at(geomname).at("unita");

    auto row = [&](int i) {
        return std::array<double, 3>{
            unita.at(3*i + 0).get<double>(),
            unita.at(3*i + 1).get<double>(),
            unita.at(3*i + 2).get<double>()
        };
    };
    const auto a1 = row(0);
    const auto a2 = row(1);
    const auto a3 = row(2);

    auto norm = [](const std::array<double, 3>& v) {
        return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    };
    auto dot = [](const std::array<double, 3>& u,
                  const std::array<double, 3>& v) {
        return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
    };
    auto safe_acos = [](double c) {
        return std::acos(std::clamp(c, -1.0, 1.0));
    };

    const double a = norm(a1);
    const double b = norm(a2);
    const double c = norm(a3);

    const double alpha = safe_acos(dot(a2, a3) / (b * c));
    const double beta  = safe_acos(dot(a1, a3) / (a * c));
    const double gamma = safe_acos(dot(a1, a2) / (a * b));

    return { a, b, c, alpha, beta, gamma };
}

void
set_triclinic_cell(std::string& rtdbstring,
                   double a, double b, double c,
                   double alpha_rad, double beta_rad, double gamma_rad)
{
    if (!(a > 0.0) || !(b > 0.0) || !(c > 0.0))
        throw std::runtime_error("set_triclinic_cell: non-positive length");
    if (!(alpha_rad > 0.0) || !(alpha_rad < M_PI) ||
        !(beta_rad  > 0.0) || !(beta_rad  < M_PI) ||
        !(gamma_rad > 0.0) || !(gamma_rad < M_PI))
        throw std::runtime_error("set_triclinic_cell: angle out of range");

    const double ca = std::cos(alpha_rad);
    const double cb = std::cos(beta_rad);
    const double cg = std::cos(gamma_rad);
    const double sg = std::sin(gamma_rad);

    // Standard triclinic cell in the "a1 along x, a2 in xy-plane" convention:
    //   a1 = (a, 0, 0)
    //   a2 = (b cg, b sg, 0)
    //   a3 = (c cb, c (ca - cb cg)/sg, c·V/sg)
    // where V = sqrt(1 - ca² - cb² - cg² + 2 ca cb cg)
    const double V2 = 1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg;
    if (V2 <= 0.0)
        throw std::runtime_error("set_triclinic_cell: degenerate angles");
    const double V = std::sqrt(V2);

    const double a3x = c * cb;
    const double a3y = c * (ca - cb * cg) / sg;
    const double a3z = c * V / sg;

    const std::array<double, 9> A_new{
        a, 0.0, 0.0,
        b * cg, b * sg, 0.0,
        a3x, a3y, a3z
    };

    json rtdb = json::parse(rtdbstring);

    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";

    json& geometry = rtdb["geometries"][geomname];

    // Read the old cell matrix directly from unita; no need to reconstruct
    // from parameters, and this is robust to whatever orientation the
    // parser emitted.
    std::array<double, 9> A_old{};
    for (int i = 0; i < 9; ++i)
        A_old[i] = geometry["unita"][i].get<double>();

    std::array<double, 9> A_inv{};
    if (!invert_3x3(A_old, A_inv))
        throw std::runtime_error("set_triclinic_cell: singular old cell");

    // Replace unita with the new cell matrix.
    for (int i = 0; i < 9; ++i)
        geometry["unita"][i] = A_new[i];

    // Rescale cartesian coords, preserving fractional positions.
    if (geometry.contains("coords") && geometry["coords"].is_array())
    {
        auto& coords = geometry["coords"];
        for (std::size_t i = 0; i + 2 < coords.size(); i += 3)
        {
            const double x = coords[i + 0].get<double>();
            const double y = coords[i + 1].get<double>();
            const double z = coords[i + 2].get<double>();

            // frac = cart · A^{-1}
            const double f1 = x*A_inv[0] + y*A_inv[3] + z*A_inv[6];
            const double f2 = x*A_inv[1] + y*A_inv[4] + z*A_inv[7];
            const double f3 = x*A_inv[2] + y*A_inv[5] + z*A_inv[8];

            // cart_new = frac · A_new
            coords[i + 0] = f1*A_new[0] + f2*A_new[3] + f3*A_new[6];
            coords[i + 1] = f1*A_new[1] + f2*A_new[4] + f3*A_new[7];
            coords[i + 2] = f1*A_new[2] + f2*A_new[5] + f3*A_new[8];
        }
    }

    // Sync nwpw.simulation_cell.unita.
    if (rtdb.contains("nwpw") &&
        rtdb["nwpw"].is_object() &&
        rtdb["nwpw"].contains("simulation_cell") &&
        rtdb["nwpw"]["simulation_cell"].is_object())
    {
        json& cell = rtdb["nwpw"]["simulation_cell"];
        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            for (int i = 0; i < 9; ++i)
                cell["unita"][i] = A_new[i];
        }
    }

    rtdbstring = rtdb.dump();
}



} // namespace pwdft
