// lattice_common.cpp
//
// Shared helpers used by the driver and by every per-system lattice
// minimizer. Nothing here is driver-specific; nothing here is
// crystal-system-specific except scale_cubic_cell, which will eventually
// be joined by scale_tetragonal_cell, scale_orthorhombic_cell, etc., or
// replaced by a generic set_lattice(...).
//
// These functions were previously static members of driver_optimize.cpp.

#include "lattice_common.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>

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

} // namespace pwdft
