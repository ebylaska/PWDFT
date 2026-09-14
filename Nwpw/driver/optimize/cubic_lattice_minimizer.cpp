// cubic_lattice_minimizer.cpp
//
// Lattice optimizer for cubic cells (a = b = c, all angles 90 degrees).
//
// Strategy: isotropic line search on the single free parameter `a`. At each
// step, evaluate E(a*(1+step)) and E(a*(1-step)), accept whichever is lower
// than the current energy, else halve the step. Converges when the step size
// falls below minimum_step or max_steps is reached.
//
// The lattice is modified in place by scaling geometry.unita, geometry.coords,
// and (if present) nwpw.simulation_cell.unita via scale_cubic_cell(). Energy
// and stress come from compute_egs_values(), which drives the electronic
// minimizer supplied as the 4th argument.
//
// Signature matches lattice_minimizer from lattice_minimizer.hpp:
//
//     int (MPI_Comm, std::string&, std::ostream&, electronic_minimizer)

#include "lattice_minimizer.hpp"
#include "lattice_common.hpp"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <string>

#include <mpi.h>

#include "json.hpp"

namespace pwdft {

using json = nlohmann::json;

int cubic_lattice_minimizer(MPI_Comm comm,
                            std::string& rtdbstring,
                            std::ostream& coutput,
                            electronic_minimizer minimizer)
{
    // -- Local diagnostics settings -----------------------------------------
    //
    // We do not have Control2 here, so decide verbosity locally. For now:
    // print the per-step summary. When the driver starts passing an oprint
    // flag (or a Control2 reference), thread it through instead.
    constexpr bool oprint = true;
    const std::string tag = "@@";

    // -- Line-search parameters ---------------------------------------------
    double step = 0.0025;
    constexpr double minimum_step     = 1.0e-5;
    constexpr double minimum_gradient = 1.0e-4;   // used only for reporting
    constexpr int    max_steps        = 25;

    // -- Initial evaluation --------------------------------------------------
    json current_result = compute_egs_values(3, comm, minimizer,
                                             rtdbstring, coutput);
    double current_energy = current_result.at("energy").get<double>();

    int    lstep = 0;
    bool   converged = false;

    // -- Optimization loop ---------------------------------------------------
    for (int istep = 0; istep < max_steps; ++istep)
    {
        // Re-read the current lattice from rtdbstring. It may have been
        // replaced with expanded_rtdb / contracted_rtdb at the end of the
        // previous iteration, so a cached copy would be stale.
        json current_json = json::parse(rtdbstring);

        const std::string geomname =
            (current_json.contains("geometry") &&
             current_json["geometry"].is_string())
                ? current_json["geometry"].get<std::string>()
                : "geometry";

        const json& current_unita =
            current_json["geometries"][geomname]["unita"];
        const double a = current_unita.at(0).get<double>();

        // Stress at the current cell, for reporting and for the gradient
        // convergence metric. Note: stress is only evaluated here, not used
        // to drive the step; the step is driven by the energy line search.
        current_result = compute_egs_values(3, comm, minimizer,
                                            rtdbstring, coutput);
        current_energy = current_result.at("energy").get<double>();
        const json& lstress = current_result.at("lstress");

        const double dE_da = lstress.at(0).get<double>();
        const double dE_db = lstress.at(1).get<double>();
        const double dE_dc = lstress.at(2).get<double>();

        if (oprint)
        {
            coutput << std::defaultfloat << std::setprecision(10)
                    << tag << "Cell step " << istep
                    << " current energy = " << current_energy
                    << " a = " << a
                    << " step = " << step
                    << '\n';
        }

        // -- Trial steps -----------------------------------------------------
        std::string expanded_rtdb   = rtdbstring;
        std::string contracted_rtdb = rtdbstring;

        scale_cubic_cell(expanded_rtdb,   1.0 + step);
        scale_cubic_cell(contracted_rtdb, 1.0 - step);

        json expanded_result   = compute_egs_values(1, comm, minimizer,
                                                    expanded_rtdb, coutput);
        json contracted_result = compute_egs_values(1, comm, minimizer,
                                                    contracted_rtdb, coutput);

        const double expanded_energy   = expanded_result.at("energy").get<double>();
        const double contracted_energy = contracted_result.at("energy").get<double>();

        if (oprint)
        {
            coutput << "\n"
                    << tag << "----------------------------------------------\n"
                    << tag << " PWDFT cubic lattice optimization\n"
                    << tag << "----------------------------------------------\n"
                    << tag << " Step        : " << istep << '\n'
                    << tag << " Energy      : "
                    << std::fixed << std::setprecision(10)
                    << current_energy << " Hartree\n"
                    << tag << " Lattice a   : "
                    << std::fixed << std::setprecision(6)
                    << a << " Bohr ("
                    << std::fixed << std::setprecision(3)
                    << a * 0.529177 << " A)\n"
                    << tag << " dE/da       : " << dE_da << '\n'
                    << tag << " dE/db       : " << dE_db << '\n'
                    << tag << " dE/dc       : " << dE_dc << '\n'
                    << tag << " Step size   : " << step << '\n';

            if (expanded_energy < current_energy &&
                expanded_energy <= contracted_energy)
            {
                coutput << tag << " Action      : Expand, accepted.\n";
            }
            else if (contracted_energy < current_energy)
            {
                coutput << tag << " Action      : Contract, accepted.\n";
            }
            else
            {
                coutput << tag << " Action      : No improvement; halve step.\n";
            }

            coutput << tag << "----------------------------------------------\n";
        }

        // -- Convergence check ----------------------------------------------
        //
        // Gradient norm is computed for reporting and for the convergence
        // criterion, but note that this line search is energy-driven: the
        // step is not derived from the gradient. If you later replace this
        // with a gradient-based (e.g. BFGS) step, this is the place to
        // change both the step computation and the convergence criterion.
        const double grad_norm =
            std::sqrt(dE_da * dE_da + dE_db * dE_db + dE_dc * dE_dc);

        if (step < minimum_step && grad_norm < minimum_gradient)
        {
            converged = true;
            break;
        }

        // -- Accept / reject -------------------------------------------------
        if (expanded_energy < current_energy &&
            expanded_energy <= contracted_energy)
        {
            rtdbstring = std::move(expanded_rtdb);

            if (oprint)
            {
                coutput << std::defaultfloat << std::setprecision(10)
                        << tag << "Accepted expansion, step = " << istep
                        << ", energy = " << expanded_energy << '\n';
            }
        }
        else if (contracted_energy < current_energy)
        {
            rtdbstring = std::move(contracted_rtdb);

            if (oprint)
            {
                coutput << std::defaultfloat << std::setprecision(10)
                        << tag << "Accepted contraction, step = " << istep
                        << ", energy = " << contracted_energy << '\n';
            }
        }
        else
        {
            step *= 0.5;

            if (oprint)
            {
                coutput << std::defaultfloat << std::setprecision(10)
                        << tag << "Rejected both directions, step = " << istep
                        << '\n';
            }

            if (step < minimum_step)
                break;
        }

        lstep = istep;
    }

    // -- Final reporting -----------------------------------------------------
    json final_result = compute_egs_values(3, comm, minimizer,
                                           rtdbstring, coutput);

    const double final_energy = final_result.at("energy").get<double>();
    const json&  final_lstress = final_result.at("lstress");

    json final_json = json::parse(rtdbstring);

    const std::string geomname =
        (final_json.contains("geometry") &&
         final_json["geometry"].is_string())
            ? final_json["geometry"].get<std::string>()
            : "geometry";

    const json& final_unita = final_json["geometries"][geomname]["unita"];
    const double final_a = final_unita.at(0).get<double>();

    if (oprint)
    {
        coutput << '\n'
                << tag << "==============================================\n"
                << tag << " PWDFT cubic lattice optimization COMPLETE\n"
                << tag << "==============================================\n"
                << tag << " Final lattice parameter (a): "
                << std::fixed << std::setprecision(6) << final_a
                << " Bohr = "
                << std::fixed << std::setprecision(3) << final_a * 0.529177
                << " A\n"
                << tag << " Minimum energy (total): "
                << std::fixed << std::setprecision(8) << final_energy
                << " Hartree\n"
                << tag << " Gradients at minimum: "
                << "dE/da = "
                << std::fixed << std::setprecision(5)
                << final_lstress.at(0).get<double>()
                << ", dE/db = " << final_lstress.at(1).get<double>()
                << ", dE/dc = " << final_lstress.at(2).get<double>() << '\n'
                << tag << " Optimization steps taken: " << lstep + 1 << '\n'
                << tag << " Status: "
                << (converged ? "Converged" : "Stopped (max steps reached)")
                << '\n'
                << tag << "==============================================\n";
    }

    return 0;
}

} // namespace pwdft
