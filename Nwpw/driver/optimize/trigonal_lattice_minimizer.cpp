// trigonal_lattice_minimizer.cpp
//
// Trigonal R setting: a = b = c, alpha = beta = gamma. Two free parameters
// (a, alpha). Coordinates x = (log a, alpha_rad).
//
// Gradients:
//     dE/dlog(a) = a * (dE/da + dE/db + dE/dc)
//     dE/dalpha  = dE/dalpha + dE/dbeta + dE/dgamma  (all equal in R)

#include "lattice_quasi_newton.hpp"

#include <array>
#include <cmath>
#include <ostream>
#include <string>

namespace pwdft {

int trigonal_lattice_minimizer(MPI_Comm comm,
                               std::string& rtdbstring,
                               std::ostream& coutput,
                               electronic_minimizer minimizer,
                               const LatticeContext& ctx)
{
    constexpr double angle_unit_factor = 1.0; // Ha/rad; use pi/180 for Ha/deg

    auto read_x = [](const std::string& rtdb) {
        const auto rl = read_rhombohedral_lattice(rtdb);
        return std::array<double, 2>{ std::log(rl.a), rl.alpha_rad };
    };

    auto write_x = [](std::string& rtdb, const std::array<double, 2>& x) {
        set_rhombohedral_cell(rtdb, std::exp(x[0]), x[1]);
    };

    auto compute_g = [angle_unit_factor](const std::array<double, 2>& x,
                                         const json& lstress) {
        const double a = std::exp(x[0]);
        const double dE_da =
            lstress.at(0).get<double>() +
            lstress.at(1).get<double>() +
            lstress.at(2).get<double>();
        const double dE_dalpha =
            lstress.at(3).get<double>() +
            lstress.at(4).get<double>() +
            lstress.at(5).get<double>();
        return std::array<double, 2>{
            a * dE_da,
            angle_unit_factor * dE_dalpha
        };
    };

    auto print_step = [](std::ostream& out, const std::string& tag, int istep,
                         const std::array<double, 2>& x,
                         const std::array<double, 2>& g,
                         double E, double trust) {
        const double a  = std::exp(x[0]);
        const double ar = x[1];
        const double ad = ar * 180.0 / M_PI;
        out << '\n'
            << tag << "----------------------------------------------\n"
            << tag << " Step          : " << istep << '\n'
            << tag << " Energy        : " << std::fixed << std::setprecision(10)
                   << E << " Hartree\n"
            << tag << " Lattice a     : " << std::fixed << std::setprecision(6)
                   << a << " Bohr ("
                   << std::fixed << std::setprecision(4)
                   << a * bohr_to_angstrom << " A)\n"
            << tag << " alpha         : " << std::fixed << std::setprecision(4)
                   << ad << " deg ("
                   << ar << " rad)\n"
            << tag << " dE/dlog(a)    : " << std::defaultfloat << std::setprecision(10)
                   << g[0] << '\n'
            << tag << " dE/dalpha     : " << g[1] << '\n'
            << tag << " Trust radius  : " << trust << '\n';
    };

    auto print_final = [](std::ostream& out, const std::string& tag,
                          const std::array<double, 2>& x, double /*E*/) {
        const double a  = std::exp(x[0]);
        const double ar = x[1];
        const double ad = ar * 180.0 / M_PI;
        out << tag << " Final lattice a: " << std::fixed << std::setprecision(6)
                 << a << " Bohr ("
                 << std::fixed << std::setprecision(4)
                 << a * bohr_to_angstrom << " A)\n"
            << tag << " Final alpha    : " << std::fixed << std::setprecision(4)
                 << ad << " deg ("
                 << ar << " rad)\n";
    };

    return run_quasi_newton<2>(rtdbstring, comm, minimizer, coutput, ctx,
                               "trigonal (R)",
                               read_x, write_x, compute_g,
                               print_step, print_final);
}

} // namespace pwdft
