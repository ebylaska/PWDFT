// triclinic_lattice_minimizer.cpp
//
// Triclinic cells: a, b, c, alpha, beta, gamma all free.
//
// Coordinates: x = (log a, log b, log c, alpha_rad, beta_rad, gamma_rad)
//
// Gradients:
//     dE/dlog(a) = a * lstress[0]
//     dE/dlog(b) = b * lstress[1]
//     dE/dlog(c) = c * lstress[2]
//     dE/dalpha  = lstress[3]
//     dE/dbeta   = lstress[4]
//     dE/dgamma  = lstress[5]
//
// Unit convention: assumes lstress[3..5] are Ha/rad. Set angle_unit_factor
// to pi/180 if they are Ha/deg.

#include "lattice_quasi_newton.hpp"

#include <array>
#include <cmath>
#include <ostream>
#include <string>

namespace pwdft {

int triclinic_lattice_minimizer(MPI_Comm comm,
                                std::string& rtdbstring,
                                std::ostream& coutput,
                                electronic_minimizer minimizer,
                                const LatticeContext& ctx)
{
    constexpr double angle_unit_factor = 1.0; // Ha/rad -> 1.0 ; Ha/deg -> pi/180

    auto read_x = [](const std::string& rtdb) {
        const auto tl = read_triclinic_lattice(rtdb);
        return std::array<double, 6>{
            std::log(tl.a),
            std::log(tl.b),
            std::log(tl.c),
            tl.alpha_rad,
            tl.beta_rad,
            tl.gamma_rad
        };
    };

    auto write_x = [](std::string& rtdb, const std::array<double, 6>& x) {
        set_triclinic_cell(rtdb,
                           std::exp(x[0]),
                           std::exp(x[1]),
                           std::exp(x[2]),
                           x[3], x[4], x[5]);
    };

    auto compute_g = [angle_unit_factor](const std::array<double, 6>& x,
                                         const json& lstress) {
        return std::array<double, 6>{
            std::exp(x[0]) * lstress.at(0).get<double>(),
            std::exp(x[1]) * lstress.at(1).get<double>(),
            std::exp(x[2]) * lstress.at(2).get<double>(),
            angle_unit_factor * lstress.at(3).get<double>(),
            angle_unit_factor * lstress.at(4).get<double>(),
            angle_unit_factor * lstress.at(5).get<double>()
        };
    };

    auto print_step = [](std::ostream& out, const std::string& tag, int istep,
                         const std::array<double, 6>& x,
                         const std::array<double, 6>& g,
                         double E, double trust) {
        const double a = std::exp(x[0]);
        const double b = std::exp(x[1]);
        const double c = std::exp(x[2]);
        const double al = x[3] * 180.0 / M_PI;
        const double be = x[4] * 180.0 / M_PI;
        const double ga = x[5] * 180.0 / M_PI;
        out << '\n'
            << tag << "----------------------------------------------\n"
            << tag << " Step          : " << istep << '\n'
            << tag << " Energy        : " << std::fixed << std::setprecision(10)
                   << E << " Hartree\n"
            << tag << " Lattice a     : " << std::fixed << std::setprecision(6)
                   << a << " Bohr ("
                   << std::fixed << std::setprecision(4)
                   << a * bohr_to_angstrom << " A)\n"
            << tag << " Lattice b     : " << std::fixed << std::setprecision(6)
                   << b << " Bohr ("
                   << std::fixed << std::setprecision(4)
                   << b * bohr_to_angstrom << " A)\n"
            << tag << " Lattice c     : " << std::fixed << std::setprecision(6)
                   << c << " Bohr ("
                   << std::fixed << std::setprecision(4)
                   << c * bohr_to_angstrom << " A)\n"
            << tag << " alpha         : " << std::fixed << std::setprecision(4)
                   << al << " deg\n"
            << tag << " beta          : " << std::fixed << std::setprecision(4)
                   << be << " deg\n"
            << tag << " gamma         : " << std::fixed << std::setprecision(4)
                   << ga << " deg\n"
            << tag << " dE/dlog(a)    : " << std::defaultfloat << std::setprecision(10)
                   << g[0] << '\n'
            << tag << " dE/dlog(b)    : " << g[1] << '\n'
            << tag << " dE/dlog(c)    : " << g[2] << '\n'
            << tag << " dE/dalpha     : " << g[3] << '\n'
            << tag << " dE/dbeta      : " << g[4] << '\n'
            << tag << " dE/dgamma     : " << g[5] << '\n'
            << tag << " Trust radius  : " << trust << '\n';
    };

    auto print_final = [](std::ostream& out, const std::string& tag,
                          const std::array<double, 6>& x, double /*E*/) {
        const double a = std::exp(x[0]);
        const double b = std::exp(x[1]);
        const double c = std::exp(x[2]);
        const double al = x[3] * 180.0 / M_PI;
        const double be = x[4] * 180.0 / M_PI;
        const double ga = x[5] * 180.0 / M_PI;
        out << tag << " Final lattice a: " << std::fixed << std::setprecision(6)
                 << a << " Bohr ("
                 << std::fixed << std::setprecision(4)
                 << a * bohr_to_angstrom << " A)\n"
            << tag << " Final lattice b: " << std::fixed << std::setprecision(6)
                 << b << " Bohr ("
                 << std::fixed << std::setprecision(4)
                 << b * bohr_to_angstrom << " A)\n"
            << tag << " Final lattice c: " << std::fixed << std::setprecision(6)
                 << c << " Bohr ("
                 << std::fixed << std::setprecision(4)
                 << c * bohr_to_angstrom << " A)\n"
            << tag << " Final alpha    : " << std::fixed << std::setprecision(4)
                 << al << " deg\n"
            << tag << " Final beta     : " << std::fixed << std::setprecision(4)
                 << be << " deg\n"
            << tag << " Final gamma    : " << std::fixed << std::setprecision(4)
                 << ga << " deg\n";
    };

    return run_quasi_newton<6>(rtdbstring, comm, minimizer, coutput, ctx,
                               "triclinic",
                               read_x, write_x, compute_g,
                               print_step, print_final);
}

} // namespace pwdft
