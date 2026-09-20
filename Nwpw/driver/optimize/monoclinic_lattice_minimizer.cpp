// monoclinic_lattice_minimizer.cpp
//
// Monoclinic cells (alpha = gamma = 90; beta free). Four free parameters:
// a, b, c, beta.
//
// Coordinates:
//     x = (log(a), log(b), log(c), beta_rad)
//
// Gradients along the monoclinic constraint:
//     dE/dlog(a) = a * dE/da
//     dE/dlog(b) = b * dE/db
//     dE/dlog(c) = c * dE/dc
//     dE/d(beta_rad) = dE/dbeta
//
// Unit convention for beta: this assumes lstress[4] is dE/dbeta in Ha/rad.
// If PWDFT reports dE/dbeta in Ha/deg, multiply by (pi/180) below.

#include "lattice_quasi_newton.hpp"

#include <array>
#include <cmath>
#include <ostream>
#include <string>

namespace pwdft {

int monoclinic_lattice_minimizer(MPI_Comm comm,
                                 std::string& rtdbstring,
                                 std::ostream& coutput,
                                 electronic_minimizer minimizer,
                                 const LatticeContext& ctx)
{
    // Set to 1.0 if lstress[4] is Ha/rad; to (M_PI/180.0) if it's Ha/deg.
    constexpr double angle_unit_factor = 1.0;

    auto read_x = [](const std::string& rtdb) {
        const auto ml = read_monoclinic_lattice(rtdb);
        return std::array<double, 4>{
            std::log(ml.a),
            std::log(ml.b),
            std::log(ml.c),
            ml.beta_rad
        };
    };

    auto write_x = [](std::string& rtdb, const std::array<double, 4>& x) {
        set_monoclinic_cell(rtdb,
                            std::exp(x[0]),
                            std::exp(x[1]),
                            std::exp(x[2]),
                            x[3]);
    };

    auto compute_g = [angle_unit_factor](const std::array<double, 4>& x,
                                         const json& lstress) {
        return std::array<double, 4>{
            std::exp(x[0]) * lstress.at(0).get<double>(),
            std::exp(x[1]) * lstress.at(1).get<double>(),
            std::exp(x[2]) * lstress.at(2).get<double>(),
            angle_unit_factor * lstress.at(4).get<double>()
        };
    };

    auto print_step = [](std::ostream& out, const std::string& tag, int istep,
                         const std::array<double, 4>& x,
                         const std::array<double, 4>& g,
                         double E, double trust) {
        const double a  = std::exp(x[0]);
        const double b  = std::exp(x[1]);
        const double c  = std::exp(x[2]);
        const double br = x[3];
        const double bd = br * 180.0 / M_PI;
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
            << tag << " beta          : " << std::fixed << std::setprecision(4)
                   << bd << " deg ("
                   << br << " rad)\n"
            << tag << " dE/dlog(a)    : " << std::defaultfloat << std::setprecision(10)
                   << g[0] << '\n'
            << tag << " dE/dlog(b)    : " << g[1] << '\n'
            << tag << " dE/dlog(c)    : " << g[2] << '\n'
            << tag << " dE/dbeta      : " << g[3] << '\n'
            << tag << " Trust radius  : " << trust << '\n';
    };

    auto print_final = [](std::ostream& out, const std::string& tag,
                          const std::array<double, 4>& x, double /*E*/) {
        const double a  = std::exp(x[0]);
        const double b  = std::exp(x[1]);
        const double c  = std::exp(x[2]);
        const double br = x[3];
        const double bd = br * 180.0 / M_PI;
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
            << tag << " Final beta     : " << std::fixed << std::setprecision(4)
                 << bd << " deg ("
                 << br << " rad)\n";
    };

    return run_quasi_newton<4>(rtdbstring, comm, minimizer, coutput, ctx,
                               "monoclinic",
                               read_x, write_x, compute_g,
                               print_step, print_final);
}

} // namespace pwdft
