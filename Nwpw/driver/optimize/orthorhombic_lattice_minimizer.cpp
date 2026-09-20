// orthorhombic_lattice_minimizer.cpp
//
// Orthorhombic cells (alpha = beta = gamma = 90; a, b, c independent).
// Free parameters: a, b, c.
//
// Coordinates x = (log a, log b, log c). Gradients:
//
//     dE/dlog(a) = a * dE/da
//     dE/dlog(b) = b * dE/db
//     dE/dlog(c) = c * dE/dc

#include "lattice_quasi_newton.hpp"

#include <array>
#include <cmath>
#include <ostream>
#include <string>

namespace pwdft {

int orthorhombic_lattice_minimizer(MPI_Comm comm,
                                   std::string& rtdbstring,
                                   std::ostream& coutput,
                                   electronic_minimizer minimizer,
                                   const LatticeContext& ctx)
{
    auto read_x = [](const std::string& rtdb) {
        const auto abc = read_orthorhombic_lattice(rtdb);
        return std::array<double, 3>{
            std::log(abc[0]),
            std::log(abc[1]),
            std::log(abc[2])
        };
    };

    auto write_x = [](std::string& rtdb, const std::array<double, 3>& x) {
        set_orthorhombic_cell(rtdb, std::exp(x[0]), std::exp(x[1]), std::exp(x[2]));
    };

    auto compute_g = [](const std::array<double, 3>& x, const json& lstress) {
        return std::array<double, 3>{
            std::exp(x[0]) * lstress.at(0).get<double>(),
            std::exp(x[1]) * lstress.at(1).get<double>(),
            std::exp(x[2]) * lstress.at(2).get<double>()
        };
    };

    auto print_step = [](std::ostream& out, const std::string& tag, int istep,
                         const std::array<double, 3>& x,
                         const std::array<double, 3>& g,
                         double E, double trust) {
        const double a = std::exp(x[0]);
        const double b = std::exp(x[1]);
        const double c = std::exp(x[2]);
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
            << tag << " dE/dlog(a)    : " << std::defaultfloat << std::setprecision(10)
                   << g[0] << '\n'
            << tag << " dE/dlog(b)    : " << g[1] << '\n'
            << tag << " dE/dlog(c)    : " << g[2] << '\n'
            << tag << " Trust radius  : " << trust << '\n';
    };

    auto print_final = [](std::ostream& out, const std::string& tag,
                          const std::array<double, 3>& x, double /*E*/) {
        const double a = std::exp(x[0]);
        const double b = std::exp(x[1]);
        const double c = std::exp(x[2]);
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
                 << c * bohr_to_angstrom << " A)\n";
    };

    return run_quasi_newton<3>(rtdbstring, comm, minimizer, coutput, ctx,
                               "orthorhombic",
                               read_x, write_x, compute_g,
                               print_step, print_final);
}

} // namespace pwdft
