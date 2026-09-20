// cubic_lattice_minimizer.cpp
//
// Cubic cells (a = b = c, all angles 90 degrees). One free parameter: a.
//
// Thin wrapper around run_quasi_newton<1>. The coordinate is x = log(a).
// The gradient along the cubic constraint is
//
//     dE/dlog(a) = a * (dE/da + dE/db + dE/dc).

#include "lattice_quasi_newton.hpp"

#include <array>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>

namespace pwdft {

namespace {

double cubic_lattice_parameter(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname =
        (rtdb.contains("geometry") && rtdb["geometry"].is_string())
            ? rtdb["geometry"].get<std::string>()
            : "geometry";
    return rtdb.at("geometries").at(geomname).at("unita").at(0).get<double>();
}

} // namespace

int cubic_lattice_minimizer(MPI_Comm comm,
                            std::string& rtdbstring,
                            std::ostream& coutput,
                            electronic_minimizer minimizer,
                            const LatticeContext& ctx)
{
    auto read_x = [](const std::string& rtdb) {
        std::array<double, 1> x{ std::log(cubic_lattice_parameter(rtdb)) };
        return x;
    };

    auto write_x = [](std::string& rtdb, const std::array<double, 1>& x) {
        const double a_old = cubic_lattice_parameter(rtdb);
        const double a_new = std::exp(x[0]);
        if (!(a_new > 0.0))
            throw std::runtime_error("cubic: non-positive a");
        scale_cubic_cell(rtdb, a_new / a_old);
    };

    auto compute_g = [](const std::array<double, 1>& x, const json& lstress) {
        const double a = std::exp(x[0]);
        const double dE_da =
            lstress.at(0).get<double>() +
            lstress.at(1).get<double>() +
            lstress.at(2).get<double>();
        return std::array<double, 1>{ a * dE_da };
    };

    auto print_step = [](std::ostream& out, const std::string& tag, int istep,
                         const std::array<double, 1>& x,
                         const std::array<double, 1>& g,
                         double E, double trust) {
        const double a = std::exp(x[0]);
        out << '\n'
            << tag << "----------------------------------------------\n"
            << tag << " Step          : " << istep << '\n'
            << tag << " Energy        : " << std::fixed << std::setprecision(10)
                   << E << " Hartree\n"
            << tag << " Lattice a     : " << std::fixed << std::setprecision(6)
                   << a << " Bohr ("
                   << std::fixed << std::setprecision(4)
                   << a * bohr_to_angstrom << " A)\n"
            << tag << " dE/dlog(a)    : " << std::defaultfloat << std::setprecision(10)
                   << g[0] << '\n'
            << tag << " Trust radius  : " << trust << '\n';
    };

    auto print_final = [](std::ostream& out, const std::string& tag,
                          const std::array<double, 1>& x, double /*E*/) {
        const double a = std::exp(x[0]);
        out << tag << " Final lattice a: " << std::fixed << std::setprecision(6)
                 << a << " Bohr ("
                 << std::fixed << std::setprecision(4)
                 << a * bohr_to_angstrom << " A)\n";
    };

    return run_quasi_newton<1>(rtdbstring, comm, minimizer, coutput, ctx,
                               "cubic",
                               read_x, write_x, compute_g,
                               print_step, print_final);
}

} // namespace pwdft
