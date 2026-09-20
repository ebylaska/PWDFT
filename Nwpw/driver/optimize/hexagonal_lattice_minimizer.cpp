// hexagonal_lattice_minimizer.cpp
//
// Hexagonal cells (a = b, gamma = 120, alpha = beta = 90; c independent).
// Free parameters: a and c.
//
// Coordinates x = (log a, log c). Gradients:
//
//     dE/dlog(a) = a * (dE/da + dE/db)
//     dE/dlog(c) = c * dE/dc

#include "lattice_quasi_newton.hpp"

namespace pwdft {

int hexagonal_lattice_minimizer(MPI_Comm comm,
                                std::string& rtdbstring,
                                std::ostream& coutput,
                                electronic_minimizer minimizer,
                                const LatticeContext& ctx)
{
    auto read_x = [](const std::string& rtdb) {
        const auto ac = read_hexagonal_lattice(rtdb);
        return std::array<double, 2>{ std::log(ac.first),
                                      std::log(ac.second) };
    };

    auto write_x = [](std::string& rtdb, const std::array<double, 2>& x) {
        set_hexagonal_cell(rtdb, std::exp(x[0]), std::exp(x[1]));
    };

    auto compute_g = [](const std::array<double, 2>& x, const json& lstress) {
        const double a = std::exp(x[0]);
        const double c = std::exp(x[1]);
        const double dE_da =
            lstress.at(0).get<double>() + lstress.at(1).get<double>();
        const double dE_dc = lstress.at(2).get<double>();
        return std::array<double, 2>{ a * dE_da, c * dE_dc };
    };

    auto print_step = [](std::ostream& out, const std::string& tag, int istep,
                     const std::array<double, 2>& x,
                     const std::array<double, 2>& g,
                     double E, double trust) {
        const double a = std::exp(x[0]);
        const double c = std::exp(x[1]);
        out << '\n'
            << tag << "----------------------------------------------\n"
            << tag << " Step          : " << istep << '\n'
            << tag << " Energy        : " << std::fixed << std::setprecision(10)
                   << E << " Hartree\n"
            << tag << " Lattice a     : " << std::fixed << std::setprecision(6)
                   << a << " Bohr ("
                   << std::fixed << std::setprecision(4)
                   << a * bohr_to_angstrom << " A)\n"
            << tag << " Lattice c     : " << std::fixed << std::setprecision(6)
                   << c << " Bohr ("
                   << std::fixed << std::setprecision(4)
                   << c * bohr_to_angstrom << " A)\n"
            << tag << " dE/dlog(a)    : " << std::defaultfloat << std::setprecision(10)
                   << g[0] << '\n'
            << tag << " dE/dlog(c)    : " << g[1] << '\n'
            << tag << " Trust radius  : " << trust << '\n';
    };

    auto print_final = [](std::ostream& out, const std::string& tag,
                      const std::array<double, 2>& x, double E) {
        const double a = std::exp(x[0]);
        const double c = std::exp(x[1]);
        out << tag << " Final lattice a: " << std::fixed << std::setprecision(6)
                 << a << " Bohr ("
                 << std::fixed << std::setprecision(4)
                 << a * bohr_to_angstrom << " A)\n"
            << tag << " Final lattice c: " << std::fixed << std::setprecision(6)
                 << c << " Bohr ("
                 << std::fixed << std::setprecision(4)
                 << c * bohr_to_angstrom << " A)\n";
    };


    return run_quasi_newton<2>(rtdbstring, comm, minimizer, coutput, ctx,
                               "hexagonal",
                               read_x, write_x, compute_g, print_step, print_final);
}

} // namespace pwdft
