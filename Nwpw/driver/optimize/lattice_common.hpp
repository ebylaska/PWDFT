// lattice_common.hpp
#pragma once

#include <array>
#include <mpi.h>
#include <ostream>
#include <string>

#include "json.hpp"
#include "lattice_minimizer.hpp"   // for electronic_minimizer

namespace pwdft {

void print_lattice_state(const nlohmann::json& rtdbjson,
                         std::ostream& coutput,
                         const std::string& tag);

bool   read_unita(const nlohmann::json& value, std::array<double, 9>& unita);
void   write_unita(nlohmann::json& value, const std::array<double, 9>& unita);
double unita_relative_difference(const std::array<double, 9>& current,
                                 const std::array<double, 9>& frozen);

void scale_cubic_cell(std::string& rtdbstring, double scale);

nlohmann::json compute_egs_values(int option,
                                  MPI_Comm comm,
                                  electronic_minimizer minimizer,
                                  std::string& rtdbstring,
                                  std::ostream& coutput);

// lattice_common.hpp (additions)
std::pair<double, double> read_tetragonal_lattice(const std::string& rtdbstring);
void set_tetragonal_cell(std::string& rtdbstring, double a_new, double c_new);

std::pair<double, double> read_hexagonal_lattice(const std::string& rtdbstring);
void set_hexagonal_cell(std::string& rtdbstring, double a_new, double c_new);

std::array<double, 3> read_orthorhombic_lattice(const std::string& rtdbstring);
void set_orthorhombic_cell(std::string& rtdbstring, double a_new, double b_new, double c_new);

struct MonoclinicLattice {
    double a;
    double b;
    double c;
    double beta_rad;
};

MonoclinicLattice read_monoclinic_lattice(const std::string& rtdbstring);

void set_monoclinic_cell(std::string& rtdbstring, double a_new, double b_new, double c_new, double beta_new_rad);

struct RhombohedralLattice {
    double a;
    double alpha_rad;
};

RhombohedralLattice read_rhombohedral_lattice(const std::string& rtdbstring);

void set_rhombohedral_cell(std::string& rtdbstring, double a_new, double alpha_new_rad);

struct TriclinicLattice {
    double a, b, c;             // Bohr
    double alpha_rad, beta_rad, gamma_rad;
};

TriclinicLattice read_triclinic_lattice(const std::string& rtdbstring);

void set_triclinic_cell(std::string& rtdbstring,
                        double a, double b, double c,
                        double alpha_rad, double beta_rad, double gamma_rad);

// 3x3 inverse; returns false if singular.
bool invert_3x3(const std::array<double, 9>& M,
                std::array<double, 9>& Minv);

} // namespace pwdft
