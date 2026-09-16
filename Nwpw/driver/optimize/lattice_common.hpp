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

} // namespace pwdft
