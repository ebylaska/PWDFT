// Copyright (C) 2025 Spglib team
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace spglib {
using array_double =
    py::array_t<double, py::array::c_style | py::array::forcecast>;
using array_int = py::array_t<int, py::array::c_style | py::array::forcecast>;
using array_uintp =
    py::array_t<uintptr_t, py::array::c_style | py::array::forcecast>;
using array_size_t =
    py::array_t<size_t, py::array::c_style | py::array::forcecast>;

// Specialized classes wrapping the numpy-like python arrays incorporating
// sanity checks like ndim and shape consistency checks.

class Lattice {
    array_double array;

   public:
    explicit Lattice(array_double &&_array);
    double (*data())[3];
    [[nodiscard]] double const (*data() const)[3];
};

class Rotations {
    array_int array;

   public:
    int const n_operations;
    explicit Rotations(array_int &&_array);
    int (*data())[3][3];
    [[nodiscard]] int const (*data() const)[3][3];
};

class Translations {
    array_double array;

   public:
    int const n_operations;
    explicit Translations(array_double &&_array);
    double (*data())[3];
    [[nodiscard]] double const (*data() const)[3];
};

class Symmetries {
   public:
    Rotations const rotations;
    Translations const translations;
    int const n_operations;
    Symmetries(Rotations &&_rotations, Translations &&_translations);
};

class Positions {
    array_double array;

   public:
    int const n_atoms;
    explicit Positions(array_double &&_array);
    double (*data())[3];
    [[nodiscard]] double const (*data() const)[3];
};

class AtomTypes {
    array_int array;

   public:
    int const n_atoms;
    explicit AtomTypes(array_int &&_array);
    int *data();
    [[nodiscard]]
    int const *data() const;
};

class Magmoms {
    array_double array;

   public:
    int const n_atoms;
    explicit Magmoms(array_double &&_array);
    double *data();
    [[nodiscard]]
    double const *data() const;
};

class Atoms {
   public:
    Positions const positions;
    AtomTypes const types;
    int const n_atoms;
    Atoms(Positions &&_positions, AtomTypes &&_types);
};

class SpglibError : public std::exception {
    std::string msg;

   public:
    SpglibError(std::string_view _msg);
    char const *what() const noexcept override;
};

py::tuple version_tuple();
py::str version_string();
py::str version_full();
py::str commit();
py::dict dataset(Lattice const &lattice, Positions const &positions,
                 AtomTypes const &atom_types, py::int_ hall_number,
                 py::float_ symprec, py::float_ angle_tolerance);
py::dict layer_dataset(Lattice const &lattice, Positions const &positions,
                       AtomTypes const &atom_types, py::int_ aperiodic_dir,
                       py::float_ symprec);
py::dict magnetic_dataset(Lattice const &lattice, Positions const &positions,
                          AtomTypes const &atom_types, array_double magmoms,
                          py::int_ tensor_rank, py::bool_ is_axial,
                          py::float_ symprec, py::float_ angle_tolerance,
                          py::float_ mag_symprec);
py::dict spacegroup_type(py::int_ hall_number);
py::dict spacegroup_type_from_symmetry(Rotations const &rotations,
                                       Translations const &translations,
                                       Lattice const &lattice,
                                       py::float_ symprec);
py::dict magnetic_spacegroup_type(py::int_ uni_number);
py::dict magnetic_spacegroup_type_from_symmetry(
    Rotations const &rotations, Translations const &translations,
    array_int time_reversals, Lattice const &lattice, py::float_ symprec);
py::int_ symmetry_from_database(Rotations &rotations,
                                Translations &translations,
                                py::int_ hall_number);
py::int_ magnetic_symmetry_from_database(Rotations &rotations,
                                         Translations &translations,
                                         array_int time_reversals,
                                         py::int_ uni_number,
                                         py::int_ hall_number);
py::tuple pointgroup(array_int rotations);
py::int_ standardize_cell(Lattice &lattice, Positions &positions,
                          array_int atom_types, py::int_ num_atom,
                          py::int_ to_primative, py::int_ no_idealize,
                          py::float_ symprec, py::float_ angle_tolerance);
py::int_ refine_cell(Lattice &lattice, Positions &positions,
                     AtomTypes &atom_types, py::int_ num_atom,
                     py::float_ symprec, py::float_ angle_tolerance);
py::int_ symmetry(Rotations &rotations, Translations &translations,
                  Lattice const &lattice, Positions const &positions,
                  AtomTypes const &atom_types, py::float_ symprec,
                  py::float_ angle_tolerance);
py::int_ symmetry_with_collinear_spin(
    Rotations &rotations, Translations &translations, array_int equiv_atoms,
    Lattice const &lattice, Positions const &positions,
    AtomTypes const &atom_types, array_double magmoms, py::float_ symprec,
    py::float_ angle_tolerance);
py::int_ symmetry_with_site_tensors(
    Rotations &rotations, Translations &translations, array_int equiv_atoms,
    Lattice &primitive_lattice, array_int spin_flips, Lattice const &lattice,
    Positions const &positions, AtomTypes const &atom_types,
    array_double tensors, py::int_ with_time_reversal, py::int_ is_axial,
    py::float_ symprec, py::float_ angle_tolerance, py::float_ mag_symprec);
py::int_ primitive(Lattice &lattice, Positions &positions,
                   AtomTypes &atom_types, py::float_ symprec,
                   py::float_ angle_tolerance);
py::int_ grid_point_from_address(array_int grid_address, array_int mesh);
py::int_ ir_reciprocal_mesh(array_int grid_address,
                            array_int grid_mapping_table, array_int mesh,
                            array_int is_shift, py::int_ is_time_reversal,
                            Lattice const &lattice, Positions const &positions,
                            AtomTypes const &atom_types, py::float_ symprec);
py::int_ ir_reciprocal_mesh(array_int grid_address,
                            array_size_t grid_mapping_table, array_int mesh,
                            array_int is_shift, py::int_ is_time_reversal,
                            Lattice const &lattice, Positions const &positions,
                            AtomTypes const &atom_types, py::float_ symprec);
py::int_ stabilized_reciprocal_mesh(array_int grid_address,
                                    array_int grid_mapping_table,
                                    array_int mesh, array_int is_shift,
                                    py::int_ is_time_reversal,
                                    Rotations const &rotations,
                                    array_double qpoints);
py::int_ stabilized_reciprocal_mesh(array_int grid_address,
                                    array_size_t grid_mapping_table,
                                    array_int mesh, array_int is_shift,
                                    py::int_ is_time_reversal,
                                    Rotations const &rotations,
                                    array_double qpoints);
void grid_points_by_rotations(array_size_t rot_grid_points,
                              array_int address_orig,
                              Rotations const &rot_reciprocal, array_int mesh,
                              array_int is_shift);
void BZ_grid_points_by_rotations(array_size_t rot_grid_points,
                                 array_int address_orig,
                                 Rotations const &rot_reciprocal,
                                 array_int mesh, array_int is_shift,
                                 array_size_t bz_map);
py::int_ BZ_grid_address(array_int bz_grid_address, array_size_t bz_map,
                         array_int grid_address, array_int mesh,
                         Lattice const &reciprocal_lattice, array_int is_shift);
py::int_ delaunay_reduce(Lattice &lattice, py::float_ symprec);
py::int_ niggli_reduce(Lattice &lattice, py::float_ eps);
py::int_ hall_number_from_symmetry(Rotations const &rotations,
                                   Translations const &translations,
                                   py::float_ symprec);
}  // namespace spglib
