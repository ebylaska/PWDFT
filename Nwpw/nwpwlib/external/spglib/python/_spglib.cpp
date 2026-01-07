// Copyright (C) 2025 Spglib team
// SPDX-License-Identifier: BSD-3-Clause

#include <pybind11/pybind11.h>

#include "py_bindings.h"

// py::mod_gil_not_used was added in pybind11 2.13 and
#if (PYBIND11_VERSION_MAJOR == 2 && PYBIND11_VERSION_MINOR >= 13) || \
    PYBIND11_VERSION_MAJOR >= 3
PYBIND11_MODULE(_spglib, module, py::mod_gil_not_used()) {
#else
PYBIND11_MODULE(_spglib, module) {
#endif
    using namespace py::literals;
    using namespace spglib;
    module.doc() = "Spglib compiled bindings.";

    // Register the cpp Spglib error with a base class of the python one.
    auto base_error = py::module_::import("spglib.error").attr("SpglibError");
    py::register_exception<SpglibError>(module, "SpglibCppError",
                                        base_error.ptr());

    py::class_<Lattice>(module, "Lattice").def(py::init<array_double>());
    py::implicitly_convertible<array_double, Lattice>();
    py::class_<Rotations>(module, "Rotations").def(py::init<array_int>());
    py::implicitly_convertible<array_int, Rotations>();
    py::class_<Translations>(module, "Translations")
        .def(py::init<array_double>());
    py::implicitly_convertible<array_double, Translations>();
    py::class_<Symmetries>(module, "Symmetries")
        .def(py::init<Rotations, Translations>());
    py::class_<Positions>(module, "Positions").def(py::init<array_double>());
    py::implicitly_convertible<array_double, Positions>();
    py::class_<AtomTypes>(module, "AtomTypes").def(py::init<array_int>());
    py::implicitly_convertible<array_int, AtomTypes>();
    py::class_<Magmoms>(module, "Magmoms").def(py::init<array_double>());
    py::implicitly_convertible<array_double, Magmoms>();
    py::class_<Atoms>(module, "Atoms").def(py::init<Positions, AtomTypes>());

    module.def("version_tuple", spglib::version_tuple, "");
    module.def("version_string", spglib::version_string, "");
    module.def("version_full", spglib::version_full, "");
    module.def("commit", spglib::commit, "");
    module.def("dataset", spglib::dataset, "");
    module.def("layer_dataset", spglib::layer_dataset, "");
    module.def("magnetic_dataset", spglib::magnetic_dataset, "");
    module.def("spacegroup_type", spglib::spacegroup_type, "");
    module.def("spacegroup_type_from_symmetry",
               spglib::spacegroup_type_from_symmetry, "");
    module.def("magnetic_spacegroup_type", spglib::magnetic_spacegroup_type,
               "");
    module.def("magnetic_spacegroup_type_from_symmetry",
               spglib::magnetic_spacegroup_type_from_symmetry, "");
    module.def("symmetry_from_database", spglib::symmetry_from_database, "");
    module.def("magnetic_symmetry_from_database",
               spglib::magnetic_symmetry_from_database, "");
    module.def("pointgroup", spglib::pointgroup, "");
    module.def("standardize_cell", spglib::standardize_cell, "");
    module.def("refine_cell", spglib::refine_cell, "");
    module.def("symmetry", spglib::symmetry, "");
    module.def("symmetry_with_collinear_spin",
               spglib::symmetry_with_collinear_spin, "");
    module.def("symmetry_with_site_tensors", spglib::symmetry_with_site_tensors,
               "");
    module.def("primitive", spglib::primitive, "");
    module.def("grid_point_from_address", spglib::grid_point_from_address, "");
    module.def(
        "ir_reciprocal_mesh",
        py::overload_cast<array_int, array_int, array_int, array_int, py::int_,
                          Lattice const &, Positions const &, AtomTypes const &,
                          py::float_>(spglib::ir_reciprocal_mesh),
        "", py::arg("grid_address"), py::arg("grid_mapping_table").noconvert(),
        py::arg("mesh"), py::arg("is_shift"), py::arg("is_time_reversal"),
        py::arg("lattice"), py::arg("positions"), py::arg("atom_types"),
        py::arg("symprec"));
    module.def("ir_reciprocal_mesh",
               py::overload_cast<array_int, array_size_t, array_int, array_int,
                                 py::int_, Lattice const &, Positions const &,
                                 AtomTypes const &, py::float_>(
                   spglib::ir_reciprocal_mesh),
               "", py::arg("grid_address"),
               py::arg("grid_mapping_table").noconvert(), py::arg("mesh"),
               py::arg("is_shift"), py::arg("is_time_reversal"),
               py::arg("lattice"), py::arg("positions"), py::arg("atom_types"),
               py::arg("symprec"));
    module.def("stabilized_reciprocal_mesh",
               py::overload_cast<array_int, array_int, array_int, array_int,
                                 py::int_, Rotations const &, array_double>(
                   spglib::stabilized_reciprocal_mesh),
               "", py::arg("grid_address"),
               py::arg("grid_mapping_table").noconvert(), py::arg("mesh"),
               py::arg("is_shift"), py::arg("is_time_reversal"),
               py::arg("rotations"), py::arg("qpoints"));
    module.def("stabilized_reciprocal_mesh",
               py::overload_cast<array_int, array_size_t, array_int, array_int,
                                 py::int_, Rotations const &, array_double>(
                   spglib::stabilized_reciprocal_mesh),
               "", py::arg("grid_address"),
               py::arg("grid_mapping_table").noconvert(), py::arg("mesh"),
               py::arg("is_shift"), py::arg("is_time_reversal"),
               py::arg("rotations"), py::arg("qpoints"));
    module.def("grid_points_by_rotations", spglib::grid_points_by_rotations,
               "");
    module.def("BZ_grid_points_by_rotations",
               spglib::BZ_grid_points_by_rotations, "");
    module.def("BZ_grid_address", spglib::BZ_grid_address, "");
    module.def("delaunay_reduce", spglib::delaunay_reduce, "");
    module.def("niggli_reduce", spglib::niggli_reduce, "");
    module.def("hall_number_from_symmetry", spglib::hall_number_from_symmetry,
               "");
}
