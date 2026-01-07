"""Test of spacegroup_type_from_symmetry."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from spglib import (
    get_spacegroup_type,
    get_spacegroup_type_from_symmetry,
    get_symmetry_dataset,
)

if TYPE_CHECKING:
    pass


@pytest.mark.parametrize("lattice_None", [True, False])
def test_spacegroup_type_from_symmetry(
    lattice_None: bool,
    crystal_data_dataset,
):
    """Test spacegroup_type_from_symmetry."""
    if lattice_None:
        lattice = None
    else:
        lattice = crystal_data_dataset["crystal_data"].cell[0]

    dataset = crystal_data_dataset["dataset"]
    spgtype = get_spacegroup_type_from_symmetry(
        dataset.rotations,
        dataset.translations,
        lattice=lattice,
        symprec=crystal_data_dataset["symprec"],
    )

    assert spgtype.number == dataset.number
    spgtype_ref = get_spacegroup_type(dataset.hall_number)
    assert spgtype == spgtype_ref


def test_spacegroup_type_from_symmetry_with_oblique_basis_vectors():
    """Test spacegroup_type_from_symmetry with oblique basis vectors."""
    # Structure derived from 1509692 in the crystallography open database.
    lattice = [
        [0.0000, 6.4345, 1.8319],
        [5.245, -6.4345, 0.0000],
        [0.0000, 0.0000, -3.6638],
    ]
    position = [
        [0.6416, 0.7350, 0.0708],
        [0.8284, 0.2350, 0.4142],
        [0.3584, 0.2650, 0.4292],
        [0.1716, 0.7650, 0.5858],
        [0.6994, 0.5334, 0.2150],
        [0.3674, 0.0334, 0.2990],
        [0.3006, 0.4666, 0.5156],
        [0.0000, 0.0000, 0.5902],
        [0.0000, 0.5000, 0.8402],
        [0.6326, 0.9666, 0.9316],
    ]
    types = [1, 1, 1, 1, 6, 6, 6, 6, 6, 6]

    dataset = get_symmetry_dataset((lattice, position, types), symprec=1e-5)
    assert dataset.hall_number == 212

    spg_type = get_spacegroup_type_from_symmetry(
        dataset.rotations,
        dataset.translations,
        lattice=lattice,
        symprec=1e-5,
    )
    assert spg_type.hall_number == 212
