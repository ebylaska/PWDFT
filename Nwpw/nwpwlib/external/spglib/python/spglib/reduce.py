"""Lattice reduction operations."""
# Copyright (C) 2015 Atsushi Togo
# This file is part of spglib.
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from . import _spglib
from .error import _set_no_error, _set_or_throw_error

__all__ = [
    "delaunay_reduce",
    "niggli_reduce",
]


def delaunay_reduce(
    lattice: ArrayLike[np.double], eps: float = 1e-5
) -> np.ndarray | None:
    r"""Run Delaunay reduction. When the search failed, `None` is returned.

    The transformation from original basis vectors
    :math:`( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )`
    to final basis vectors :math:`( \mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )` is achieved by linear
    combination of basis vectors with integer coefficients without
    rotating coordinates. Therefore the transformation matrix is obtained
    by :math:`\mathbf{P} = ( \mathbf{a} \; \mathbf{b} \; \mathbf{c} ) ( \mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )^{-1}` and the matrix
    elements have to be almost integers.

    The algorithm is found in the international tables for crystallography volume A.

    Parameters
    ----------
    lattice: ndarray, (3, 3)
        Lattice parameters in the form of

        .. code-block::

            [[a_x, a_y, a_z],
                [b_x, b_y, b_z],
                [c_x, c_y, c_z]]

    eps: float
        Tolerance parameter, but unlike `symprec` the unit is not a length.
        Tolerance to check if volume is close to zero or not and
        if two basis vectors are orthogonal by the value of dot
        product being close to zero or not.

    Returns
    -------
    delaunay_lattice: ndarray, (3, 3)
        Reduced lattice parameters are given as a numpy 'double' array:

        .. code-block::

            [[a_x, a_y, a_z],
             [b_x, b_y, b_z],
             [c_x, c_y, c_z]]

    Notes
    -----
    .. versionadded:: 1.9.4

    """  # noqa: E501
    _set_no_error()

    delaunay_lattice = np.array(np.transpose(lattice), dtype="double", order="C")
    try:
        _spglib.delaunay_reduce(delaunay_lattice, float(eps))
    except Exception as exc:
        _set_or_throw_error(exc)
        return None
    return np.array(np.transpose(delaunay_lattice), dtype="double", order="C")


def niggli_reduce(
    lattice: ArrayLike[np.double], eps: float = 1e-5
) -> np.ndarray | None:
    r"""Run Niggli reduction. When the search failed, ``None`` is returned.

    The transformation from original basis vectors :math:`( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )` to final basis vectors :math:`( \mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )` is achieved by linear
    combination of basis vectors with integer coefficients without
    rotating coordinates. Therefore the transformation matrix is obtained
    by :math:`\mathbf{P} = ( \mathbf{a} \; \mathbf{b} \; \mathbf{c} ) ( \mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )^{-1}` and the matrix
    elements have to be almost integers.

    The algorithm detail is found at https://atztogo.github.io/niggli/ and
    the references are there in.

    Parameters
    ----------
    lattice: ndarray
        Lattice parameters in the form of

        .. code-block::

            [[a_x, a_y, a_z],
            [b_x, b_y, b_z],
            [c_x, c_y, c_z]]

    eps: float
        Tolerance parameter, but unlike `symprec` the unit is not a length.
        This is used to check if difference of norms of two basis
        vectors is close to zero or not and if two basis vectors are
        orthogonal by the value of dot product being close to zero or
        not.
        The detail is shown at https://atztogo.github.io/niggli/.

    Returns
    -------
    niggli_lattice: ndarray, (3, 3)
        if the Niggli reduction succeeded:
            Reduced lattice parameters are given as a numpy 'double' array:

            .. code-block::

                [[a_x, a_y, a_z],
                [b_x, b_y, b_z],
                [c_x, c_y, c_z]]

        otherwise None is returned.

    Notes
    -----
    .. versionadded:: 1.9.4

    """  # noqa: E501
    _set_no_error()

    niggli_lattice = np.array(np.transpose(lattice), dtype="double", order="C")
    try:
        _spglib.niggli_reduce(niggli_lattice, float(eps))
    except Exception as exc:
        _set_or_throw_error(exc)
        return None
    return np.array(np.transpose(niggli_lattice), dtype="double", order="C")
