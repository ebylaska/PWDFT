#ifndef PWDFT_ECKART_PROJECTOR_HPP
#define PWDFT_ECKART_PROJECTOR_HPP

#include <vector>

namespace pwdft {

/**
 * Build translation-only constraint vectors (periodic solids).
 *
 * N       : number of atoms
 * masses  : atomic masses (amu), size N
 * V       : output constraint matrix (ndof × 3), column-major
 *           V[:,0]=Tx, V[:,1]=Ty, V[:,2]=Tz
 *
 * Note: vectors are mass-weighted and orthonormal (within numerical eps).
 */
void build_translation_constraints(
    int N,
    const std::vector<double>& masses,
    std::vector<double>& V);

/**
 * Build molecular Eckart constraint vectors (translations + rotations).
 *
 * N       : number of atoms
 * coords  : Cartesian coordinates (bohr), size 3N
 * masses  : atomic masses (amu), size N
 * V       : output constraint matrix (ndof × 6), column-major
 * m_out   : number of independent constraints returned (5 or 6)
 *
 * Notes:
 *  - coordinates are internally shifted to COM before building rotations
 *  - vectors are mass-weighted and orthonormalized by Gram-Schmidt
 */
void build_molecular_constraints(
    int N,
    const std::vector<double>& coords,
    const std::vector<double>& masses,
    std::vector<double>& V,
    int& m_out);

/**
 * Apply Eckart projector to Hessian:
 *   H <- (I - P) H (I - P),   where P = V V^T and V has orthonormal columns.
 *
 * ndof : dimension of H
 * m    : number of constraint vectors (columns of V)
 * V    : constraint matrix (ndof × m), column-major
 * H    : Hessian (ndof × ndof), column-major (in/out)
 */
void apply_projector(
    int ndof,
    int m,
    std::vector<double>& V,
    std::vector<double>& H);

} // namespace pwdft

#endif
