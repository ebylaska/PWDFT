#include "eckart_projector.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "blas.h" // *_PWDFT macros

namespace pwdft {

// ============================================================
// Translation-only constraints (solids)
// V is ndof x 3 (column-major), columns are mass-weighted translations
// ============================================================
void build_translation_constraints(int N,
                                  const std::vector<double> &masses,
                                  std::vector<double> &V) {
  int ndof = 3 * N;
  V.assign(ndof * 3, 0.0);

  double M = 0.0;
  for (int a = 0; a < N; ++a)
    M += masses[a];

  // mass-weighted translations: V(3a+d, d) = sqrt(m_a / M)
  for (int a = 0; a < N; ++a) {
    const double w = std::sqrt(masses[a] / M);
    for (int d = 0; d < 3; ++d) {
      const int row = 3 * a + d;
      V[row + d * ndof] = w; // column-major
    }
  }
}

// ============================================================
// Molecular constraints (translations + rotations)
// V is ndof x 6 (column-major). Returns m_out = number of independent
// constraints (usually 6; 5 for linear molecules).
// ============================================================
void build_molecular_constraints(int N, const std::vector<double> &coords,
                                const std::vector<double> &masses,
                                std::vector<double> &V, int &m_out) {
  int ndof = 3 * N;
  V.assign(ndof * 6, 0.0);

  // BLAS-friendly scalars
  int one = 1;
  double rone = 1.0;
  double rzero = 0.0;

  // ---- center of mass
  double com[3] = {0.0, 0.0, 0.0};
  double M = 0.0;
  for (int a = 0; a < N; ++a) {
    M += masses[a];
    com[0] += masses[a] * coords[3 * a + 0];
    com[1] += masses[a] * coords[3 * a + 1];
    com[2] += masses[a] * coords[3 * a + 2];
  }
  com[0] /= M;
  com[1] /= M;
  com[2] /= M;

  // ---- translations (mass-weighted; orthonormal by construction)
  for (int a = 0; a < N; ++a) {
    const double w = std::sqrt(masses[a] / M);
    for (int d = 0; d < 3; ++d) {
      const int row = 3 * a + d;
      V[row + d * ndof] = w;
    }
  }

  // ---- rotations about x,y,z through COM:
  // v = sqrt(m_a) * (r_a - com) × e_axis
  for (int axis = 0; axis < 3; ++axis) {
    const int col = 3 + axis;
    for (int a = 0; a < N; ++a) {
      const double x = coords[3 * a + 0] - com[0];
      const double y = coords[3 * a + 1] - com[1];
      const double z = coords[3 * a + 2] - com[2];

      double wx = 0.0, wy = 0.0, wz = 0.0;
      if (axis == 0) {
        wy = -z;
        wz = y;
      } // about x
      if (axis == 1) {
        wx = z;
        wz = -x;
      } // about y
      if (axis == 2) {
        wx = -y;
        wy = x;
      } // about z

      const double w = std::sqrt(masses[a]);
      V[(3 * a + 0) + col * ndof] = w * wx;
      V[(3 * a + 1) + col * ndof] = w * wy;
      V[(3 * a + 2) + col * ndof] = w * wz;
    }
  }

  // ---- Gram-Schmidt orthonormalize columns of V (6 columns),
  //      drop dependent ones (linear molecule -> 5 constraints)
  int m = 6;

  for (int i = 0; i < 6; ++i) {
    // subtract projections on previous vectors
    for (int j = 0; j < i; ++j) {
      const double dot =
          DDOT_PWDFT(ndof, &V[j * ndof], one, &V[i * ndof], one);
      double alpha = -dot;
      DAXPY_PWDFT(ndof, alpha, &V[j * ndof], one, &V[i * ndof],
                  one); // Vi -= dot*Vj
    }

    const double norm2 =
        DDOT_PWDFT(ndof, &V[i * ndof], one, &V[i * ndof], one);

    if (norm2 < 1.0e-10) {
      std::fill(&V[i * ndof], &V[(i + 1) * ndof], 0.0);
      --m; // dependent constraint vector
    } else {
      double inv = rone / std::sqrt(norm2);
      DSCAL_PWDFT(ndof, inv, &V[i * ndof], one);
    }
  }

  // We keep ndof×6 storage; caller uses first m columns (m_out).
  m_out = m;
}

// ============================================================
// Apply projector: H <- (I-P) H (I-P), P = V V^T
// V: ndof×m (stored in a larger array is OK; we use first m columns)
// H: ndof×ndof (column-major)
// ============================================================
void apply_projector(int ndof, int m, 
                     std::vector<double> &V,
                     std::vector<double> &H) {
  if (m <= 0)
    return;

  int one = 1; // for daxpy increments if needed
  double rone = 1.0;
  double rzero = 0.0;
  double mrone = -1.0;

  std::vector<double> P(ndof * ndof, 0.0);
  std::vector<double> W1(ndof * ndof, 0.0);
  std::vector<double> W2(ndof * ndof, 0.0);

  // P = V V^T  (V is ndof×m)
  DGEMM_PWDFT((char *)"N", (char *)"T", ndof, ndof, m, rone, V.data(), ndof,
             V.data(), ndof, rzero, P.data(), ndof);

  // W1 = P H
  DGEMM_PWDFT((char *)"N", (char *)"N", ndof, ndof, ndof, rone, P.data(), ndof,
             H.data(), ndof, rzero, W1.data(), ndof);

  // W2 = H P
  DGEMM_PWDFT((char *)"N", (char *)"N", ndof, ndof, ndof, rone, H.data(), ndof,
             P.data(), ndof, rzero, W2.data(), ndof);

  // H = H - W1 - W2
  {
    int n = ndof * ndof;
    int one = 1;
    DAXPY_PWDFT(n, mrone, W1.data(), one, H.data(), one);
    DAXPY_PWDFT(n, mrone, W2.data(), one, H.data(), one);
  }

  // W1 = P W2 = P H P
  DGEMM_PWDFT((char *)"N", (char *)"N", ndof, ndof, ndof, rone, P.data(), ndof,
             W2.data(), ndof, rzero, W1.data(), ndof);

  // H = H + W1
  {
    int n = ndof * ndof;
    int one = 1;
    DAXPY_PWDFT(n, rone, W1.data(), one, H.data(), one);
  }
}

} // namespace pwdft
