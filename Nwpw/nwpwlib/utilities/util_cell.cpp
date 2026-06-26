
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include "units.hpp"

/** @ingroup nwpw_utilities */
namespace pwdft {


// C++ equivalent of Fortran:
//   subroutine cell_unita_abc_abg(unita,lattice)
//   real*8 unita(3,3), lattice(6)
//
// Conventions mirrored exactly:
// - unita(i,j): jth column is jth lattice vector; i is Cartesian component
// - unita is provided in Fortran column-major layout in a flat array unita_f[9]:
//     unita_f[(i-1) + 3*(j-1)] <-> unita(i,j)
// - lattice(1:6) = (a,b,c, alpha,beta,gamma) in radians

static inline double clamp_cos(double x)
{
  // Protect dacos against slight roundoff putting x outside [-1,1]
  return std::max(-1.0, std::min(1.0, x));
}

/*****************************************
 *                                       *
 *        util_cell_unita_abc_abg        *
 *                                       *
 *****************************************/
void util_cell_unita_abc_abg(const double unita_f[9], double lattice[6])
{
   auto U = [&](int i, int j) -> double { return unita_f[(i-1) + 3*(j-1)]; };
 
   // a = sqrt( unita(1,1)^2 + unita(2,1)^2 + unita(3,1)^2 )
   const double a = std::sqrt(U(1,1)*U(1,1) + U(2,1)*U(2,1) + U(3,1)*U(3,1));
 
   // b = sqrt( unita(1,2)^2 + unita(2,2)^2 + unita(3,2)^2 )
   const double b = std::sqrt(U(1,2)*U(1,2) + U(2,2)*U(2,2) + U(3,2)*U(3,2));
 
   // c = sqrt( unita(1,3)^2 + unita(2,3)^2 + unita(3,3)^2 )
   const double c = std::sqrt(U(1,3)*U(1,3) + U(2,3)*U(2,3) + U(3,3)*U(3,3));
 
   // d2 = |unita(:,2) - unita(:,3)|^2
   double d2 = (U(1,2)-U(1,3))*(U(1,2)-U(1,3)) +
               (U(2,2)-U(2,3))*(U(2,2)-U(2,3)) +
               (U(3,2)-U(3,3))*(U(3,2)-U(3,3));
   double alpha = (b*b + c*c - d2)/(2.0*b*c);
   alpha = std::acos(clamp_cos(alpha));
 
   // d2 = |unita(:,3) - unita(:,1)|^2
   d2 = (U(1,3)-U(1,1))*(U(1,3)-U(1,1)) +
        (U(2,3)-U(2,1))*(U(2,3)-U(2,1)) +
        (U(3,3)-U(3,1))*(U(3,3)-U(3,1));
   double beta = (c*c + a*a - d2)/(2.0*c*a);
   beta = std::acos(clamp_cos(beta));
 
   // d2 = |unita(:,1) - unita(:,2)|^2
   d2 = (U(1,1)-U(1,2))*(U(1,1)-U(1,2)) +
        (U(2,1)-U(2,2))*(U(2,1)-U(2,2)) +
        (U(3,1)-U(3,2))*(U(3,1)-U(3,2));
   double gamma = (a*a + b*b - d2)/(2.0*a*b);
   gamma = std::acos(clamp_cos(gamma));
 
   // lattice(1:6) = a,b,c,alpha,beta,gamma
   lattice[0] = a;
   lattice[1] = b;
   lattice[2] = c;
   lattice[3] = alpha;
   lattice[4] = beta;
   lattice[5] = gamma;
}


/*****************************************
 *                                       *
 *        util_cell_abc_abg_unita        *
 *                                       *
 *****************************************/
/*  
   The Jth column of this matrix contains the jth lattice vector
   and    aij is the ith Cartesian component of the jth lattice vector.
*/
// C++ equivalent of Fortran subroutine:
//   cell_abc_abg_unita(lattice, lattice_a)
//
// Fortran assumptions (exactly mirrored here):
//   lattice(1:3) = a,b,c   (in atomic units)
//   lattice(4:6) = alpha,beta,gamma (in radians)
//   lattice_a(3,3): "Jth column contains jth lattice vector"
//   lattice_a(i,j) = ith Cartesian component of jth lattice vector
//
// Memory/layout note:
//   Fortran lattice_a(3,3) is column-major.
//   We store into a flat C++ array lattice_a_f[9] in *Fortran column-major order*:
//     lattice_a_f[(i-1) + 3*(j-1)]  <->  lattice_a(i,j)
//
// So element mapping is:
//   (1,1)->0 (2,1)->1 (3,1)->2  (1,2)->3 (2,2)->4 (3,2)->5  (1,3)->6 (2,3)->7 (3,3)->8

inline double deter3_fortran_gmat(const double gmat_f[9])
{
  // gmat_f is 3x3 in Fortran column-major: gmat(i,j) at (i-1)+3*(j-1)
  const auto G = [&](int i, int j) -> double { return gmat_f[(i-1) + 3*(j-1)]; };

  // Same determinant expansion as the Fortran code conceptually does via deter3(gmat)
  return
      G(1,1) * (G(2,2)*G(3,3) - G(3,2)*G(2,3))
    - G(2,1) * (G(1,2)*G(3,3) - G(3,2)*G(1,3))
    + G(3,1) * (G(1,2)*G(2,3) - G(2,2)*G(1,3));
}


void util_cell_abc_abg_unita(const double lattice[6], double lattice_a_f[9])
{
   // --- Fortran local variables mirrored ---
   double cdist[3], cang[3];
   double gmat_f[9] = {0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0};
 
   // cdist(1:3) = lattice(1:3)
   cdist[0] = lattice[0];
   cdist[1] = lattice[1];
   cdist[2] = lattice[2];
 
   // cang(1:3) = lattice(4:6) = alpha,beta,gamma  (already radians)
   cang[0] = lattice[3]; // alpha
   cang[1] = lattice[4]; // beta
   cang[2] = lattice[5]; // gamma
 
   // Helpers to access gmat in Fortran column-major
   auto GM = [&](int i, int j) -> double& { return gmat_f[(i-1) + 3*(j-1)]; };
 
   // ---- build the metrical matrix gmat exactly like Fortran ----
   // do i=1,3; gmat(i,i)=cdist(i)**2; end do
   for (int i = 1; i <= 3; ++i) {
     GM(i,i) = cdist[i-1] * cdist[i-1];
   }
 
   // iang=3
   // do i=1,3
   //   do j=i+1,3
   //     gmat(i,j)=cdist(i)*cdist(j)*dcos(cang(iang))
   //     gmat(j,i)=gmat(i,j)
   //     iang=iang-1
   //   end do
   // end do
   int iang = 3;
   for (int i = 1; i <= 3; ++i) {
     for (int j = i + 1; j <= 3; ++j) {
       const double val = cdist[i-1] * cdist[j-1] * std::cos(cang[iang-1]);
       GM(i,j) = val;
       GM(j,i) = val;
       --iang;
     }
   }
 
   // vol = dsqrt(deter3(gmat))
   const double det = deter3_fortran_gmat(gmat_f);
   const double vol = std::sqrt(std::max(0.0, det)); // keep non-negative like robust code
 
   // c1=dcos(cang(1)); c2=dcos(cang(2)); c3=dcos(cang(3)); s3=dsin(cang(3))
   const double c1 = std::cos(cang[0]); // cos(alpha)
   const double c2 = std::cos(cang[1]); // cos(beta)
   const double c3 = std::cos(cang[2]); // cos(gamma)
   const double s3 = std::sin(cang[2]); // sin(gamma)
 
   // ---- generate lattice_a exactly like Fortran ----
   // lattice_a(1,1) = cdist(1)*s3
   // lattice_a(1,2) = 0
   // lattice_a(1,3) = cdist(3)*(c2-c1*c3)/s3
   // lattice_a(2,1) = cdist(1)*c3
   // lattice_a(2,2) = cdist(2)
   // lattice_a(2,3) = cdist(3)*c1
   // lattice_a(3,1) = 0
   // lattice_a(3,2) = 0
   // lattice_a(3,3) = vol/(cdist(1)*cdist(2)*s3)
 
   auto A = [&](int i, int j) -> double& { return lattice_a_f[(i-1) + 3*(j-1)]; };
 
   A(1,1) = cdist[0] * s3;
   A(1,2) = 0.0;
   A(1,3) = cdist[2] * (c2 - c1*c3) / s3;
 
   A(2,1) = cdist[0] * c3;
   A(2,2) = cdist[1];
   A(2,3) = cdist[2] * c1;
 
   A(3,1) = 0.0;
   A(3,2) = 0.0;
   A(3,3) = vol / (cdist[0] * cdist[1] * s3);
}



// Accessor for a 3x3 stored in Fortran column-major as flat[9]
static inline double  Aget(const double* A_f, int i, int j) { return A_f[(i-1) + 3*(j-1)]; }
static inline double& Aref(double* A_f, int i, int j)       { return A_f[(i-1) + 3*(j-1)]; }

// --- Fortran invert3(vector, vecinv) ---
// NOTE: your Fortran invert3 declares vector(9) but caller passes vec(3,3).
// This is relying on contiguous storage equivalence. We do the same: invert a 3x3.
static void invert3_fortran(const double vector_f[9], double vecinv_f[9])
{
  // Keep the exact Fortran formulas/indexing (1..9) as written.
  // vector(1)=A(1,1), vector(2)=A(2,1), vector(3)=A(3,1),
  // vector(4)=A(1,2), ..., vector(9)=A(3,3)  (column-major).
  const double* v = vector_f;

  const double deter =
      (v[4]*v[8] - v[5]*v[7]) * v[0] +   // (vector(5)*vector(9)-vector(6)*vector(8))*vector(1)
      (v[5]*v[6] - v[3]*v[8]) * v[1] +   // (vector(6)*vector(7)-vector(4)*vector(9))*vector(2)
      (v[3]*v[7] - v[4]*v[6]) * v[2];    // (vector(4)*vector(8)-vector(5)*vector(7))*vector(3)

  vecinv_f[0] = (v[4]*v[8] - v[5]*v[7]) / deter;
  vecinv_f[1] = (v[7]*v[2] - v[8]*v[1]) / deter;
  vecinv_f[2] = (v[1]*v[5] - v[2]*v[4]) / deter;
  vecinv_f[3] = (v[5]*v[6] - v[3]*v[8]) / deter;
  vecinv_f[4] = (v[8]*v[0] - v[6]*v[2]) / deter;
  vecinv_f[5] = (v[2]*v[3] - v[0]*v[5]) / deter;
  vecinv_f[6] = (v[3]*v[7] - v[4]*v[6]) / deter;
  vecinv_f[7] = (v[6]*v[1] - v[7]*v[0]) / deter;
  vecinv_f[8] = (v[0]*v[4] - v[1]*v[3]) / deter;
}


/*****************************************
 *                                       *
 *        util_cell_lattice_gradient     *
 *                                       *
 *****************************************/
//void util_cell_unita_abc_abg(const double unita_f[9], double lattice[6])
// This is what your Fortran does: amat(ii,jj) = lattice_unita(ii,jj)
// Here we pass in unita_f explicitly (Fortran column-major).
// If in your code lattice_unita is global state, pass that state in.
void util_cell_lattice_gradient(const double stressin_f[9],
                                const double unita_f[9],
                                double lstress[6])
{
   double amat_f[9];
   double stress_f[9];
 
   double vec_f[9], vecinv_f[9];
   double length[3];
 
   // Initialize lstress(1:6)=0
   for (int ii = 1; ii <= 6; ++ii) lstress[ii-1] = 0.0;
 
   // stress(jj,ii)=stressin(jj,ii); zero small
   for (int ii = 1; ii <= 3; ++ii) {
     for (int jj = 1; jj <= 3; ++jj) {
       double v = Aget(stressin_f, jj, ii);
       if (std::abs(v) < 1.0e-10) v = 0.0;
       Aref(stress_f, jj, ii) = v;
     }
   }
 
   // amat(ii,jj)=unita(ii,jj); zero small
   for (int jj = 1; jj <= 3; ++jj) {
     for (int ii = 1; ii <= 3; ++ii) {
       double v = Aget(unita_f, ii, jj);
       if (std::abs(v) < 1.0e-10) v = 0.0;
       Aref(amat_f, ii, jj) = v;
     }
   }
 
   // call lattice_abc_abg(a,b,c,alpha,beta,gamma) ; degrees -> radians
   double lattice[6];
   util_cell_unita_abc_abg(unita_f,lattice);
 
   double a=lattice[0];
   double b=lattice[1];
   double c=lattice[2];
   double alpha=lattice[3];
   double beta=lattice[4];
   double gamma=lattice[5];
 
   alpha = alpha * units::PI / 180.0;
   beta  = beta  * units::PI / 180.0;
   gamma = gamma * units::PI / 180.0;
 
   length[0] = a;
   length[1] = b;
   length[2] = c;
 
   const double c1 = std::cos(alpha);
   const double c2 = std::cos(beta);
   const double c3 = std::cos(gamma);
   const double s1 = std::sin(alpha);
   const double s2 = std::sin(beta);
   const double s3 = std::sin(gamma);
 
   double coser[7] = {0,0,0,0,0,0,0}; // index 4..6 used
   double siner[7] = {0,0,0,0,0,0,0};
   coser[4] = c1; coser[5] = c2; coser[6] = c3;
   siner[4] = s1; siner[5] = s2; siner[6] = s3;
 
   // vol = (a*b*c)*sqrt(1-(c1^2+c2^2+c3^2)+2*c1*c2*c3)
   // (Computed but never used below in your Fortran snippet; kept for fidelity.)
   const double vol = (a*b*c) * std::sqrt( 1.0 - (c1*c1 + c2*c2 + c3*c3) + (2.0*c1*c2*c3));
   (void)vol;
 
   // lstress(1:3): projection for a,b,c
   // do ii=1,3
   //   do jj=1,3
   //     dum = amat(jj,ii) * (stress(jj,ii)/length(ii))
   //     lstress(ii) += dum
   for (int ii = 1; ii <= 3; ++ii) {
     for (int jj = 1; jj <= 3; ++jj) {
       const double dum = Aget(amat_f, jj, ii) * (Aget(stress_f, jj, ii) / length[ii-1]);
       lstress[ii-1] += dum;
     }
   }
 
   // Angle components (4..6)
   for (int ii = 1; ii <= 3; ++ii) {          // a,b,c
     for (int iangle = 4; iangle <= 6; ++iangle) { // alpha,beta,gamma
       if (ii != (iangle - 3)) { // d(a)/d(alpha)=0 etc.
 
         // kk = 9 - (iangle + ii)
         // jj = 6 - (kk + ii)
         const int kk = 9 - (iangle + ii);
         const int jj = 6 - (kk + ii);
 
         // vec(1,ll) = amat(ll,kk)
         // vec(2,ll) = amat(ll,ii)
         // vec(3,ll) = amat(ll,jj)
         for (int ll = 1; ll <= 3; ++ll) {
           Aref(vec_f, 1, ll) = Aget(amat_f, ll, kk);
           Aref(vec_f, 2, ll) = Aget(amat_f, ll, ii);
           Aref(vec_f, 3, ll) = Aget(amat_f, ll, jj);
         }
 
         invert3_fortran(vec_f, vecinv_f);
 
         // do i=1,3
         //   vec(i,1) = vecinv(i,1)*(-length(ii)*length(kk)*siner(iangle))
         // enddo
         for (int i = 1; i <= 3; ++i) {
           Aref(vec_f, i, 1) =
               Aget(vecinv_f, i, 1) * (-length[ii-1] * length[kk-1] * siner[iangle]);
         }
 
         // do i=1,3
         //   lstress(iangle) += vec(i,1)*stress(i,ii)/2
         // enddo
         for (int i = 1; i <= 3; ++i) {
           lstress[iangle-1] += Aget(vec_f, i, 1) * Aget(stress_f, i, ii) / 2.0;
         }
       }
     }
   }
}




}
