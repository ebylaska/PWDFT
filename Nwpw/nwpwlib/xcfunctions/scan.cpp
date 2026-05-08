/* scan.cpp
  Author - Eric Bylaska
*/

#include <cmath>

namespace pwdft {


static constexpr double pi = 3.14159265358979311599;
static constexpr double thrd = (1.0 / 3.0);
static constexpr double twthrd = (2.0 / 3.0);
static constexpr double frthrd = (4.0 / 3.0);
static constexpr double fvthrd = (5.0 / 3.0);
static constexpr double etthrd = (8.0 / 3.0);
static constexpr double tol = 1.0e-18;
static constexpr double ETA = 1.0e-20;
static constexpr double a1 = 4.9479;
static constexpr double b3 = 0.5;
static constexpr double c1x = 0.667;
static constexpr double c2x = 0.8;
static constexpr double dx_x = 1.24;
static constexpr double muAK = (10.0 / 81.0);
static constexpr double K1 = 0.065;
static constexpr double h0x = 1.174;
static constexpr double gamma_scan = 0.03109069086965489503494086371273;
static constexpr double beta0 = 0.06672455060314922;
static constexpr double beta1 = 0.1;
static constexpr double beta2 = 0.1778;
static constexpr double b1c = 0.0285764;
static constexpr double b2c = 0.0889;
static constexpr double b3c = 0.125541;
static constexpr double c1c = 0.64;
static constexpr double c2c = 1.5;
static constexpr double dc = 0.7;
static constexpr double dxc = 2.3631;
static constexpr double xi = 0.12802585262625815;
static constexpr double thr1 = 0.996;
static constexpr double thr2 = 1.004;

/************************************************
 *                                              *
 *                pw91c_zeta                    *
 *                                              *
 ************************************************/
/**
 * @brief Computes the $\zeta$-dependent contribution to the PW91 exchange functional.
 * 
 * This function calculates the spin-polarization dependent part of the 
 * Perdew-Wang (PW91) exchange energy density. It uses a numerical tolerance 
 * to ensure stability in the limit where the spin-polarization $\zeta$ 
 * approaches zero, preventing singularities.
 *
 * @param tol   Numerical tolerance threshold to maintain stability near $\zeta = 0$.
 * @param s     Scaling parameter related to the local density.
 * @param zeta  The spin-polarating parameter (representing spin-up/spin-down ratio).
 * 
 * @note The function relies on external state/pointers to update the functional density.
 */
static void pw91c_zeta(double tol, double s, double zeta, double *fz, double *dfdz) 
{
   /* local variables */
   double small;
   double omz;
   double opz;
   double omz_13;
   double opz_13;
 
   small = tol;
   *fz = -2.0;
   *dfdz = 0.0;
 
   omz = 1.0 - zeta;
   opz = 1.0 + zeta;
 
   if (omz > small) 
   {
      //omz_13 = std::pow(omz, thrd);
      omz_13 = std::cbrt(omz);

      *fz = *fz + omz * omz_13;
      *dfdz = *dfdz - omz_13;
   }
 
   if (opz > small) 
   {
      //opz_13 = std::pow(opz, thrd);
      opz_13 = std::cbrt(opz);
      *fz = *fz + opz * opz_13;
      *dfdz = *dfdz + opz_13;
   }
 
   *fz = *fz * s;
   *dfdz = frthrd * (*dfdz) * s;
}


/************************************************
 *                                              *
 *                pw91c_rs                      *
 *                                              *
 ************************************************/
/**
 * @brief Computes the PW91 exchange energy density and its $r_s$ derivative.
 * 
 * This function implements the core of the PW91 exchange functional expansion 
 * with respect to the Wigner-Seitz radius ($r_s$). It utilizes the $q$-parameter 
 * expansion method to evaluate the energy density $v(r_s)$ and its first 
 * derivative $\frac{dv}{dr_s}$ through intermediate $q_0, q_1,$ and $q_2$ 
 * polynomials.
 * 
 * @param a, a1, b1, b2, b3, b4  The fundamental PW91 functional coefficients.
 * @param rs                     The Wigner-Seitz radius (dimensionless or scaled).
 * @param v                      [out] The computed exchange energy density.
 * @param dvdr                   [out] The derivative of the density w.r.t $r_s$.
 */
static void pw91c_rs(double a, double a1, double b1, double b2, double b3, double b4, double rs, double *v, double *dvdr) 
{
   /* local variables */
   double q0;
   double rs_12;
   double rs_32;
   double q1;
   double q2;
   double dq0_drs;
   double dq1_drs;
   double dq2_drs;
 
   q0 = -2.0 * a * (1.0 + a1 * rs);
 
   rs_12 = std::sqrt(rs);
   rs_32 = rs * rs_12;
 
   q1 = 2.0 * a * (b1 * rs_12 + b2 * rs + b3 * rs_32 + b4 * rs * rs);
   q2 = std::log(1.0 + 1.0 / q1);
 
   *v = q0 * q2;
 
   dq0_drs = -2.0 * a * a1;
   dq1_drs = a * (b1 / rs_12 + 2.0 * b2 + 3.0 * b3 * rs_12 + 4.0 * b4 * rs);
   dq2_drs = -dq1_drs / (q1 + q1 * q1);
 
   *dvdr = dq0_drs * q2 + q0 * dq2_drs;
}


/************************************************
 *                                              *
 *                gen_PW91_c_rz                 *
 *                                              *
 ************************************************/
/**
 * @brief Computes the PW91 correlation energy density and its derivatives.
 * 
 * This function serves as the top-level assembly routine for the PW91 
 * correlation functional. It integrates results from the spin-polarization 
 * kernel (pw91c_zeta) and the Wigner-Seitz radius kernels (pw91c_rs) 
 * to compute the final energy density and its partial derivatives.
 *
 * @param tol         Numerical stability tolerance for $\zeta \to 0$ limit.
 * @param rs          The Wigner-Seitz radius ($r_s$).
 * @param zeta        The spin-polarization parameter ($\zeta$).
 * @param pwc         [out] Pointer to the computed correlation energy density.
 * @param dpwc_drs    [out] Pointer to the derivative $\partial (E_c) / \partial r_s$.
 * @param dpwc_dzeta  [out] Pointer to the derivative $\partial (E_c) / \partial \zeta$.
 * 
 * @note This function relies on the correctness of the underlying 
 *       $\zeta$ and $r_s$ expansion kernels.
 */
static void gen_PW91_c_rz(double tol, double rs, double zeta, double *pwc, double *dpwc_drs, double *dpwc_dzeta) 
{
   /* local variables */
   double zeta2;
   double zeta3;
   double zeta4;
   double eu;
   double deu_drs;
   double ep;
   double dep_drs;
   double alpham;
   double dam_drs;
   double gz;
   double hz;
   double dgz;
   double dhz;
   double fzeta;
   double df_dzeta;
 
   /* constants from the original Fortran */
   static constexpr double eps0c1 = 0.03109070;
   static constexpr double eps0c2 = 0.21370;
   static constexpr double eps0c3 = 7.5957;
   static constexpr double eps0c4 = 3.5876;
   static constexpr double eps0c5 = 1.6382;
   static constexpr double eps0c6 = 0.49294;
 
   static constexpr double eps1c1 = 0.01554535;
   static constexpr double eps1c2 = 0.20548;
   static constexpr double eps1c3 = 14.1189;
   static constexpr double eps1c4 = 6.1977;
   static constexpr double eps1c5 = 3.3662;
   static constexpr double eps1c6 = 0.62517;
 
   static constexpr double epsc1 = 0.01688686394;
   static constexpr double epsc2 = 0.11125;
   static constexpr double epsc3 = 10.3570;
   static constexpr double epsc4 = 3.6231;
   static constexpr double epsc5 = 0.88026;
   static constexpr double epsc6 = 0.49671;
 
   double fzzI = 9.0 * (std::pow(2.0, thrd) - 1.0) / 4.0;
   double gammaI = 1.0 / (2.0 * std::pow(2.0, thrd) - 2.0);
 
   pw91c_zeta(tol, gammaI, zeta, &fzeta, &df_dzeta);
 
   pw91c_rs(eps0c1, eps0c2, eps0c3, eps0c4, eps0c5, eps0c6, rs, &eu, &deu_drs);
 
   pw91c_rs(eps1c1, eps1c2, eps1c3, eps1c4, eps1c5, eps1c6, rs, &ep, &dep_drs);
 
   pw91c_rs(epsc1, epsc2, epsc3, epsc4, epsc5, epsc6, rs, &alpham, &dam_drs);
 
   zeta2 = zeta * zeta;
   zeta3 = zeta2 * zeta;
   zeta4 = zeta3 * zeta;
 
   gz = fzeta * zeta4;
   hz = fzzI * (fzeta - gz);
   dgz = df_dzeta * zeta4 + 4.0 * fzeta * zeta3;
   dhz = fzzI * (df_dzeta - dgz);
 
   *pwc = eu * (1.0 - gz) + ep * gz - alpham * hz;
   *dpwc_drs = deu_drs * (1.0 - gz) + dep_drs * gz - dam_drs * hz;
   *dpwc_dzeta = (ep - eu) * dgz - alpham * dhz;
}






/****************************************
 *					*
 *	    gen_SCAN_unrestricted	*
 *					*
 ****************************************/
/*
     This function returns the SCAN exchange-correlation
   energy density, xce, and its derivatives with respect
   to nup, ndn, |grad nup|, |grad ndn|, and |grad n|.

    Entry - n2ft3d     : number of grid points
            dn_in(*,2) : spin densites nup and ndn
            agr_in(*,3): |grad nup|, |grad ndn|, and |grad n|
            tau_in(*,2): |grad nup|, |grad ndn|, and |grad n|
            x_parameter: scale parameter for exchange
            c_parameter: scale parameter for correlation

    Exit - xce(*)  : PBE96 energy density
         - fn(*,2) : d(n*xce)/dnup, d(n*xce)/dndn
         - fdn(*,3): d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|
                     d(n*xce)/d|grad n|
           fdtau(*,2)
*/
void gen_SCAN_unrestricted(const int n2ft3d, double *dn_in, double *agr_in, double *tau_in,
                           const double x_parameter, const double c_parameter, 
                           double *xce, double *fn, double *fdn, double *fdtau) 
{
}

/****************************************
 *					*
 *	    gen_SCAN_restricted	        *
 *					*
 ****************************************/
/*
   This routine calculates the SCAN exchange-correlation
   potential(xcp) and energy density(xce).

   Entry - n2ft3d     : number of grid points
           rho_in(*) :  density (nup+ndn)
           agr_in(*): |grad rho_in|
           tau_in(*):  tau
           x_parameter: scale parameter for exchange
           c_parameter: scale parameter for correlation

     Exit  - xce(n2ft3d) : PBE96 exchange correlation energy density
             fn(n2ft3d)  : d(n*xce)/dn
             fdn(n2ft3d) : d(n*xce/d|grad n|
             fdtau()
*/

void gen_SCAN_restricted(const int n2ft3d, double *rho_in, double *agr_in, double *tau_in,
                         const double x_parameter, const double c_parameter,
                         double *xce, double *fn, double *fdn, double *fdtau) 
{
}


} // namespace pwdft
