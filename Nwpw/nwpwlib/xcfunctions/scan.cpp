/* scan.cpp
  Author - Eric Bylaska
*/

#include <cmath>
#include <cstddef>
#include <algorithm>

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
static constexpr double gamma = 0.03109069086965489503494086371273;
static constexpr double dx    = 1.24;
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
void gen_SCAN_unrestricted(const int n2ft3d, double *rho_in, double *agr_in, double *tau_in,
                           const double x_parameter, const double c_parameter, 
                           double *xce, double *fn, double *fdn, double *fdtau) 
{

    // SCAN constants
    double b1, b2, b4;
    constexpr double delta = 1.0e-3;

    b2 = std::sqrt(5913.0/405000.0);
    b1 = (511.0/13500.0)/(2.0*b2);
    b4 = muAK*muAK/K1 - 1606.0/18225.0 - b1*b1;

    const double Cx   = (-0.75)*std::pow(3.0/pi, thrd);
    const double P23  = std::pow(3.0*pi*pi, twthrd);
    const double P13  = std::pow(3.0/(4.0*pi), thrd);
    const double P23t = std::pow(3.0*pi*pi/16.0, twthrd);

    for (int i = 0; i < n2ft3d; ++i) 
    {
        // inputs + ETA
        double nup   = rho_in[i] + ETA;
        double agrup = agr_in[i] + ETA;
        double tauup = tau_in[i] + ETA;

        double ndn   = rho_in[n2ft3d + i] + ETA;
        double agrdn = agr_in[n2ft3d + i] + ETA;
        double taudn = tau_in[n2ft3d + i] + ETA;

        // ========= SCAN Exchange: UP =========
        double n   = 2.0*nup;
        double agr = 2.0*agrup;
        double tau = 2.0*tauup;

        double n_13 = std::pow(n, thrd);
        double n_53 = n_13*n_13*n;
        double n_83 = n_53*n;
        double agr2 = agr*agr;

        double p       = agr2/(4.0*P23*n_83);
        double p_14    = std::sqrt(std::sqrt(p));
        double dp_dn   = -etthrd*p/n;
        double dp_dagr =  2.0*p/agr;

        double tauW = 0.125*agr2/n;
        double tauU = 0.3*P23*n_53;
        double z    = tauW/tau;
        double dz_dn   = -z/n;
        double dz_dagr =  2.0*z/agr;
        double dz_dtau = -z/tau;
        (void)dz_dn; (void)dz_dagr; (void)dz_dtau; // not used later in this exchange part
        double z2 = z*z; (void)z2;

        double alpha = (tau - tauW)/tauU;

        double dalpha_dn=0.0, dalpha_dagr=0.0, dalpha_dtau=0.0;
        double thr0 = agr2/(8.0*tau);
        double dthr0_dagr = 2.0*thr0/agr;
        double dthr0_dtau = -thr0/tau;
        (void)thr0; (void)dthr0_dagr; (void)dthr0_dtau; // only used in commented Fortran branch
        double tauUn = tauU*n;

        if (alpha < 0.0) {
            alpha = 0.0;
            dalpha_dn = dalpha_dagr = dalpha_dtau = 0.0;
        } else {
            double falpha       = tau*n - agr2/8.0;
            double dfalpha_dn   = tau;
            double dfalpha_dagr = -agr/4.0;
            double dfalpha_dtau = n;

            alpha = falpha/tauUn;
            dalpha_dn   = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn;
            dalpha_dagr = dfalpha_dagr/tauUn;
            dalpha_dtau = dfalpha_dtau/tauUn;
        }

        double oma  = 1.0 - alpha;
        double oma2 = oma*oma;

        double exp1 = std::exp(-b4*p/muAK);
        double exp2 = std::exp(-b3*oma2);

        double x1 = muAK*p*(1.0 + (b4*p/muAK)*exp1);
        double x2 = b1*p + b2*oma*exp2;
        double x  = x1 + x2*x2;

        double denh1x = K1 + x;
        double numh1x = denh1x + K1*x;
        double h1x    = numh1x/denh1x;

        double dx1_dp    = muAK + b4*p*exp1*(2.0 - p*b4/muAK);
        double dx2_dp    = b1;
        double dx_dp     = dx1_dp + 2.0*x2*dx2_dp;
        double dx_dalpha = 2.0*b2*exp2*x2*(2.0*b3*oma2 - 1.0);

        double dh1x_dx     = std::pow(K1/denh1x, 2.0);
        double dh1x_dp     = dh1x_dx*dx_dp;
        double dh1x_dalpha = dh1x_dx*dx_dalpha;

        double gx = 1.0, dgx_dp = 0.0;
        if (p_14 < 0.002) {
            gx = 1.0; dgx_dp = 0.0;
        } else {
            double exp3 = std::exp(-a1/p_14);
            gx     = 1.0 - exp3;
            dgx_dp = -0.25*a1*exp3/(p*p_14);
        }

        double fxa = 0.0, dfxa_dalpha = 0.0;
        if (alpha <= thr1) {
            double exp4 = std::exp(-c1x*alpha/oma);
            fxa = exp4;
            dfxa_dalpha = -c1x*exp4/oma2;
        } else if (alpha < thr2) {
            fxa = 0.0;
            dfxa_dalpha = 0.0;
        } else {
            double exp5 = std::exp(c2x/oma);
            fxa = -dx*exp5;
            dfxa_dalpha = -dx*c2x*exp5/oma2;
        }

        double Fx = (h1x + fxa*(h0x - h1x))*gx;

        double dFx_dp =
            dgx_dp*(h1x + fxa*(h0x - h1x)) +
            gx*dh1x_dp*(1.0 - fxa);

        double dFx_dalpha =
            gx*(dh1x_dalpha + dfxa_dalpha*(h0x - h1x) - fxa*dh1x_dalpha);

        double dFx_dn   = dFx_dalpha*dalpha_dn   + dFx_dp*dp_dn;
        double dFx_dagr = dFx_dalpha*dalpha_dagr + dFx_dp*dp_dagr;
        double dFx_dtau = dFx_dalpha*dalpha_dtau;

        double ex0  = Cx*n_13;
        double nex0 = n*ex0;

        double eupx     = ex0*Fx;
        double fnupx    = nex0*dFx_dn + frthrd*eupx;
        double fdnupx   = nex0*dFx_dagr;
        double fdtauupx = nex0*dFx_dtau;

        // ========= SCAN Exchange: DOWN =========
        n   = 2.0*ndn;
        agr = 2.0*agrdn;
        tau = 2.0*taudn;

        n_13 = std::pow(n, thrd);
        n_53 = n_13*n_13*n;
        n_83 = n_53*n;
        agr2 = agr*agr;

        p       = agr2/(4.0*P23*n_83);
        p_14    = std::sqrt(std::sqrt(p));
        dp_dn   = -etthrd*p/n;
        dp_dagr =  2.0*p/agr;

        tauW = 0.125*agr2/n;
        tauU = 0.3*P23*n_53;
        z    = tauW/tau;
        z2   = z*z; (void)z2;

        alpha = (tau - tauW)/tauU;

        dalpha_dn=0.0; dalpha_dagr=0.0; dalpha_dtau=0.0;
        thr0 = agr2/(8.0*tau);
        dthr0_dagr = 2.0*thr0/agr;
        dthr0_dtau = -thr0/tau;
        (void)thr0; (void)dthr0_dagr; (void)dthr0_dtau;

        tauUn = tauU*n;

        if (alpha < 0.0) {
            alpha = 0.0;
            dalpha_dn = dalpha_dagr = dalpha_dtau = 0.0;
        } else {
            double falpha       = tau*n - agr2/8.0;
            double dfalpha_dn   = tau;
            double dfalpha_dagr = -agr/4.0;
            double dfalpha_dtau = n;

            alpha = falpha/tauUn;
            dalpha_dn   = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn;
            dalpha_dagr = dfalpha_dagr/tauUn;
            dalpha_dtau = dfalpha_dtau/tauUn;
        }

        oma  = 1.0 - alpha;
        oma2 = oma*oma;

        exp1 = std::exp(-b4*p/muAK);
        exp2 = std::exp(-b3*oma2);

        x1 = muAK*p*(1.0 + (b4*p/muAK)*exp1);
        x2 = b1*p + b2*oma*exp2;
        x  = x1 + x2*x2;

        denh1x = K1 + x;
        numh1x = denh1x + K1*x;
        h1x    = numh1x/denh1x;

        dx1_dp    = muAK + b4*p*exp1*(2.0 - p*b4/muAK);
        dx2_dp    = b1;
        dx_dp     = dx1_dp + 2.0*x2*dx2_dp;
        dx_dalpha = 2.0*b2*exp2*x2*(2.0*b3*oma2 - 1.0);

        dh1x_dx     = std::pow(K1/denh1x, 2.0);
        dh1x_dp     = dh1x_dx*dx_dp;
        dh1x_dalpha = dh1x_dx*dx_dalpha;

        gx = 1.0; dgx_dp = 0.0;
        if (p_14 < 0.002) {
            gx = 1.0; dgx_dp = 0.0;
        } else {
            double exp3 = std::exp(-a1/p_14);
            gx     = 1.0 - exp3;
            dgx_dp = -0.25*a1*exp3/(p*p_14);
        }

        fxa = 0.0; dfxa_dalpha = 0.0;
        if (alpha <= thr1) {
            double exp4 = std::exp(-c1x*alpha/oma);
            fxa = exp4;
            dfxa_dalpha = -c1x*exp4/oma2;
        } else if (alpha < thr2) {
            fxa = 0.0;
            dfxa_dalpha = 0.0;
        } else {
            double exp5 = std::exp(c2x/oma);
            fxa = -dx*exp5;
            dfxa_dalpha = -dx*c2x*exp5/oma2;
        }

        Fx = (h1x + fxa*(h0x - h1x))*gx;

        dFx_dp =
            dgx_dp*(h1x + fxa*(h0x - h1x)) +
            gx*dh1x_dp*(1.0 - fxa);

        dFx_dalpha =
            gx*(dh1x_dalpha + dfxa_dalpha*(h0x - h1x) - fxa*dh1x_dalpha);

        dFx_dn   = dFx_dalpha*dalpha_dn   + dFx_dp*dp_dn;
        dFx_dagr = dFx_dalpha*dalpha_dagr + dFx_dp*dp_dagr;
        dFx_dtau = dFx_dalpha*dalpha_dtau;

        ex0  = Cx*n_13;
        nex0 = n*ex0;

        double ednx     = ex0*Fx;
        double fndnx    = nex0*dFx_dn + frthrd*ednx;
        double fdndnx   = nex0*dFx_dagr;
        double fdtaudnx = nex0*dFx_dtau;

        // spin-combined exchange energy density
        n = nup + ndn;
        double ex = (eupx*nup + ednx*ndn)/n;

        // ========= SCAN Correlation =========
        agr = agr_in[2*n2ft3d + i] + ETA;
        tau = tauup + taudn;

        n_13 = std::pow(n, thrd);
        n_53 = n*n_13*n_13;
        n_83 = n*n_53;
        double n2 = n*n;
        agr2 = agr*agr;

        double rs     = P13/n_13;
        double drs_dn = -thrd*rs/n;
        double rs_12  = std::sqrt(rs);

        double zeta = (nup - ndn)/n;
        double dzeta_dnup = 0.0, dzeta_dndn = 0.0;
        if (std::abs(nup - ndn) < tol) {
            zeta = 0.0;
            dzeta_dnup = 0.0;
            dzeta_dndn = 0.0;
        } else {
            dzeta_dnup =  2.0*ndn/n2;
            dzeta_dndn = -2.0*nup/n2;
        }

        p       = agr2/(4.0*P23*n_83);
        dp_dn   = -etthrd*p/n;
        dp_dagr =  2.0*p/agr;

        double opz    = 1.0 + zeta;
        double omz    = 1.0 - zeta;
        double opz_23 = std::pow(opz, twthrd);
        double omz_23 = std::pow(omz, twthrd);

        double phi  = 0.5*(opz_23 + omz_23);
        double phi2 = phi*phi;
        double phi3 = phi2*phi;

        double dphi_dzeta = 0.0;
        if (omz < tol) {
            dphi_dzeta =  0.5*twthrd*(opz_23/opz);
        } else if (opz < tol) {
            dphi_dzeta = -0.5*twthrd*(omz_23/omz);
        } else {
            dphi_dzeta = 0.5*twthrd*(opz_23/opz - omz_23/omz);
        }

        double beta      = beta0*(1.0 + beta1*rs)/(1.0 + beta2*rs);
        double dbeta_drs = beta0*(beta1 - beta2)/std::pow(1.0 + beta2*rs, 2.0);

        double ecLDA1, decLDA1_drs, decLDA1_dzeta;
        gen_PW91_c_rz(tol, rs, zeta, &ecLDA1, &decLDA1_drs, &decLDA1_dzeta);

        double w1fac     = ecLDA1/(gamma*phi3);
        double expw1     = std::exp(-w1fac);
        double w1        = expw1 - 1.0;
        double dw1_drs   = -expw1*decLDA1_drs/(gamma*phi3);
        double dw1_dzeta = (3.0*w1fac*dphi_dzeta/phi - decLDA1_dzeta/(gamma*phi3))*expw1;

        double A        = beta/(gamma*w1);
        double dA_drs   = dbeta_drs/(gamma*w1) - A*dw1_drs/w1;
        double dA_dzeta = -A*dw1_dzeta/w1;

        double t2        = P23t*p/(phi2*rs);
        double dt2_drs   = -t2/rs;
        double dt2_dzeta = -2.0*t2*dphi_dzeta/phi;
        double dt2_dp    = P23t/(phi2*rs);

        double At2     = A*t2;
        double dengAt2 = 1.0 + 4.0*At2;
        double gAt2    = 1.0/(std::sqrt(std::sqrt(dengAt2)));

        double tmp1        = 1.0 + w1*(1.0 - gAt2);
        double dtmp1_drs   = dw1_drs*(1.0 - gAt2) + gAt2*w1*(t2*dA_drs + A*dt2_drs)/dengAt2;
        double dtmp1_dzeta = dw1_dzeta*(1.0 - gAt2) + gAt2*w1*(t2*dA_dzeta + A*dt2_dzeta)/dengAt2;
        double dtmp1_dp    = gAt2*w1*A*dt2_dp/dengAt2;

        double H1        = gamma*phi3*std::log(tmp1);
        double dH1_drs   = gamma*phi3*dtmp1_drs/tmp1;
        double dH1_dzeta = 3.0*H1*dphi_dzeta/phi + gamma*phi3*dtmp1_dzeta/tmp1;
        double dH1_dp    = gamma*phi3*dtmp1_dp/tmp1;

        double ec1        = ecLDA1 + H1;
        double dec1_drs   = decLDA1_drs + dH1_drs;
        double dec1_dzeta = decLDA1_dzeta + dH1_dzeta;
        double dec1_dp    = dH1_dp;

        double ecLDA0      = -b1c/(1.0 + b2c*rs_12 + b3c*rs);
        double decLDA0_drs = b1c*(b3c + 0.5*b2c/rs_12)/std::pow(1.0 + b2c*rs_12 + b3c*rs, 2.0);

        double dxz        = 0.5*(std::pow(opz, frthrd) + std::pow(omz, frthrd));
        double ddxz_dzeta = 0.5*frthrd*(std::pow(opz, thrd) - std::pow(omz, thrd));
        double zeta12     = std::pow(zeta, 12.0);
        double omz12      = 1.0 - zeta12;
        double zeta11     = std::pow(zeta, 11.0);

        double gc        = (1.0 - dxc*(dxz - 1.0))*omz12;
        double dgc_dzeta = -dxc*ddxz_dzeta*omz12 - 12.0*zeta11*(1.0 - dxc*(dxz - 1.0));

        double w0fac   = ecLDA0/b1c;
        double expw0   = std::exp(-w0fac);
        double w0      = expw0 - 1.0;
        double dw0_drs = -decLDA0_drs*expw0/b1c;

        double ginf     = 1.0/(std::sqrt(std::sqrt(1.0 + 4.0*xi*p)));
        double dginf_dp = -xi*ginf/(1.0 + 4.0*xi*p);

        double tmp0      = 1.0 + w0*(1.0 - ginf);
        double dtmp0_drs = dw0_drs*(1.0 - ginf);
        double dtmp0_dp  = -w0*dginf_dp;

        double H0      = b1c*std::log(tmp0);
        double dH0_drs = b1c*dtmp0_drs/tmp0;
        double dH0_dp  = b1c*dtmp0_dp/tmp0;

        double ec0        = (ecLDA0 + H0)*gc;
        double dec0_drs   = gc*(decLDA0_drs + dH0_drs);
        double dec0_dzeta = dgc_dzeta*(ecLDA0 + H0);
        double dec0_dp    = gc*dH0_dp;

        double ds        = 0.5*(std::pow(opz, fvthrd) + std::pow(omz, fvthrd));
        double dds_dzeta = 0.5*fvthrd*(std::pow(opz, twthrd) - std::pow(omz, twthrd));

        tauW = 0.125*agr2/n;
        tauU = 0.3*P23*ds*n_53;
        z    = tauW/tau;
        dz_dn   = -z/n;
        dz_dagr =  2.0*z/agr;
        dz_dtau = -z/tau;
        (void)dz_dn; (void)dz_dagr; (void)dz_dtau;
        z2 = z*z; (void)z2;

        // alpha for correlation: uses smooth switching (the non-commented Fortran branch)
        double dalpha_dnup=0.0; 
        double dalpha_dndn=0.0; 
        dalpha_dagr=0.0; 
        double dalpha_dtauup=0.0; 
        double dalpha_dtaudn=0.0;

        thr0 = agr2/(8.0*tau);
        dthr0_dagr = 2.0*thr0/agr;
        dthr0_dtau = -thr0/tau;

        double delta2 = delta*delta;
        double diff  = n - thr0;
        double diff2 = diff*diff;
        double diff3 = diff2*diff;
        double dod   = diff/delta;
        double dod2  = dod*dod;

        tauUn = tauU*n;

        if (n <= thr0) {
            alpha = 0.0;
        } else if (n < (thr0 + delta)) {
            double falpha       = tau*(-diff3/delta2 + 2.0*diff2/delta);
            double dfalpha_dn   = tau*(-3.0*dod2 + 4.0*dod);
            double dfalpha_dagr = tau*(3.0*dod2 - 4.0*dod)*dthr0_dagr;
            double dfalpha_dtau = falpha/tau + tau*(3.0*dod2 - 4.0*dod)*dthr0_dtau;

            alpha = falpha/tauUn;
            dalpha_dnup = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn - (alpha/ds)*dds_dzeta*dzeta_dnup;
            dalpha_dndn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn - (alpha/ds)*dds_dzeta*dzeta_dndn;
            dalpha_dagr = dfalpha_dagr/tauUn;
            dalpha_dtauup = dfalpha_dtau/tauUn;
            dalpha_dtaudn = dalpha_dtauup;
        } else { // n >= thr0 + delta
            double falpha       = tau*n - agr2/8.0;
            double dfalpha_dn   = tau;
            double dfalpha_dagr = -agr/4.0;
            double dfalpha_dtau = n;

            alpha = falpha/tauUn;
            dalpha_dnup = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn - (alpha/ds)*dds_dzeta*dzeta_dnup;
            dalpha_dndn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn - (alpha/ds)*dds_dzeta*dzeta_dndn;
            dalpha_dagr = dfalpha_dagr/tauUn;
            dalpha_dtauup = dfalpha_dtau/tauUn;
            dalpha_dtaudn = dalpha_dtauup;
        }

        oma  = 1.0 - alpha;
        oma2 = oma*oma;

        double fca = 0.0, dfca_dalpha = 0.0;
        if (alpha <= thr1) {
            double exp6 = std::exp(-c1c*alpha/oma);
            fca = exp6;
            dfca_dalpha = -c1c*exp6/oma2;
        } else if (alpha < thr2) {
            fca = 0.0;
            dfca_dalpha = 0.0;
        } else {
            double exp7 = std::exp(c2c/oma);
            fca = -dc*exp7;
            dfca_dalpha = -dc*c2c*exp7/oma2;
        }

        double dec1_dnup = dec1_drs*drs_dn + dec1_dzeta*dzeta_dnup + dec1_dp*dp_dn;
        double dec1_dndn = dec1_drs*drs_dn + dec1_dzeta*dzeta_dndn + dec1_dp*dp_dn;
        double dec1_dagr = dec1_dp*dp_dagr;

        double dec0_dnup = dec0_drs*drs_dn + dec0_dzeta*dzeta_dnup + dec0_dp*dp_dn;
        double dec0_dndn = dec0_drs*drs_dn + dec0_dzeta*dzeta_dndn + dec0_dp*dp_dn;
        double dec0_dagr = dec0_dp*dp_dagr;

        double dfca_dnup   = dfca_dalpha*dalpha_dnup;
        double dfca_dndn   = dfca_dalpha*dalpha_dndn;
        double dfca_dagr   = dfca_dalpha*dalpha_dagr;
        double dfca_dtauup = dfca_dalpha*dalpha_dtauup;
        double dfca_dtaudn = dfca_dalpha*dalpha_dtaudn;

        double ec = ec1 + fca*(ec0 - ec1);

        double fnupc = n*(dec1_dnup + dfca_dnup*(ec0 - ec1) + fca*(dec0_dnup - dec1_dnup)) + ec;
        double fndnc = n*(dec1_dndn + dfca_dndn*(ec0 - ec1) + fca*(dec0_dndn - dec1_dndn)) + ec;

        double fdnupc = 0.0;
        double fdndnc = 0.0;
        double fdnc   = n*(dec1_dagr + dfca_dagr*(ec0 - ec1) + fca*(dec0_dagr - dec1_dagr));

        double fdtauupc = n*dfca_dtauup*(ec0 - ec1);
        double fdtaudnc = n*dfca_dtaudn*(ec0 - ec1);

        // ========= final mix =========
        xce[i]           = x_parameter*ex     + c_parameter*ec;

        fn[i]           = x_parameter*fnupx  + c_parameter*fnupc;
        fn[n2ft3d + i]  = x_parameter*fndnx  + c_parameter*fndnc;

        fdn[i]             = x_parameter*fdnupx + c_parameter*fdnupc;
        fdn[n2ft3d + i]    = x_parameter*fdndnx + c_parameter*fdndnc;
        fdn[2*n2ft3d + i]  = c_parameter*fdnc;

        fdtau[i]          = x_parameter*fdtauupx + c_parameter*fdtauupc;
        fdtau[n2ft3d + i] = x_parameter*fdtaudnx + c_parameter*fdtaudnc;
    }
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
// C++ translation of Fortran subroutine gen_SCAN_restricted
// Notes:
//  - Arrays are assumed 0-based in C++ (i = 0..n2ft3d-1).
//  - This code calls gen_PW91_c_rz(...), which you must provide.
//  - Uses double precision throughout.


// You must implement/provide this (Fortran: call gen_PW91_c_rz(tol,rs,zeta,...))

void gen_SCAN_restricted(const int n2ft3d,
                         double *rho_in, double *agr_in, double *tau_in,
                         const double x_parameter, const double c_parameter,
                         double* xce, double* fn, double* fdn, double* fdtau)
{
    // constants

    // density cutoff parameters

    // SCAN constants
    double b1, b2;
    double b4;

    constexpr double delta = 5.0e-5;

    b2 = std::sqrt(5913.0/405000.0);
    b1 = (511.0/13500.0)/(2.0*b2);
    b4 = muAK*muAK/K1 - 1606.0/18225.0 - b1*b1;

    const double Cx   = (-0.75)*std::pow(3.0/pi, thrd);
    const double P23  = std::pow(3.0*pi*pi, twthrd);
    const double P13  = std::pow(3.0/(4.0*pi), thrd);
    const double P23t = std::pow(3.0*pi*pi/16.0, twthrd);

    for (int i = 0; i < n2ft3d; ++i)
    {
        double n   = rho_in[i] + ETA;
        double agr = agr_in[i] + ETA;
        double tau = 2.0*tau_in[i] + ETA;

        double n_13 = std::pow(n, thrd);
        double n_53 = n_13*n_13*n;
        double n_83 = n_53*n;
        double agr2 = agr*agr;

        double p       = agr2/(4.0*P23*n_83);
        double p_14    = std::sqrt(std::sqrt(p));
        double dp_dn   = -etthrd*p/n;
        double dp_dagr =  2.0*p/agr;

        double tauW = 0.125*agr2/n;
        double tauU = 0.3*P23*n_53;
        double z    = tauW/tau;
        double dz_dn   = -z/n;
        double dz_dagr =  2.0*z/agr;
        double dz_dtau = -z/tau;

        double z2 = z*z;

        double alpha = 0.0, dalpha_dn = 0.0, dalpha_dagr = 0.0, dalpha_dtau = 0.0;

        double thr0 = agr2/(8.0*tau);
        double dthr0_dagr = 2.0*thr0/agr;
        double dthr0_dtau = -thr0/tau;

        double delta2 = delta*delta;
        double diff  = n - thr0;
        double diff2 = diff*diff;
        double diff3 = diff2*diff;
        double dod   = diff/delta;
        double dod2  = dod*dod;
        double tauUn = tauU*n;

        if (n <= thr0)
        {
            alpha = 0.0; dalpha_dn = 0.0; dalpha_dagr = 0.0; dalpha_dtau = 0.0;
        }
        else if (n > thr0 && n < (thr0 + delta))
        {
            double falpha       = tau*(-diff3/delta2 + 2.0*diff2/delta);
            double dfalpha_dn   = tau*(-3.0*dod2 + 4.0*dod);
            double dfalpha_dagr = tau*(3.0*dod2 - 4.0*dod)*dthr0_dagr;
            double dfalpha_dtau = falpha/tau + tau*(3.0*dod2 - 4.0*dod)*dthr0_dtau;

            alpha       = falpha/tauUn;
            dalpha_dn   = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn;
            dalpha_dagr = dfalpha_dagr/tauUn;
            dalpha_dtau = dfalpha_dtau/tauUn;
        }
        else // n >= thr0 + delta
        {
            double falpha       = tau*n - agr2/8.0;
            double dfalpha_dn   = tau;
            double dfalpha_dagr = -agr/4.0;
            double dfalpha_dtau = n;

            alpha       = falpha/tauUn;
            dalpha_dn   = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn;
            dalpha_dagr = dfalpha_dagr/tauUn;
            dalpha_dtau = dfalpha_dtau/tauUn;
        }

        double oma  = 1.0 - alpha;
        double oma2 = oma*oma;

        // ***** SCAN Exchange *****
        double exp1 = std::exp(-b4*p/muAK);
        double exp2 = std::exp(-b3*oma2);

        double x1 = muAK*p*(1.0 + (b4*p/muAK)*exp1);
        double x2 = b1*p + b2*oma*exp2;
        double x  = x1 + x2*x2;

        double denh1x = K1 + x;
        double numh1x = denh1x + K1*x;
        double h1x    = numh1x/denh1x;

        double dx1_dp    = muAK + b4*p*exp1*(2.0 - p*b4/muAK);
        double dx2_dp    = b1;
        double dx_dp     = dx1_dp + 2.0*x2*dx2_dp;
        double dx_dalpha = 2.0*b2*exp2*x2*(2.0*b3*oma2 - 1.0);

        double dh1x_dx     = std::pow(K1/denh1x, 2.0);
        double dh1x_dp     = dh1x_dx*dx_dp;
        double dh1x_dalpha = dh1x_dx*dx_dalpha;

        double gx = 1.0, dgx_dp = 0.0;
        if (p_14 < 0.002)
        {
            gx = 1.0; dgx_dp = 0.0;
        }
        else
        {
            double exp3 = std::exp(-a1/p_14);
            gx     = 1.0 - exp3;
            dgx_dp = -0.25*a1*exp3/(p*p_14);
        }

        double fxa = 0.0, dfxa_dalpha = 0.0;
        if (alpha <= thr1)
        {
            double exp4 = std::exp(-c1x*alpha/oma);
            fxa = exp4;
            dfxa_dalpha = -c1x*exp4/oma2;
        }
        else if (alpha > thr1 && alpha < thr2)
        {
            fxa = 0.0;
            dfxa_dalpha = 0.0;
        }
        else // alpha >= thr2
        {
            double exp5 = std::exp(c2x/oma);
            fxa = -dx*exp5;
            dfxa_dalpha = -dx*c2x*exp5/oma2;
        }

        double Fx = (h1x + fxa*(h0x - h1x))*gx;

        double dFx_dp     = dgx_dp*(h1x + fxa*(h0x - h1x))
                          + gx*dh1x_dp*(1.0 - fxa);
        double dFx_dalpha = gx*(dh1x_dalpha + dfxa_dalpha*(h0x - h1x) - fxa*dh1x_dalpha);

        double dFx_dn   = dFx_dalpha*dalpha_dn   + dFx_dp*dp_dn;
        double dFx_dagr = dFx_dalpha*dalpha_dagr + dFx_dp*dp_dagr;
        double dFx_dtau = dFx_dalpha*dalpha_dtau;

        double ex0  = Cx*n_13;
        double nex0 = n*ex0;

        double ex     = ex0*Fx;
        double fnx    = nex0*dFx_dn   + frthrd*ex;
        double fdnx   = nex0*dFx_dagr;
        double fdtaux = nex0*dFx_dtau;

        // ***** SCAN Correlation *****
        double zeta = 0.0;

        double rs     = P13/n_13;
        double drs_dn = -thrd*rs/n;
        double rs_12  = std::sqrt(rs);

        double beta      = beta0*(1.0 + beta1*rs)/(1.0 + beta2*rs);
        double dbeta_drs = beta0*(beta1 - beta2)/std::pow(1.0 + beta2*rs, 2.0);

        double ecLDA1, decLDA1_drs, decLDA1_dzeta;
        gen_PW91_c_rz(tol, rs, zeta, &ecLDA1, &decLDA1_drs, &decLDA1_dzeta);

        double w1fac   = ecLDA1/gamma;
        double expw1   = std::exp(-w1fac);
        double w1      = expw1 - 1.0;
        double dw1_drs = -expw1*decLDA1_drs/gamma;

        double A      = beta/(gamma*w1);
        double dA_drs = dbeta_drs/(gamma*w1) - A*dw1_drs/w1;

        double t2      = P23t*p/rs;
        double dt2_drs = -t2/rs;
        double dt2_dp  =  P23t/rs;

        double At2     = A*t2;
        double dengAt2 = 1.0 + 4.0*At2;
        double gAt2    = 1.0/(std::sqrt(std::sqrt(dengAt2)));

        double tmp1      = 1.0 + w1*(1.0 - gAt2);
        double dtmp1_drs = dw1_drs*(1.0 - gAt2)
                         + gAt2*w1*(t2*dA_drs + A*dt2_drs)/dengAt2;
        double dtmp1_dp  = gAt2*w1*A*dt2_dp/dengAt2;

        double H1      = gamma*std::log(tmp1);
        double dH1_drs = gamma*dtmp1_drs/tmp1;
        double dH1_dp  = gamma*dtmp1_dp/tmp1;

        double ec1      = ecLDA1 + H1;
        double dec1_drs = decLDA1_drs + dH1_drs;
        double dec1_dp  = dH1_dp;

        double ecLDA0      = -b1c/(1.0 + b2c*rs_12 + b3c*rs);
        double decLDA0_drs =  b1c*(b3c + 0.5*b2c/rs_12)
                            /std::pow(1.0 + b2c*rs_12 + b3c*rs, 2.0);

        double w0fac   = ecLDA0/b1c;
        double expw0   = std::exp(-w0fac);
        double w0      = expw0 - 1.0;
        double dw0_drs = -decLDA0_drs*expw0/b1c;

        double ginf     = 1.0/(std::sqrt(std::sqrt(1.0 + 4.0*xi*p)));
        double dginf_dp = -xi*ginf/(1.0 + 4.0*xi*p);

        double tmp0      = 1.0 + w0*(1.0 - ginf);
        double dtmp0_drs = dw0_drs*(1.0 - ginf);
        double dtmp0_dp  = -w0*dginf_dp;

        double H0      = b1c*std::log(tmp0);
        double dH0_drs = b1c*dtmp0_drs/tmp0;
        double dH0_dp  = b1c*dtmp0_dp/tmp0;

        double ec0      = ecLDA0 + H0;
        double dec0_drs = decLDA0_drs + dH0_drs;
        double dec0_dp  = dH0_dp;

        double fca = 0.0, dfca_dalpha = 0.0;
        if (alpha <= thr1)
        {
            double exp6 = std::exp(-c1c*alpha/oma);
            fca = exp6;
            dfca_dalpha = -c1c*exp6/oma2;
        }
        else if (alpha >= thr1 && alpha <= thr2)
        {
            fca = 0.0;
            dfca_dalpha = 0.0;
        }
        else // alpha > thr2
        {
            double exp7 = std::exp(c2c/oma);
            fca = -dc*exp7;
            dfca_dalpha = -dc*c2c*exp7/oma2;
        }

        double dec1_dn   = dec1_drs*drs_dn + dec1_dp*dp_dn;
        double dec1_dagr = dec1_dp*dp_dagr;
        double dec0_dn   = dec0_drs*drs_dn + dec0_dp*dp_dn;
        double dec0_dagr = dec0_dp*dp_dagr;

        double dfca_dn   = dfca_dalpha*dalpha_dn;
        double dfca_dagr = dfca_dalpha*dalpha_dagr;
        double dfca_dtau = dfca_dalpha*dalpha_dtau;

        double ec     = ec1 + fca*(ec0 - ec1);
        double fnc    = n*(dec1_dn + dfca_dn*(ec0 - ec1) + fca*(dec0_dn - dec1_dn)) + ec;
        double fdnc   = n*(dec1_dagr + dfca_dagr*(ec0 - ec1) + fca*(dec0_dagr - dec1_dagr));
        double fdtauc = n*dfca_dtau*(ec0 - ec1);

        xce[i]   = x_parameter*ex     + c_parameter*ec;
        fn[i]    = x_parameter*fnx    + c_parameter*fnc;
        fdn[i]   = x_parameter*fdnx   + c_parameter*fdnc;
        fdtau[i] = x_parameter*fdtaux + c_parameter*fdtauc;
    }
}


} // namespace pwdft
