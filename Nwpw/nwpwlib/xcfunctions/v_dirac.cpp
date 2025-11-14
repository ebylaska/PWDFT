// v_dirac.cpp
#include <cmath>

namespace pwdft {

/**
 * C++ port of Fortran:
 *
 *   subroutine v_dirac(ispin,n2ft3d,dn,xcp,xce,x)
 *   real*8 dn(n2ft3d,2), xcp(n2ft3d,2), xce(n2ft3d,2), x(n2ft3d)
 *
 * Arrays are stored in Fortran column-major order.
 * We keep that layout:
 *   dn(k,ms)  -> dn[k + ms*n2ft3d], ms = 0,1
 *   xcp(k,ms) -> xcp[k + ms*n2ft3d]
 *   xce(k,ms) -> xce[k + ms*n2ft3d]
 */
void v_dirac(const int ispin,
             const int n2ft3d,
             double *dn,   // n2ft3d*2
             double *xcp,  // n2ft3d*2
             double *xce,  // n2ft3d*2
             double *x)    // n2ft3d (scratch)
{
    // ---- constants (from Fortran PARAMETERs) ----
    const double one3rd  = 1.0 / 3.0;
    const double for3rd  = 4.0 / 3.0;
    const double one6th  = 1.0 / 6.0;
    const double dncut   = 1.0e-30;

    // Vosko et al parameters (kept for completeness, not used directly below)
    const double ap  = 3.109070e-02;
    const double af  = 1.554530e-02;
    const double x0p = -1.049800e-01;
    const double x0f = -3.250000e-01;
    const double bp  = 3.727440e+00;
    const double bf  = 7.060420e+00;
    const double cp  = 1.293520e+01;
    const double cf  = 1.805780e+01;

    // Precomputed constants from Vosko parameters
    const double xp   = -4.581653e-01;
    const double xf   = -5.772521e-01;
    const double qp   =  6.151991e+00;
    const double qf   =  4.730927e+00;
    const double xx0p =  1.255491e+01;
    const double xx0f =  1.586879e+01;

    const double cp1  =  3.109070e-02;
    const double cf1  =  1.554530e-02;
    const double cp2  =  9.690228e-04;
    const double cf2  =  2.247860e-03;
    const double cp3  =  1.049800e-01;
    const double cf3  =  3.250000e-01;
    const double cp4  =  3.878329e-02;
    const double cf4  =  5.249122e-02;
    const double cp5  =  3.075995e+00;
    const double cf5  =  2.365463e+00;
    const double cp6  =  1.863720e+00;
    const double cf6  =  3.530210e+00;

    const double dp1  =  6.218140e-02;
    const double df1  =  3.109060e-02;
    const double dp2  =  1.938045e-03;
    const double df2  =  4.495720e-03;
    const double dp3  =  1.049800e-01;
    const double df3  =  3.250000e-01;
    const double dp4  = -3.205972e-02;
    const double df4  = -1.779316e-02;
    const double dp5  = -1.192972e-01;
    const double df5  = -1.241661e-01;
    const double dp6  =  1.863720e+00;
    const double df6  =  3.530210e+00;
    const double dp7  =  9.461748e+00;
    const double df7  =  5.595417e+00;

    const double fc   =  1.923661e+00;
    const double fd   =  2.564881e+00;
    const double crs  =  7.876233e-01;


    const double pi = 4.0 * std::atan(1.0);

    // Convenience lambdas for column-major indexing
    auto idx = [n2ft3d](int k, int ms) -> int {
        return k + ms * n2ft3d; // k: 0..n2ft3d-1, ms: 0 or 1
    };

    // ---- 1. sqrt of Wigner radius: x(k) = crs / rho^(1/6) ----
    // rho = dn(k,1) + dn(k,ispin) + dncut
    // Fortran is 1-based; here we are 0-based.
    // ispin can be 1 or 2, but arrays have 2 columns.
    const int ms_spin = (ispin == 1) ? 0 : 1;

    for (int k = 0; k < n2ft3d; ++k)
    {
        double rho = dn[idx(k,0)] + dn[idx(k,ms_spin)] + dncut;
        // rho**(1/6)
        double rho_16 = std::pow(rho, one6th);
        x[k] = crs / rho_16;
    }

    // ---- 2. Paramagnetic exchange energy & potential ----
    // xce(k,1) = xp / x(k)^2
    // xcp(k,1) = (4/3) * xce(k,1)
    for (int k = 0; k < n2ft3d; ++k)
    {
        double val = xp / (x[k] * x[k]);
        xce[idx(k,0)] = val;
        xcp[idx(k,0)] = for3rd * val;
    }

    // ---- 3. Return if spin-restricted Dirac (ispin == 1) ----
    if (ispin == 1)
    {
        return;
    }

    // ---- 4. Ferromagnetic exchange energy & potential ----
    // xce(k,2) = xf / x(k)^2
    // xcp(k,2) = (4/3) * xce(k,2)
    for (int k = 0; k < n2ft3d; ++k)
    {
        double val = xf / (x[k] * x[k]);
        xce[idx(k,1)] = val;
        xcp[idx(k,1)] = for3rd * val;
    }

    // ---- 5. Spin-polarized exchange potential ----
    for (int k = 0; k < n2ft3d; ++k)
    {
        double rho = dn[idx(k,0)] + dn[idx(k,1)] + dncut;
        double zup = 2.0 * dn[idx(k,0)] / rho;
        double zdw = 2.0 * dn[idx(k,1)] / rho;

        // f = (zup*zup^(1/3) + zdw*zdw^(1/3) - 2)*fc
        double zup13 = std::pow(zup, one3rd);
        double zdw13 = std::pow(zdw, one3rd);
        double f = (zup * zup13 + zdw * zdw13 - 2.0) * fc;

        double xcp_p = xcp[idx(k,0)];
        double xcp_f = xcp[idx(k,1)];
        double xce_p = xce[idx(k,0)];
        double xce_f = xce[idx(k,1)];

        // mix paramagnetic/ferromagnetic exchange potentials
        xcp[idx(k,0)] = (1.0 - f) * xcp_p + f * xcp_f;

        // df = (zup^(1/3) - zdw^(1/3)) * (xce_f - xce_p) * fd
        double df = (zup13 - zdw13) * (xce_f - xce_p) * fd;

        // new xcp(k,2) and xcp(k,1)
        xcp[idx(k,1)] = xcp[idx(k,0)] - zup * df;
        xcp[idx(k,0)] = xcp[idx(k,0)] + zdw * df;

        // xce(k,1) = xce_p + f*(xce_f - xce_p)
        xce[idx(k,0)] = xce_p + f * (xce_f - xce_p);
    }

}

} // namespace pwdft

