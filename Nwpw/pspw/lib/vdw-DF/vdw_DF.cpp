
#include <cmath>

#include "vdw.hpp"
#include "Control2.hpp"
#include "compressed_io.hpp"
#include "util.hpp"

#include "NwpwLibraryVdwConfig.hpp"

#include "v_exc.hpp"
#include "v_dirac.hpp"
#include "vdw_DF.hpp"

#include <iostream>
#include "iofmt.hpp"
#include <cstring>


/*
------------------------------------------------------------------------------
  NOTE ON k-POINT INDEPENDENCE OF THE VDW-DF FUNCTIONAL

  The nonlocal van der Waals correlation (vdW-DF) depends only on the 
  *total* electron density rho(r) and related local response quantities
  (e.g., grad rho and q0).  No orbital- or k-resolved terms enter the 
  functional.  Therefore, vdW-DF is formally independent of the Bloch 
  vector k, and the vdW energy does not depend on Brillouin-zone sampling.

  Practically, the vdW kernel is evaluated using densities in real space
  and standard FFT machinery.  Although the intermediate arrays are complex
  (because FFTs are complex-valued in general), the physical quantities
  entering vdW-DF are real and contain no explicit k dependence.

  Consequence: vdW-DF can be evaluated identically for Gamma-only and
  k-point runs, using the same kernel tables and the same real-space
  density (summed over k if needed).  No special handling of k-points
  is required inside the vdW routines.
------------------------------------------------------------------------------
*/


namespace pwdft {

//private functions here

#include <cstddef>

/**********************************************
 *                                            *
 *           vdw_DF_init_poly (C++)           *
 *                                            *
 **********************************************/
static void vdw_DF_init_poly(
    int n,
    const double *x,   // length n  (qmesh)
    double *ya,        // length n*n, column-major: ya(i,j) = ya[j*n + i]
    double *y2a,       // length n*n, same layout
    double *utmp       // length n (workspace)
)
{
    const double yp1 = 0.0;
    const double ypn = 0.0;

    // ya(i,j) = δ_ij, Fortran: do i=1,n; do j=1,n
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            ya[j * n + i] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Fortran: do j=1,n; call nwpw_spline(x, ya(1,j), n, yp1, ypn, y2a(1,j), utmp)
    // Column j is contiguous at &ya(j,1) → &ya[j*n]
    for (int j = 0; j < n; ++j)
    {
        util_spline(
            x,                 // xa
            &ya[j * n],        // ya(1,j)
            n,
            yp1,
            ypn,
            &y2a[j * n],       // y2a(1,j)
            utmp
        );
    }
}


/**************************************
 *                                    *
 *            vdw_DF::init_poly       *
 *                                    *
 **************************************/
/**
 * @brief Initialize cubic–spline interpolation tables for the vdW kernel.
 *
 * Constructs the spline basis ya,y2a over the q–mesh so that subsequent
 * polynomial evaluations poly(i,x) are O(1).  Must be called after reading
 * qmesh from vdw_kernels.dat.  The result is independent of k–points.
 */
void vdw_DF::init_poly()
{
    std::vector<double> utmp(Nqs);
    vdw_DF_init_poly(Nqs, qmesh, ya, ya2, utmp.data());
}


/**************************************
 *                                    *
 *            vdw_DF::poly            *
 *                                    *
 **************************************/
/**
 * @brief Evaluate the i-th polynomial p(x) and its derivative dp(x) for q0.
 *
 * Implements the cubic-spline evaluation of the vdW kernel “shape functions”
 * (Román-Pérez & Soler table).  Used in both θ and u-function construction.
 *
 * @param i   1-based kernel index (1…Nqs)
 * @param x   local q0 value
 * @param p   output spline value
 * @param dp  output derivative dp/dx
 */
inline void vdw_DF::poly(int i, double x, double &p, double &dp)
{
    // i is 1-based in Fortran, so convert to 0-based in C++:
    int idx = i - 1;

    if ((x >= qmin) && (x <= qmax))
    {
        // Fortran: nx = dsqrt(x/dbl_mb(qmesh(1)+Nqs-1)) * Nqs
        // qmesh(1) = qmesh[0]
        // qmesh(Nqs) = qmesh[Nqs-1]

        double denom = qmesh[Nqs - 1];   // same as (qmesh(1)+Nqs-1)
        double nx_f  = std::sqrt(x / denom) * Nqs;
        int nx = static_cast<int>(nx_f);

        if (nx < 1)       nx = 1;
        if (nx > Nqs - 1) nx = Nqs - 1;

        // Offsets for this i:
        // ya(1 + Nqs*(i-1))  becomes  ya[idx*Nqs + 0]
        // y2a(1 + Nqs*(i-1)) becomes  y2a[idx*Nqs + 0]
        const double *ya_i  = ya  + idx*Nqs;
        const double *y2a_i = ya2 + idx*Nqs;

        // Evaluate cubic spline and derivative
        p  = util_splint(qmesh, ya_i,  y2a_i, Nqs, nx, x);
        dp = util_dsplint(qmesh, ya_i, y2a_i, Nqs, nx, x);
    }
    else
    {
        p  = 0.0;
        dp = 0.0;
    }
}

/**************************************
 *                                    *
 *       vdw_DF::generate_ufunc       *
 *                                    *
 **************************************/
/**
 * @brief Contract θ(G,j) with tabulated kernel φ to form u(G,i).
 *
 * Performs the Román-Pérez/Soler convolution in reciprocal space using the
 * spline-interpolated kernel φ and the complex θ(G,j).  Complex arithmetic
 * arises from FFTs only; the vdW functional itself has no explicit k-dependence.
 *
 * @param ufunc complex u(G,i) array, accumulated over (i,j) pairs.
 */
void vdw_DF::generate_ufunc(const int nk1,
                           const int Nqs,
                           const double *gphi,              // [nk1]
                           const double *phi,               // [nk1 * 2 * Npairs], Npairs = Nqs*(Nqs+1)/2
                           const int npack0,
                           const int nfft3d,
                           const double *Gpack,             // [npack0]
                           const int *nxpack,               // [npack0]
                           const std::complex<double> *theta, // [nfft3d * Nqs]
                           std::complex<double> *ufunc)       // [nfft3d * Nqs]
{
    // ---- parallel info (OpenMP-style + 2D decomp over j) ----
    const int taskid_j  = myparall->taskid_j();   // like Parallel2d_taskid_j()
    const int np_j      = myparall->np_j();       // like Parallel2d_np_j()

    // ---- zero ufunc (all G, all i) ----
    std::fill(ufunc, ufunc + static_cast<std::size_t>(nfft3d) * Nqs,
              std::complex<double>(0.0, 0.0));

    // Helper lambdas for indexing (Fortran: ufunc(k,i), theta(k,j))
    auto U = [=](int k, int i) -> std::complex<double>& {
        // k, i are 0-based; layout is (k, i) with leading dimension nfft3d
        return ufunc[i * nfft3d + k];
    };
    auto TH = [=](int k, int j) -> const std::complex<double>& {
        return theta[j * nfft3d + k];
    };

    // Number of (i,j) pairs: Nqs*(Nqs+1)/2, visited in the same order as Fortran
    int pcount = 0;
    int indx   = 0;   // 0-based; Fortran indx = 1..Npairs

    for (int j=0; j<Nqs; ++j)          // Fortran j=1..Nqs
    {
        for (int i=0; i<=j; ++i)       // Fortran i=1..j
        {
            ++indx;                        // Fortran indx = indx + 1
            const int owner = pcount;      // pcount modulo np_j

            if (owner == taskid_j)
            {
                // Pointer to φ(:, :, indx) in C layout:
                // Fortran: phi(nk1,2,Npairs)
                // => C: phi[ (indx-1)*2*nk1 + p*nk1 + k ], p=0,1 ; k=0..nk1-1
                const double *phi_block = phi + static_cast<std::size_t>(indx - 1) * 2 * nk1;

                // Loop over packed G indices, split amongst threads
                for (int k=0; k<npack0; ++k)   // Fortran k=tid+1,npack0,nthr
                {
                    const double g = Gpack[k];

                    // nxpack(k) gives an index into gphi, 1-based  fortran indx
                    const int klo = nxpack[k];      // 1-based, i.e. 1..nk
                    const int khi = klo + 1;        // 1-based, i.e., 2..nk+1

                    const int klo0 = klo - 1;       // 0-based, i.e.  0..nk-1
                    const int khi0 = khi - 1;       // 0-based, i.e.  1..nk

                    const double h = gphi[khi0] - gphi[klo0];
                    const double a = (gphi[khi0] - g) / h;
                    const double b = (g - gphi[klo0]) / h;

                    const double *phi1 = phi_block + 0 * nk1;  // φ(:,1,indx)
                    const double *phi2 = phi_block + 1 * nk1;  // φ(:,2,indx)

                    const double f =
                        a * phi1[klo0]
                      + b * phi1[khi0]
                      + ((a*a*a - a) * phi2[klo0]
                      +  (b*b*b - b) * phi2[khi0]) * h * h / 6.0;

                    //if ((k==0)) std::cout << "taskid_j=" << taskid_j << " k=" << k << " i=" << i << " j=" << j  << " f=" << f << " U=" << Efmt(13,6) << U(k,i) << " klo0=" << klo0 << " khi0=" << khi0  << std::endl;

                    // Update ufunc(k,i) and ufunc(k,j), k is 0-based
                    U(k, i) += TH(k, j) * f;
                    if (i != j)
                        U(k, j) += TH(k, i) * f;
                }
            }

            // advance pair counter and wrap over np_j (same as mod(pcount+1,np_j))
            pcount = (pcount + 1) % np_j;
        }
    }

    // OpenMP barrier in Fortran: !$OMP BARRIER
    // If you compile with OpenMP, you can add:
    // #pragma omp barrier

    // Global sum over all 2D-j tasks: D1dB_Vector_SumAll(2*Nqs*nfft3d, ufunc)
    // Treat complex array as contiguous doubles.
    const int len = 2 * Nqs * nfft3d;  // 2 for real+imag
    myparall->Vector_SumAll(2,len, reinterpret_cast<double*>(ufunc));
}



/**************************************
 *                                    *
 *     vdw_DF::generate_rho           *
 *                                    *
 **************************************/
/**
 * @brief Build total electron density ρ(r) from spin densities dn.
 *
 * Adds spin channels plus a numerical cutoff.  Required prior to θ(G,j).
 *
 * @param dn  input density dn(r,spin)
 * @param rho output total density (real space)
 */
void vdw_DF::generate_rho(const int ispin,
                         const int n2ft3d,
                         const double* dn,   // dn[n2ft3d][ispin]
                         double* rho)
{
    const double dncut = 1.0e-30;

    // Loop i 
    for (int i=0; i<n2ft3d; ++i)
    {
        // dn is laid out as dn(k, spin)
        double dn1 = dn[i +         0*n2ft3d];    // dn(i,1)
        double dn2 = dn[i + (ispin-1)*n2ft3d];    // dn(i,2)

        rho[i] = dn1 + dn2 + dncut;
    }

}

// C++ version of vdw_DF_Generate_potentials
/**************************************
 *                                    *
 *     vdw_DF::generate_potentials    *
 *                                    *
 **************************************/
/**
 * @brief Contract u(r,i) with spline-weighted q0 response to form
 *        nonlocal correlation energy density and potentials.
 *
 * Takes real-space u(r,i) (after inverse FFT), multiplies by p_j(q0) and
 * derivatives, and accumulates nonlocal energy exc(r), potential fn(r,spin),
 * and density-derivative fdn(r,spin).  Fully consistent with the Fortran driver.
 *
 * @note vdW-DF potentials depend only on ρ and gradients; no explicit k terms.
 */
void vdw_DF::generate_potentials(int Nqs,
                                 int nfft3d,
                                 int ispin,
                                 int n2ft3d,
                                 double *ufunc,     // [n2ft3d * Nqs] real-space ufunc(r,j)
                                 double *q0,        // xce
                                 double *drho,      // xcp
                                 double *ddrho,     // xxe
                                 double *tmpexc,    // xce + n2ft3d
                                 double *tmpfn,     // xxp
                                 double *tmpfdn,    // rho
                                 double *exc,       // output exc(n2ft3d)
                                 double *fn,        // output fn(n2ft3d, ispin)
                                 double *fdn)
{
    // ---- 2D decomposition in j (same logic as generate_theta_g) ----
    int taskid_j = myparall->taskid_j();
    int np_j     = myparall->np_j();

    int base = Nqs / np_j;
    int rem  = Nqs % np_j;

    int nj, jstart;
    if (taskid_j < rem)
    {
        nj     = base + 1;
        jstart = taskid_j * (base + 1);                // 0-based
    }
    else
    {
        nj     = base;
        jstart = rem * (base + 1) + (taskid_j - rem) * base;
    }

    // ---- zero local accumulators ----
    std::fill(tmpexc, tmpexc + n2ft3d, 0.0);
    for (int ms = 0; ms < ispin; ++ms)
    {
        std::fill(tmpfn  + ms * n2ft3d, tmpfn  + (ms + 1) * n2ft3d, 0.0);
        std::fill(tmpfdn + ms * n2ft3d, tmpfdn + (ms + 1) * n2ft3d, 0.0);
    }

    // NOTE: In the original Fortran, there is a call:
    //   Grsm_gh_fftb(nfft3d, nj, ufunc(1,jstart))
    // which does a G->R back FFT on ufunc.
    // Here we assume ufunc is already in real space (or that the FFT
    // is handled elsewhere). If you later wire a proper G->R FFT,
    // that call should go here.

    if (nj > 0)
    {
        //mygrid->nbngh_fftb(0,nj,ufunc+jstart*n2ft3d);
        for (int j=0; j<nj; ++j) 
        {
           int jj = jstart + j;
           mygrid->c_unpack(0,ufunc+jj*n2ft3d);
           mygrid->cr_pfft3b(0,ufunc+jj*n2ft3d);
        }
       
        // ---- main loops: i = real-space grid, j = q index block ----
        for (int i=0; i<n2ft3d; ++i)
        {
            double q0i = q0[i];

            for (int jj = 0; jj < nj; ++jj)
            {
                int j = jstart + jj;   // 0-based index for this processor’s j

                double pj, dpj;
                // poly() uses Fortran-style 1-based i, so pass (j+1)
                // and returns pj, dpj
                // (poly is the C++ version of vdw_DF_poly)
                // signature: void vdw_DF::poly(int i, double x, double &p, double &dp);
                // We'll call it through a lambda capturing i, but since this
                // helper is free-standing, we assume a global poly-like function
                // exists or is provided by the caller.
                // ==> Here we just *declare* it and assume the definition
                //     below in this file:
                //extern void vdw_DF_poly_wrapper(int i, double x, double &p, double &dp);
                poly(j+1, q0i, pj, dpj);

                double u = ufunc[i + j * n2ft3d];  // ufunc(r,i,j)

                tmpexc[i] += 0.5 * pj * u;

                for (int ms = 0; ms < ispin; ++ms)
                {
                    double dr  = drho [i + ms * n2ft3d];
                    double ddr = ddrho[i + ms * n2ft3d];

                    tmpfn [i + ms * n2ft3d] += (pj + dpj * dr) * u;
                    tmpfdn[i + ms * n2ft3d] += (dpj * ddr) * u;
                }
            }
        }
    }

    // ---- MPI reductions (sum across all ranks) ----
    myparall->Vector_SumAll(2, n2ft3d,        tmpexc);
    myparall->Vector_SumAll(2, ispin*n2ft3d,  tmpfn);
    myparall->Vector_SumAll(2, ispin*n2ft3d,  tmpfdn);

    // ---- accumulate into global exc, fn, fdn ----
    for (int i=0; i<n2ft3d; ++i)
        exc[i] += tmpexc[i];

    for (int ms=0; ms<ispin; ++ms)
    {
        double *fn_ms  = fn  + ms * n2ft3d;
        double *fdn_ms = fdn + ms * n2ft3d;
        double *tmpfn_ms  = tmpfn  + ms * n2ft3d;
        double *tmpfdn_ms = tmpfdn + ms * n2ft3d;

        for (int i=0; i<n2ft3d; ++i)
        {
            fn_ms [i] += tmpfn_ms [i];
            fdn_ms[i] += tmpfdn_ms[i];
        }
    }
}



/**************************************
 *                                    *
 *     vdw_DF::generate_theta_g       *
 *                                    *
 **************************************/
/**
 * @brief Construct the nonlocal density combination θ(G,j).
 *
 * Builds q0(r) and θ(r,j) = ρ(r)*p_j(q0(r)) on real space, and FFTs to G.
 * All exchange-correlation screening from vdW-DF is handled here.  The result
 * depends only on total density and gradients, not on Bloch k.
 *
 * @param Nqs    number of kernel points
 * @param ispin  number of spin channels (1 or 2)
 * @param theta  complex θ(G,j) array [nfft3d*Nqs]
 */
void vdw_DF::generate_theta_g(
        int Nqs,
        int nfft3d,
        int ispin,
        int n2ft3d,
        double Zab,
        double qmin,
        double qmax,
        const double *rho,
        const double *agr,
        double *vxc,
        double *exc,
        double *vxx,
        double *exx,
        double *theta)
{
    // ---- get parallel info ----
    //int tid     = myparall->threadid();
    //int nthr    = myparall->nthreads();

    int taskid_j = myparall->taskid_j();
    int np_j     = myparall->np_j();

    int nx = mygrid->nx;
    int ny = mygrid->ny;
    int nz = mygrid->nz;
    double scal1 = 1.0 / double(nx * ny * nz);

    // ---- block decomposition over j ----
    int base = Nqs / np_j;
    int rem  = Nqs % np_j;

    int nj, jstart;
    if (taskid_j < rem) 
    {
       nj = base + 1;
       jstart = taskid_j * (base + 1);
    } 
    else 
    {
       nj = base;
       jstart = rem * (base + 1) + (taskid_j - rem) * base;
    }

    // zero theta
    std::fill(theta, theta + size_t(n2ft3d)*Nqs, 0.0);

    if (nj > 0)
    {
        const double pi = M_PI;
        const double Cf  = pow(3.0 * pi * pi, 1.0/3.0);
        const double A   = -4.0 * pi / 3.0;
        const double Cxi = Zab / (36.0 * Cf);

        const double onethird = 1.0/3.0;
        const double frthrd   = 4.0/3.0;
        const double elthrd   = 11.0/3.0;
        const double dncut    = 1e-12;

        // ---- main loop over real-space grid ----
        for (int i=0; i<n2ft3d; ++i)
        {
           double rh = rho[i];

           double q0sat = 0.0;

           if (rh >= dncut)
           {
              double temp = (agr[i] / pow(rh, frthrd));
              double xi  = Cxi * temp * temp;
              double dxi = 2.0 * Cxi * agr[i] * pow(1.0 / pow(rh, frthrd), 2);

              double exgga = exx[i] - xi * exx[i];
              double vxgga = vxx[i] - xi * (vxx[i] - elthrd * exx[i]);

              if ((vxgga - exgga) < 0.0) 
              {
                 xi = 0.0;
                 dxi = 0.0;
              }

              double q0 = A * (exc[i] - exx[i] * xi);

              double dq0drho[2] = {0.0, 0.0};
              double dq0ddrho   = 0.0;

              if (q0 < qmin) 
              {
                 q0 = qmin;
              } 
              else 
              {
                 dq0drho[0] = A*((vxc[i] - exc[i]) - xi*(vxx[i] - elthrd*exx[i]))/rh;
                 if (ispin==2)
                     dq0drho[1] = A*((vxc[i+n2ft3d] - exc[i])
                                   - xi*(vxx[i+n2ft3d] - elthrd*exx[i]))/rh;

                 dq0ddrho = -A * exx[i] * dxi;
              }

              double tsum  = 0.0;
              double dtsum = 0.0;
              for (int k = 1; k <= 12; k++) 
              {
                 tsum  += std::pow(q0/qmax, k) / double(k);
                 dtsum += std::pow(q0/qmax, k-1);
              }
              q0sat = qmax * (1.0 - std::exp(-tsum));
              double dq0satdq0 = std::exp(-tsum) * dtsum;

              exc[i] = q0sat;
              for (int ms = 0; ms < ispin; ++ms)
                 vxc[i + ms*n2ft3d] = rh * dq0satdq0 * dq0drho[ms];

              exx[i] = rh * dq0satdq0 * dq0ddrho;
           }
           else
           {
              q0sat = qmax;
              exc[i] = qmax;
              for (int ms=0; ms<ispin; ++ms)
                  vxc[i + ms*n2ft3d] = 0.0;
              exx[i] = 0.0;
           }

           // ---- now fill θ(i,j) for this processor's j-block ----
           for (int j=0; j<nj; j++)
           {
              int jj = jstart + j;
              double pj, dpj;
              poly(jj+1, q0sat, pj, dpj);
              theta[i + jj*n2ft3d] = rho[i] * pj;
           }
        }


        // ---- transform θ(r,j) → θ(G,j) ----
        //for (int j=0; j<nj; j++) 
        //{
        //   int jj = jstart + j;
        //   mygrid->r_zero_ends(theta + jj*n2ft3d);
           //std::cout << "jj=" << jj << " rho*theta_r=" << Efmt(13,8) << mygrid->rr_dot(rho,theta+jj*n2ft3d) << std::endl;
        //   mygrid->rc_pfft3f(0,theta+jj*n2ft3d);
        //   mygrid->c_pack(0,theta+jj*n2ft3d);
        //   mygrid->c_pack_SMul(0,scal1,theta+jj*n2ft3d);
        //}
        mygrid->nh_zero_ends(nj,theta +jstart*n2ft3d);
        mygrid->nbnhg_fftf(0,nj,theta +jstart*n2ft3d);
        mygrid->nbnhg_SMul(0,nj,scal1,theta +jstart*n2ft3d);
    }

    // ---- global reduce sum across processors ----
    myparall->Vector_SumAll(2, Nqs*n2ft3d, theta);
}



/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/**
 * @brief Constructor for the `d3db` class.
 *
 * This constructor initializes an instance of the `d3db` class with the given parameters.
 *
 * @param inparall A pointer to a `Parallel` object.
 * @param inmaptype An integer specifying the mapping type.
 * @param nx The number of grid points in the x-direction.
 * @param ny The number of grid points in the y-direction.
 * @param nz The number of grid points in the z-direction.
 */
vdw_DF::vdw_DF(Pneb *inmygrid, Control2 &control, bool is_vdw2)
{
   mygrid   = inmygrid;
   myparall = mygrid->d3db::parall;
   has_vdw = true;

   nfft3d = mygrid->nfft3d;
   npack0 = mygrid->npack(0);
   n2ft3d = mygrid->n2ft3d;

   bool oprint = (myparall->is_master() && control.print_level("medium"));


   const std::string nwpw_vdw_qmesh = std::string(Nwpw_LIBRARYVDW_Default) + "/vdw_qmesh.dat";
   //const char *nwpw_libraryps = Nwpw_LIBRARYPS_Default + "/VDW/vdw_qmesh.dat";

   char datafile[256];
   strcpy(datafile, "vdw_kernels.dat");
   control.add_permanent_dir(datafile);
   //std::cout << "nwpw_vdw = " << nwpw_vdw_qmesh << " datafile=" << datafile << std::endl;

   int ifound = cfileexists(datafile);
   if (ifound == 0)
   {
      if (oprint) std::cout << "Generating VDW kernel filename:" << datafile << std::endl;
      //vdw_DF_kernel_gen_data(datafile)
      vdw_DF_kernel_gen_data(myparall,datafile,nwpw_vdw_qmesh.c_str());
   }

   if (myparall->is_master())
   {
      openfile(5,datafile,"r");
      iread(5,&Nqs,1);
      iread(5,&nk,1);
      dread(5,&kmax,1);
   }
   myparall->Brdcst_iValue(0,MASTER,&Nqs);
   myparall->Brdcst_iValue(0,MASTER,&nk);
   myparall->Brdcst_Values(0,MASTER,1,&kmax);
   nk1 = nk + 1;


   theta  = new (std::nothrow) double[Nqs*n2ft3d]();
   ufunc  = new (std::nothrow) double[Nqs*n2ft3d]();

   qmesh  = new (std::nothrow) double[Nqs]();
   ya     = new (std::nothrow) double[Nqs*Nqs]();
   ya2    = new (std::nothrow) double[Nqs*Nqs]();
   gphi   = new (std::nothrow) double[nk1]();
   phi    = new (std::nothrow) double[nk1*Nqs*(Nqs+1)]();

   xcp    = new (std::nothrow) double[2*n2ft3d]();
   xce    = new (std::nothrow) double[2*n2ft3d]();
   xxp    = new (std::nothrow) double[2*n2ft3d]();
   xxe    = new (std::nothrow) double[2*n2ft3d]();
   rho    = new (std::nothrow) double[2*n2ft3d]();
   Gpack  = new (std::nothrow) double[npack0]();
   nxpack = new (std::nothrow) int[npack0]();


   if (myparall->is_master())
   {
      dread(5,qmesh,Nqs);
      dread(5,phi,nk1*Nqs*(Nqs+1));
   }
   myparall->Brdcst_Values(0, MASTER, Nqs, qmesh);
   myparall->Brdcst_Values(0, MASTER, nk1*Nqs*(Nqs+1), phi);

    
   init_poly();   // <<<=== Required (Fortran does it right after reading qmesh)

   double dk = kmax/double(nk);
   for (int k=0; k<=nk; ++k)
      gphi[k] = k*dk;

   double *Gx = mygrid->Gpackxyz(0,0);
   double *Gy = mygrid->Gpackxyz(0,1);
   double *Gz = mygrid->Gpackxyz(0,2);
   for (int k=0; k<npack0; ++k)
   {
      double gg = std::sqrt(Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k]);
      Gpack[k] = gg;
         
      int nx = gg/dk;
      nxpack[k] = util_splint_nx(gphi,nx,gg,nk1);
      //std::cout << k << " " << nx << " " << gg << " " << nxpack[k] << std::endl;
   }



   if (is_vdw2)
      Zab = -1.887;
   else
      Zab = -0.8491;

   qmax = qmesh[Nqs-1];
   qmin = 0.0;
}



/********************************
 *                              *
 *       vdw_DF::evaluate       *
 *                              *
 ********************************/
/**
 * @brief Top-level vdW-DF driver: compute nonlocal correlation and potentials.
 *
 * Performs the full pipeline:
 *   (1) LDA pieces, (2) total density, (3) θ(G,j), (4) u(G,i),
 *   (5) exc, fn, and fdn.
 *
 * The entire vdW-DF evaluation is independent of Brillouin-zone sampling:
 * only total real-space density (summed over k, if present) enters.
 */
void vdw_DF::evaluate(int ispin, const double *dn, const double *agr,
                      double *exc, double *fn, double *fdn)
{
    // 1. LDA pieces
    v_exc(ispin, n2ft3d, const_cast<double*>(dn), xcp, xce, rho);
    //std::cout << "VXCA = " << Ffmt(13,9) << mygrid->rr_dot(dn,xce) << std::endl;
    //std::cout << "VXCAA= " << Ffmt(13,9) << mygrid->rr_dot(dn,xcp) << std::endl;
    if (ispin == 1)
        v_dirac(ispin, n2ft3d, const_cast<double*>(dn), xxp, xxe, rho);
    //std::cout << "VXCB = " << Ffmt(13,9) << mygrid->rr_dot(dn,xxp) << std::endl;

    // 2. build rho
    generate_rho(ispin, n2ft3d, dn, rho);
    mygrid->r_zero_ends(rho);
    //std::cout << "RHOA = " << Ffmt(13,9) << mygrid->rr_dot(rho,rho) << std::endl;

    // 3. theta(G)
    generate_theta_g(Nqs, nfft3d, ispin, n2ft3d, Zab, qmin, qmax, rho, agr, xcp, xce, xxp, xxe, theta);

  // 2. IMPORTANT: generate_theta_g modifies xcp,xce,xxp,xxe
  //    Now compute the theta-contractions EXACTLY as Fortran does:
  
  //double dum1 = mygrid->rr_dot(rho, xcp);
  //double dum2 = mygrid->rr_dot(rho, xce);
  //std::cout << "theta xcp = " << Ffmt(13,9) << dum1 << std::endl;
  //std::cout << "theta xce = " << Ffmt(13,9) << dum2 << std::endl;
  
  //dum1 = mygrid->rr_dot(rho, xxp);
  //dum2 = mygrid->rr_dot(rho, xxe);
  //std::cout << "theta xxp = " << Ffmt(13,9) << dum1 << std::endl;
  //std::cout << "theta xxe = " << Ffmt(13,9) << dum2 << std::endl;
    //std::cout << "theta xcp = " << Ffmt(13,9) << mygrid->rr_dot(rho,xcp) << std::endl;
    //std::cout << "theta xce = " << Ffmt(13,9) << mygrid->rr_dot(rho,xce) << std::endl;
    //std::cout << "theta xxp = " << Ffmt(13,9) << mygrid->rr_dot(rho,xxp) << std::endl;
    //std::cout << "theta xxe = " << Ffmt(13,9) << mygrid->rr_dot(rho,xxe) << std::endl;

  //dum1 = mygrid->cc_pack_dot(0,theta, theta);
  //std::cout << "theta*theta = " << Efmt(13,9) << dum1 << std::endl;
  //for (auto jj=0; jj<Nqs; ++jj)
 // {
  //   dum1 = mygrid->cc_pack_dot(0,theta+jj*n2ft3d, theta+jj*n2ft3d);
   //  std::cout << "jj=" << jj << " theta*theta = " << Efmt(13,9) << dum1 << std::endl;
 // }



    // 4. ufunc(G,i)
    generate_ufunc(nk1, Nqs, gphi, phi, npack0, nfft3d, Gpack, nxpack,
                   reinterpret_cast<const std::complex<double>*>(theta), 
                   reinterpret_cast<std::complex<double>*>(ufunc));
    //std::cout << "rho*ufunc = " << Ffmt(13,9) << mygrid->rr_dot(rho,ufunc) << std::endl;

    //std::cout << "pre xce = " << Ffmt(13,9) << mygrid->rr_dot(rho,xce) << std::endl;

    // 5. exc, fn, fdn
    generate_potentials(Nqs, nfft3d, ispin, n2ft3d, ufunc, xce, xcp, xxe,
                        xce+n2ft3d, xxp, rho, exc, fn, fdn);
    //std::cout << "xce = " << Ffmt(13,9) << mygrid->rr_dot(rho,xce) << std::endl;
}



}



