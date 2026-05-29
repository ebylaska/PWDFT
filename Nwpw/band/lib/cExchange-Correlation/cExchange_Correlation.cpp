/*
 *
 *   gga = -1 off
 *
 *   **** LDAs ****
 *   gga = 0 vosko,
 *   gga = 1-9 (reserved for other lda's)
 *
 *   **** GGAs ****
 *   gga = 10 pbe96, pbe
 *   gga = 11 blyp
 *   gga = 12 revpbe
 *   gga = 13 pbesol
 *   gga = 14 hser remainder   (-0.25*Ex(w,pbe,sr) + Ex(pbe) + Ec(pbe))
 *   gga = 15 b3lypr remainder
 *   gga = 16 BEEF
 *   gga = 17 XBEEF-CPBE
 *   gga = 17-99 (reserved for other gga's)
 *
 *   **** hybrids ****
 *   gga = 100-109 (reserved for lda hybrids)
 *   gga = 110  pbe0
 *   gga = 111  ?????
 *   gga = 112  revpbe0
 *   gga = 113  bnl
 *   gga = 114  hse
 *   gga = 115  b3lyp
 *   gga = 116-199 (reserved for hybrids)
 *   gga = 200 hartree-fock
 *
 *   **** meta ggas and hybrid metagga ****
 *   gga = 300 vs98
 *   gga = 301 tpss03
 *   gga = 302 scan
 *   gga = 303 pkzb
 *   gga = 304 m06-l
 *   gga = 305  m06
 *   gga = 306  m06-2x
 *
 *   **** control dispersion block ****
 *   has_disp = false
 *   is_grimme2 = false
 *   has_vdw  = false
 *   is_vdw2 = false
 *   options_disp = ''
 *  
 *   std::string xc_name,options_disp;
 *   int gga;
 *   bool use_lda, use_gga, use_mgga;
 *
 *   bool has_disp = false;
 *   bool has_vdw  = false;
 *   bool is_grimme2,is_vdw2;
 */

#include "cExchange_Correlation.hpp"
#include "v_exc.hpp"
#include "v_cwexc.hpp"
#include "v_cmexc.hpp"
#include "cvdw_DF.hpp"
#include <algorithm>
#include "parsestring.hpp"

namespace pwdft {

/* Constructors */

/*******************************************
 *                                         *
 *        cXC_Operator::cXC_Operator       *
 *                                         *
 *******************************************/

cXC_Operator::cXC_Operator(Cneb *mygrid, Control2 &control) 
{
    mycneb = mygrid;
    xc_name = control.xc_name();
 
    const std::string name = mystring_lowercase(xc_name);
 
    auto add_disp = [&](const std::string& func_flag) {
       if (!has_disp) return;
       // ensure exactly one -func … prefix
       if (!func_flag.empty()) {
          // always keep a leading space before subsequent flags
          options_disp = "-func " + func_flag + (options_disp.empty() ? "" : " " + options_disp);
       }
    };
 
 
    gga = 0;
 
    // set the grimme and vdw options
    if      (mystring_contains(name, "-grimme2")) {has_disp = true; is_grimme2 = true;  options_disp = "-old -noprint";}
    else if (mystring_contains(name, "-grimme3")) {has_disp = true; is_grimme2 = false; options_disp = "-zero -noprint";}
    else if (mystring_contains(name, "-grimme4")) {has_disp = true; is_grimme2 = false; options_disp = "-bj -num -noprint";}
    else if (mystring_contains(name, "-grimme5")) {has_disp = true; is_grimme2 = false; options_disp = "-zerom -noprint";}
    else if (mystring_contains(name, "-grimme6")) {has_disp = true; is_grimme2 = false; options_disp = "-bjm -num -noprint";}
 
    if      (mystring_contains(name, "-vdw2"))    {has_vdw = true; is_vdw2 = true; }
    else if (mystring_contains(name, "-vdw"))     {has_vdw = true; is_vdw2 = false; }
 
    // set the gga options
 
    // ---- XC family: choose longest/specific first; first match wins ----
    // Hybrids first (specific → general)
    if      (mystring_contains(name, "revpbe0")) {gga = 112; add_disp("revpbe0");}
    else if (mystring_contains(name, "pbe0"))    {gga = 110; add_disp("pbe0");}
    else if (mystring_contains(name, "hse"))     {gga = 114; add_disp("hse06");}
    else if (mystring_contains(name, "bnl"))     {gga = 113; add_disp("hse06");}  // if that’s really what you want
    else if (mystring_contains(name, "b3lypr"))  {gga = 115; add_disp("b3-lyp");} // treat b3lypr as hybrid remainder path
    else if (mystring_contains(name, "blyp0"))   {gga = 111; add_disp("b3-lyp");} // if this alias is desired
 
    // Non-hybrid GGAs (specific → general)
    else if (mystring_contains(name, "revpbe"))     {gga = 12; add_disp("revpbe");}
    else if (mystring_contains(name, "pbesol"))     {gga = 13; add_disp("pbesol");}
    else if (mystring_contains(name, "pbe96"))      {gga = 10; add_disp("pbe");}
    else if (mystring_contains(name, "xbeef-cpbe")) {gga = 17; add_disp("pbesol");} // confirm
    else if (mystring_contains(name, "beef"))       {gga = 16; add_disp("pbesol");} // confirm
    else if (mystring_contains(name, "blyp"))       {gga = 11; add_disp("b-lyp");}
    else if (mystring_contains(name, "pbe"))        {gga = 10; add_disp("pbe");}
 
    // LDA family
    else if (mystring_contains(name, "vosko") || mystring_contains(name, "lda")) {gga = 0;}
 
    // HF exact exchange - dispersion func flag generally not needed/used here
    if (mystring_contains(name, "hartree-fock") || name == "hf" || name.find(" hf ") != std::string::npos) {gga = 200;}
 
    // Meta-GGAs - future options
    if      (mystring_contains(name, "scan"))   gga = 302;
    else if (mystring_contains(name, "tpss03")) gga = 301;
    else if (mystring_contains(name, "vs98"))   gga = 300;
    else if (mystring_contains(name, "pkzb"))   gga = 303;
    else if (mystring_contains(name, "m06-2x")) gga = 306;
    else if (mystring_contains(name, "m06-l"))  gga = 304;
    else if (mystring_contains(name, "m06"))    gga = 305;
 
    // HF exact exchange - dispersion func flag generally not needed/used here
    if (mystring_contains(name, "hartree-fock") || name == "hf" || name.find(" hf ") != std::string::npos) {gga = 200;}
 
    // Meta-GGAs - future options
    if      (mystring_contains(name, "scan"))   gga = 302;
    else if (mystring_contains(name, "tpss03")) gga = 301;
    else if (mystring_contains(name, "vs98"))   gga = 300;
    else if (mystring_contains(name, "pkzb"))   gga = 303;
    else if (mystring_contains(name, "m06-2x")) gga = 306;
    else if (mystring_contains(name, "m06-l"))  gga = 304;
    else if (mystring_contains(name, "m06"))    gga = 305;
 
    use_lda = false;
    use_gga = false;
    use_mgga = false;
    if (gga == 0) 
    {
       use_lda = true;
       xtmp = new double[mycneb->ispin * mycneb->nfft3d];
    }
    if ((gga >= 10) && (gga < 100)) 
    {
       use_gga = true;
       if (mycneb->ispin == 1) 
       {
          rho = new double[mycneb->n2ft3d]; // real
         
          grx = new double[mycneb->n2ft3d]; // complex
          gry = new double[mycneb->n2ft3d]; // complex
          grz = new double[mycneb->n2ft3d]; // complex
         
          agr = new double[mycneb->n2ft3d]; // real|complex
          fn  = new double[mycneb->n2ft3d]; // real|complex
          fdn = new double[mycneb->n2ft3d]; // real|complex
       } 
       else 
       {
          rho = new double[2 * mycneb->n2ft3d]; // real
         
          grx = new double[3 * mycneb->n2ft3d]; // complex
          gry = new double[3 * mycneb->n2ft3d]; // complex
          grz = new double[3 * mycneb->n2ft3d]; // complex
         
          agr = new double[3 * mycneb->n2ft3d]; // real|complex
          fn  = new double[2 * mycneb->n2ft3d]; // real|complex
          fdn = new double[3 * mycneb->n2ft3d]; // real|complex
       }
    }
    if ((gga >= 300))
    {
        use_mgga = true;
        if (mycneb->ispin == 1) 
        {
            rho = new double[mycneb->n2ft3d];
            
            grx = new double[mycneb->n2ft3d];
            gry = new double[mycneb->n2ft3d];
            grz = new double[mycneb->n2ft3d];
            
            agr = new double[mycneb->n2ft3d];
            fn  = new double[mycneb->n2ft3d];
            fdn = new double[mycneb->n2ft3d];
          
            tau    = new double[mycneb->n2ft3d];
            dfdtau = new double[mycneb->n2ft3d];
        }  
        else
        {
            rho = new double[2 * mycneb->n2ft3d];
            
            grx = new double[3 * mycneb->n2ft3d];
            gry = new double[3 * mycneb->n2ft3d];
            grz = new double[3 * mycneb->n2ft3d];
            
            agr = new double[3 * mycneb->n2ft3d];
            fn  = new double[2 * mycneb->n2ft3d];
            fdn = new double[3 * mycneb->n2ft3d];
          
            tau    = new double[2 * mycneb->n2ft3d];
            dfdtau = new double[2 * mycneb->n2ft3d];
        }
    }

 
    if (has_vdw)
    {
       myvdw = new cvdw_DF(mygrid,control,is_vdw2);
    }
 
  
    // std::cout << "xc_name =" << xc_name << std::endl;
}

/*******************************************
 *                                         *
 *        cXC_Operator::v_exc_all           *
 *                                         *
 *******************************************/
void cXC_Operator::v_exc_all(int ispin, double *dn, double *xcp, double *xce) 
{
   if (use_lda) 
   {
      v_exc(ispin, mycneb->nfft3d, dn, xcp, xce, xtmp);

      //std::cout << "dn=" << dn[0] << " " << dn[1] << std::endl;
      //std::cout << "xcp=" << xcp[0] << " " << xcp[1] << std::endl;
      //double sumall = 0.0;
      //for (auto i=0; i<mycneb->nfft3d; ++i)
      // {
      //    std::cout << "i=" << i << " dnall=" << dn[i] << " xcp=" << xcp[i] << std::endl;
      //    sumall += dn[i];
      // }
      // std::cout << "sumall=" << sumall << std::endl;
   } 
   else if (use_gga) 
   {
       v_cwexc(gga, mycneb, myvdw, dn, 1.0, 1.0, xcp, xce, rho, grx, gry, grz, agr, fn, fdn);
   } 
   else if (use_mgga) 
   {
       std::cout << "v_cmexc into" << std::endl;
       v_cmexc(gga, mycneb, myvdw, dn, tau, 1.0, 1.0, xcp, xce, rho, grx, gry, grz, agr, fn, fdn, dfdtau);
       std::cout << "v_cmexc OUT " << std::endl;
   }
}



/*******************************************
 *                                         *
 *        XC_Operator::gga_gen_tau         *
 *                                         *
 *******************************************/
void cXC_Operator::gga_gen_tau(const int ispin, const int neq[2], const double *psi)
{
    std::cout << "GERA" << std::endl;

    // 1. Check if Meta-GGA is active (mapped from use_mgga)
    if (!this->use_mgga) {
        return;
    }

    // 2. Setup orbital ranges
    // Fortran: n1(1)=1, n2(1)=neq(1)... (Adjusted for 0-based C++)
    int n1[2], n2[2];
    n1[0] = 0;
    n2[0] = neq[0];
    n1[1] = neq[0];
    n2[1] = neq[0] + neq[1];

    // 3. Scaling factor
    // lattice_omega() must be a member of XC_Operator or Pneb
    double scal2 = 0.5 / this->mycneb->lattice->omega();

    std::cout << "GERB" << std::endl;

    // 4. Prepare temporary buffer for dpsi
    // We use std::vector for the temporary FFT buffer to ensure RAII (auto-cleanup)
    // The size is 2 * nfft3d to account for the complex nature (real/imag)
    //size_t nfft3d = this->mycneb->nfft3d;
    size_t npack2 = 2*this->mycneb->CGrid::npack1_max();
    std::cout << "GERC" << std::endl;
    size_t n2ft3d = this->mycneb->n2ft3d;
    int    nbrillq = this->mycneb->nbrillq;
    std::cout << "GERD" << std::endl;
    std::vector<double> dpsi_tmp(n2ft3d, 0.0);
    std::cout << "GERE" << std::endl;
    double *dpsi = dpsi_tmp.data();

    std::cout << "GERF, nbrillq=" << nbrillq <<  std::endl;
    int ishift = (neq[0]+neq[1])*npack2;

    std::cout << "GERF, ispin=" << ispin << " n2ft3d="<< n2ft3d <<   std::endl;
    this->mycneb->r_nzero(ispin,tau);

    // 5. Main Computation Loop
    for (int nbq=0; nbq<nbrillq; ++nbq)
    {
    std::cout << "GERG, nbq=" << nbq << std::endl;
        double weight =  this->mycneb->pbrill_weight(nbq);
    std::cout << "GERG, nbq=" << nbq << " weight=" << weight <<   std::endl;
        int nbq1 = nbq+1;
        for (int ms=0; ms<ispin; ++ms)
        {
            // Offset for the tau array for the current spin
            size_t tau_offset = ms*n2ft3d;

            for (int n=n1[ms]; n<n2[ms]; ++n)
            {
                for (int xyz=0; xyz<3; ++xyz)
                {
                    //STEP A: Compute Gradient in Reciprocal Space
                    double *gxyz = this->mycneb->Gpackxyz(nbq1,xyz);
                    this->mycneb->tcr_pack_iMul_unpack_fft(nbq1, gxyz, psi + n*npack2 + nbq*ishift, dpsi);
                
                
                    // STEP B: sqr and add
                    this->mycneb->rr_addsqr(dpsi, this->tau + tau_offset);
                }
            }

            this->mycneb->r_SMul(weight*scal2,tau+tau_offset);
            //this->mycneb->r_zero_ends(tau+tau_offset);
        }
    }

    // 7. Final Step: Sum tau across all spins into a single array
    this->mycneb->c3db::parall->Vector_SumAll(2, nbrillq*ispin*n2ft3d,tau);
}



/*******************************************
 *                                         *
 *        cXC_Operator::meta_gga_Hpsik      *
 *                                         *
 *******************************************/
void cXC_Operator::meta_gga_Hpsik(const int ispin, const int neq[2], const double *psi, double *hpsi)
{
    // 1. Check if Meta-GGA is active (mapped from use_mgga) 
    if (!this->use_mgga) {
        return;
    }

    // 2. Setup orbital ranges             
    // Fortran: n1(1)=1, n2(1)=neq(1)... (Adjusted for 0-based C++)
    int n1[2], n2[2];                      
    n1[0] = 0;
    n2[0] = neq[0];
    n1[1] = neq[0];
    n2[1] = neq[0] + neq[1];

    double scal1 = 1.0 / ((double)((mycneb->nx) * (mycneb->ny) * (mycneb->nz)));
    double scal2 = 1.0 / this->mycneb->lattice->omega();
    double scal = 0.5*scal1;

    size_t npack2 = 2*this->mycneb->CGrid::npack1_max();
    size_t n2ft3d = this->mycneb->n2ft3d;
    size_t nbrillq = this->mycneb->nbrillq;
    std::vector<double> dpsi_tmp(n2ft3d, 0.0);
    double *dpsi = dpsi_tmp.data();

    int ishift = (neq[0]+neq[1])*npack2;

    for (int nbq=0; nbq<nbrillq; ++nbq)
    {
        double weight =  this->mycneb->pbrill_weight(nbq);
        int nbq1 = nbq+1;

        for (int ms=0; ms<ispin; ++ms) 
        {
            size_t tau_offset = ms*n2ft3d;
            for (int n=n1[ms]; n<n2[ms]; ++n) 
            {
                for (int xyz=0; xyz<3; ++xyz) 
                {
                    double *gxyz = this->mycneb->Gpackxyz(nbq1, xyz);
                    this->mycneb->tcr_pack_iMul_unpack_fft(nbq1, gxyz, psi + n*npack2 + nbq*ishift, dpsi);
 
                    this->mycneb->rr_Mul(dfdtau+tau_offset, dpsi);
 
                    this->mycneb->rc_pfft3f(nbq1, dpsi);
                    this->mycneb->c_pack(nbq1,dpsi);
                    this->mycneb->c_pack_SMul(nbq1, scal, dpsi);
 
                    this->mycneb->tc_pack_iMul(nbq1, gxyz, dpsi);
                    this->mycneb->cc_pack_Sum2(nbq1, dpsi, hpsi+n*npack2 + nbq+ishift);
                }
 
            }
        }
    }

}



/*******************************************
 *                                         *
 *        XC_Operator::meta_gga_pxc        *
 *                                         *
 *******************************************/
/**
 * @brief Computes the energy contribution from the meta-GGA functional.
 * 
 * This function calculates the integral of the product of the functional 
 * derivative (df/dtau) and the kinetic energy density (tau) over the 
 * 3D volume, weighted by the volume of a single voxel (dV).
 *
 * @param ispin The number of spin components (1 for restricted, 2 for polarized).
 * @arg ne      Array containing the number of electrons for each spin.
 * @arg psi     Pointer to the wavefunction data.
 * @return The calculated energy contribution.
 */
double  cXC_Operator::meta_gga_pxc(const int ispin, const int ne[2], const double *psi) 
{
    size_t n2ft3d = this->mycneb->n2ft3d;

    //Calculate the volume element (dV)
    double total_voxels = static_cast<double>(this->mycneb->nx * this->mycneb->ny * this->mycneb->nz);
    double dV = this->mycneb->lattice->omega() / total_voxels;


    //Sum the energy contribution across all active spins
    double pmeta = 0.0;
    for (int ms=0; ms<ispin; ++ms) 
    {
        size_t tau_offset = ms*n2ft3d;
        pmeta += dV*this->mycneb->rr_dot(dfdtau + tau_offset, tau + tau_offset);
    }
    if (ispin==1) pmeta = pmeta + pmeta;

    return pmeta;
}





} // namespace pwdft
