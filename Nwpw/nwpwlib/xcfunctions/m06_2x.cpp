/* m06_2x.cpp
  Author - Eric Bylaska
*/

#include <cmath>

namespace pwdft {


/****************************************
 *					*
 *	    gen_M06_2x_unrestricted	*
 *					*
 ****************************************/
/*
     This function returns the M06_2x exchange-correlation
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
void gen_M06_2x_unrestricted(const int n2ft3d, double *dn_in, double *agr_in, double *tau_in,
                             const double x_parameter, const double c_parameter, 
                             double *xce, double *fn, double *fdn, double *fdtau) 
{
}

/****************************************
 *					*
 *	    gen_M06_2x_restricted       *
 *					*
 ****************************************/
/*
   This routine calculates the M06_2x exchange-correlation
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

void gen_M06_2x_restricted(const int n2ft3d, double *rho_in, double *agr_in, double *tau_in,
                           const double x_parameter, const double c_parameter,
                           double *xce, double *fn, double *fdn, double *fdtau) 
{
}


} // namespace pwdft
