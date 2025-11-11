#ifndef _VDW_DF_HPP_
#define _VDW_DF_HPP_
/* vdw_DF.hpp
   Author - Eric Bylaska

*/

#include "PGrid.hpp"
#include "Parallel.hpp"
#include "Control2.hpp"

namespace pwdft {

/**
 * @class d1db
 * @brief Container for operations on distributed 3D blocks.
 *
 * The `d1db` class is designed to handle operations related to distributed 3D blocks.
 * It provides functionality for performing matrix multiplications and other mathematical
 * operations on distributed data structures. This class is often used in conjunction
 * with the `Parallel` and `gdevice2` classes to manage parallelism and device operations.
 */

class vdw_DF  {

   PGrid    *mygrid;
   Parallel *myparall;


   //initialize r and k grids
   //r-grid 
   //int nr  = 2048;
   //int nr1 = nr+1;
   //double rmax = 100.0;
   //double dr = rmax / static_cast<double>(nr);

   // kgrid - maximum g=64 and gg=4096 ... 
   //int  nk   = 1024;
   //int  nk1  = nk+1;
   //double kmax = 64.0;
   //double dk = kmax / static_cast<double>(nr);

   //double *g, *phir, *phir0, *rphi, *sphi, *utmp, *xtmp;
   //double *phik, *phik0, *phik2, *a, *a2;

   //double *aweights, *cos_a, *sina, *nu, *nu1, *Wab;

   int Nqs, nk, nk1;
   int nfft3d,n2ft3d,npack0;
   double kmax;


   double *qmesh, *ya, *ya2, *gphi, *phi, *theta, *ufunc;
   double *xcp, *xce, *xxp, *xxe, *rho, *Gpack, *nxpack;


   /* constructor */
   vdw_DF(PGrid *, Control2 &);
 
   /* destructor */

    /**
     * @brief Destructor for the vdw_DF class.
     */
   ~vdw_DF() {
       delete [] qmesh;
       delete [] ya;
       delete [] ya2;
       delete [] gphi;
       delete [] phi;
       delete [] theta;
       delete [] ufunc;
       delete [] xcp;
       delete [] xce;
       delete [] xxp;
       delete [] xxe;
       delete [] rho;
       delete [] Gpack;
       delete [] nxpack;
    }



};

} // namespace pwdft

#endif
                 
