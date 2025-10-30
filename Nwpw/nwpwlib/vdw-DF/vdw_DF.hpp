#ifndef _VDW_DF_HPP_
#define _VDW_DF_HPP_
/* vdw_DF.hpp
   Author - Eric Bylaska

*/


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

int Nqs;

   //initialize r and k grids
   //r-grid 
   int nr  = 2048;
   int nr1 = nr+1;
   double rmax = 100.0;
   double dr = rmax / static_cast<double>(nr);

   // kgrid - maximum g=64 and gg=4096 ... 
   int  nk   = 1024;
   int  nk1  = nk+1;
   double kmax = 64.0;
   double dk = kmax / static_cast<double>(nr);

   double *g, *phir, *phir0, *rphi, *sphi, *utmp, *xtmp;
   double *phik, *phik0, *phik2, *a, *a2

   double *aweights, *cos_a, *sina, *nu, *nu1, *Wab,
   double *qmesh


 /* constructor */
   vdw_DF(Parallel *, const int, const int, const int, const int);
 
   /* destructor */

    /**
     * @brief Destructor for the vdw_DF class.
     */
   ~vdw_DF();







};

} // namespace pwdft

#endif
                 
