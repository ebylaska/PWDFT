#ifndef _cVDW_DF_HPP_
#define _cVDW_DF_HPP_
/* cvdw_DF.hpp
   Author - Eric Bylaska

*/

//#include "PGrid.hpp"
#include "Pneb.hpp"
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

class cvdw_DF  {

private: 

   Pneb     *mygrid;
   Parallel *myparall;


   bool has_vdw = false;
   int Nqs, nk, nk1;
   int nfft3d,n2ft3d,npack0;
   double kmax;

   double qmin, qmax, Zab;   // <<< ADD THESE


   double *theta_c, *ufunc_c;
   double *theta, *ufunc;
   double *qmesh, *ya, *ya2, *gphi, *phi;
   double *xcp, *xce, *xxp, *xxe, *rho, *Gpack;
   int *nxpack;

   // private functions
   void init_poly();
   void poly(int, double, double &, double &);
   void generate_rho(const int, const int, const double *, double *);
   void generate_theta_g(int, int, int, int, double, double, double, const double *, const double *,
                         double *, double *, double *, double *, double *);
   void generate_ufunc(const int, const int, const double *, const double *, const int, const int,
                       const double *, const int *, const std::complex<double> *, std::complex<double> *);
   void generate_potentials(int, int, int, int, double *, double *, double *, double *, 
                            double *, double *, double *, double *, double *, double *);



public:

   /* constructor */
   cvdw_DF(Pneb *, Control2 &, bool);
 
   /* destructor */

    /**
     * @brief Destructor for the vdw_DF class.
     */
   ~cvdw_DF() {
       delete [] theta_c;
       delete [] ufunc_c;
       delete [] theta;
       delete [] ufunc;
       delete [] qmesh;
       delete [] ya;
       delete [] ya2;
       delete [] gphi;
       delete [] phi;
       delete [] xcp;
       delete [] xce;
       delete [] xxp;
       delete [] xxe;
       delete [] rho;
       delete [] Gpack;
       delete [] nxpack;
    }



   bool exist(){return has_vdw;}
   //void v_vdw(int, int, double *, double *, double * double *, double * double *);
   void evaluate(int, const double *, const double *, double *,  double *, double *);


};

} // namespace pwdft

#endif
                 
