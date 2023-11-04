#ifndef _C1dB_HPP_
#define _C1dB_HPP_
/* c1db.hpp
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#include "Mapping1.hpp"
#include "Parallel.hpp"
#include "gdevice2.hpp"

namespace pwdft {

/**
 * @class c1db
 * @brief Container for operations on distributed 3D blocks.
 *
 * The `c1db` class is designed to handle operations related to distributed 3D blocks.
 * It provides functionality for performing matrix multiplications and other mathematical
 * operations on distributed data structures. This class is often used in conjunction
 * with the `Parallel` and `gdevice2` classes to manage parallelism and device operations.
 */

class c1db : public Mapping1 {

int request1_indx,request2_indx;

public:

  Parallel *parall;

  /* constructor */
  c1db(Parallel *, const int, const int, int *);

  /* destructor */
  ~c1db() {  
      if (parall->np_j()>1) {
         parall->aend(request1_indx);
         parall->aend(request2_indx);
      }
   }

  /* SUMMA matrix mutiplications */
  void CMatrix_dgemm1(Parallel *parall, gdevice2 *,
                      int, int, int, int,
                      double, double *, int, int *, int *,
                              double *, int, int *, int *,
                      double, double *, int, int *, int *,
                      double *, double *);

  void CMatrix_dgemm2(Parallel *parall, gdevice2 *,
                      int, int, int, int,
                      double, double *, int, int *, int *,
                              double *, int, int *, int *,
                      double, double *, int, int *, int *,
                      double *, double *);

  void CMatrix_dgemm3(Parallel *parall, gdevice2 *,
                      int, int, int, int,
                      double, double *, int, int *, int *,
                              double *, int, int *, int *,
                      double, double *, int, int *, int *,
                      double *, double *);

  void CMatrix_dgemm2c(Parallel *, gdevice2 *,
                          int, int, int, int,
                          double *, double *, int, int *, int *, int *,
                          double *, int, int *, int *,
                          double  *, double *);

  void CMatrix_dgemm1_rot2(Parallel *, gdevice2 *,
                           int, int, int, int,
                           double,
                           double *, int, int *, int *,
                           double *, int, int *, int *,
                           double,
                           double *, int, int *, int *,
                           double *, double *, double *, double *);

};

} // namespace pwdft

#endif
