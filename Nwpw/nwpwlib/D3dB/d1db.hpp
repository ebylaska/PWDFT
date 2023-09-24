#ifndef _D1dB_HPP_
#define _D1dB_HPP_
/* d1db.hpp
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#include "Mapping1.hpp"
#include "Parallel.hpp"
#include "gdevice2.hpp"

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

class d1db : public Mapping1 {

int request1_indx,request2_indx;

public:

  Parallel *parall;

  /* constructor */
  d1db(Parallel *, const int, const int, int *);

  /* destructor */
  ~d1db() {  
      if (parall->np_j()>1) {
         parall->aend(request1_indx);
         parall->aend(request2_indx);
      }
   }

  /* SUMMA matrix mutiplications */
  void DMatrix_dgemm1(Parallel *parall, gdevice2 *,
                      int, int, int, int,
                      double, double *, int, int *, int *,
                              double *, int, int *, int *,
                      double, double *, int, int *, int *,
                      double *, double *);

  void DMatrix_dgemm2(Parallel *parall, gdevice2 *,
                      int, int, int, int,
                      double, double *, int, int *, int *,
                              double *, int, int *, int *,
                      double, double *, int, int *, int *,
                      double *, double *);

  void DMatrix_dgemm3(Parallel *parall, gdevice2 *,
                      int, int, int, int,
                      double, double *, int, int *, int *,
                              double *, int, int *, int *,
                      double, double *, int, int *, int *,
                      double *, double *);

  void DMatrix_dgemm2c(Parallel *, gdevice2 *,
                          int, int, int, int,
                          double *, double *, int, int *, int *, int *,
                          double *, int, int *, int *,
                          double  *, double *);

  void DMatrix_dgemm1_rot2(Parallel *, gdevice2 *,
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
