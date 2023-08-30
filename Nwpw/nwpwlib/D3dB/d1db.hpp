#ifndef _D1dB_HPP_
#define _D1dB_HPP_
/* d1db.h
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#include "Mapping1.hpp"
#include "Parallel.hpp"
#include "gdevice2.hpp"

namespace pwdft {

class d1db : public Mapping1 {

public:
  Parallel *parall;

  /* constructor */
  d1db(Parallel *, const int, const int, int *);

  /* destructor */
  ~d1db() {}

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
};

} // namespace pwdft

#endif
