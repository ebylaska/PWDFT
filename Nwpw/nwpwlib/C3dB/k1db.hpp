#ifndef _K1dB_HPP_
#define _K1dB_HPP_
/* k1db.hpp
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#include "Mapping1k.hpp"
#include "Parallel.hpp"
#include "gdevice2.hpp"

namespace pwdft {

/**
 * @class k1db
 * @brief Container for operations on distributed k point blocks.
 *
 * The `k1db` class is designed to handle operations related to distributed k point blocks.
 * It provides functionality for performing matrix multiplications and other mathematical
 * operations on distributed data structures. This class is often used in conjunction
 * with the `Parallel` and `gdevice2` classes to manage parallelism and device operations.
 */

class k1db : public Mapping1k {

int request1_indx,request2_indx;

public:

  Parallel *parall;

  /* constructor */
  k1db(Parallel *, const int, const int);

  /* destructor */
  ~k1db() {  
      if (parall->np_k()>1) {
         //parall->aend(request1_indx);
         //parall->aend(request2_indx);
      }
   }


};

} // namespace pwdft

#endif
