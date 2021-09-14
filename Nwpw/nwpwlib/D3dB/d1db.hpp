#ifndef _D1dB_H_
#define _D1dB_H_
/* d1db.h
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#include	"Parallel.hpp"
#include	"Mapping1.hpp"

namespace pwdft {
using namespace pwdft;

class d1db : public Mapping1 {


public:
        Parallel  *parall;

        /* constructor */
	d1db(Parallel *, const int, const int, int *);

        /* destructor */
        ~d1db() {}

};

}

#endif
