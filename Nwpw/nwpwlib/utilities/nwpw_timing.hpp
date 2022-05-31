#ifndef _NWPW_TIMING_HPP_
#define _NWPW_TIMING_HPP_

#include        "nwpw_timers.hpp"

namespace pwdft {


extern void nwpw_timing_start(const int);
extern void nwpw_timing_end(const int);
extern void nwpw_timing_print_final(int);


class nwpw_timing_function {
   int id;

public:

   /* constructor */
   nwpw_timing_function(const int i) {
      id = i;
      nwpw_timing_start(id);
   }

   /* destructor */
   ~nwpw_timing_function() {nwpw_timing_end(id);}
};


}
#endif
