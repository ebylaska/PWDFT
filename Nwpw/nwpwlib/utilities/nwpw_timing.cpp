

#include "nwpw_timers.hpp"

namespace pwdft {

static nwpw_timers mytimer;

void nwpw_timing_start(const int i) { mytimer.start_timer(i); }
void nwpw_timing_end(const int i) { mytimer.end_timer(i); }
void nwpw_timing_print_final(int count, std::ostream &coutput) {
  mytimer.print_final(count, coutput);
}

} // namespace pwdft
