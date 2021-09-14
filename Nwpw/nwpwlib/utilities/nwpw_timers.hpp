#ifndef _nwpw_timers_HPP_
#define _nwpw_timers_HPP_
/* nwpw_timers.hpp
   Author - Eric Bylaska
        this class is used to keep track of timings.


*   1 - total FFT 
*   2 - total dot products
*   3 - lagrange multipliers
*   4 - exchange correlation
*   5 - local pseudopotentials
*   6 - non-local pseudopotentials
*   7 - hartree potentials
*   8 - structure factors
*   9 - masking and packing

*   10 - geodesic time
*   11 - gen psi_r and dn
*   12 - allocating memory from stack
*   13 - miscellaneous steepest descent update
*   15 - ffm_dgemm
*   16 - fmf_dgemm
*   17 - m_diagonalize
*   18 - mmm_dgemm

*
*   20 - phase factors
*   21 - ewald /ion-ion

*   22 - tredq
*   23 - getdiags
*   24 - tqliq
*   25 - eigsrt

*   30 - queue fft
*   31 - queue fft serial
*   32 - queue fft parallel
*   33 - HFX

*   34 - paw gaussian integrals
*   35 - paw atomic coulomb
*   36 - paw atomic xc
*   37 - paw gen dEmult/dQlm
*   38 - paw gen dElocal/dQlm
*   39 - paw cmp operations

*   40 - qmmm LJ
*   41 - qmmm residual Q

*   42 - MATHIAS InnerLoop
*   43 - MATHIAS Phaze
*   44 - MATHIAS Pipelined FFTs
*   45 - MATHIAS Lagrange
*   46 - MATHIAS Exch Corr
*   47 - MATHIAS Hpsi


*   50 - io time

*   52 - HFX localization
*   53 - HFX DM columns
*   54 - HFX DM Cholesky
*   55 - re-gridding
*
*   70 - timer0
*   71 - timer1
*   72 - timer2
*   73 - timer3
*

*/



#include <string>
#include <cstdio>
#include <iomanip> 
#include <iostream>
#include <chrono>

namespace pwdft {
using namespace pwdft;

#define MASTER  0
#define	nwpw_tim_max	80


class nwpw_timers {


   std::chrono::high_resolution_clock::time_point start[nwpw_tim_max+1];
   double times[nwpw_tim_max+1];

public: 

   /* constructor */
   nwpw_timers() {
      for (auto i=0; i<nwpw_tim_max; ++i)
         times[i] = 0.0;
      start[0] = std::chrono::high_resolution_clock::now();
   }

   void start_timer(const int i) {
       start[i] = std::chrono::high_resolution_clock::now();
   }
   void end_timer(const int i) {
       auto stop = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> deltatime = stop-start[i];
       times[i] += (double) deltatime.count();
   }


   void print_timer(std::string msg, double time, int counter, double ttime) {
      if (time>1.0e-9) printf(" %s %10.6le %10.6le %10.2lf%%\n",msg.c_str(),time,time/((double) counter),100*time/ttime);
      //std::cout << msg << time << " " << time/((double) counter) << " " << std::setprecision(2) << 100*time/ttime  << "%" << std::endl;
   }

   void print_final(int counter) {
       auto stop = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> deltatime = stop-start[0];
       times[0] = (double) deltatime.count();

       std::cout << " Time spent doing      total        step             percent" << std::endl;
       print_timer("total time           ", times[0],counter,times[0]);
       print_timer("total FFT time       ", times[1],counter,times[0]);
       print_timer("total dot products   ", times[2],counter,times[0]);
       print_timer("lagrange multipliers ", times[3],counter,times[0]);
       print_timer("exchange correlation ", times[4],counter,times[0]);
       print_timer("local potentials     ", times[5],counter,times[0]);
       print_timer("non-local potentials ", times[6],counter,times[0]);
       print_timer("hartree potentials   ", times[7],counter,times[0]);
       print_timer("structure factors    ", times[8],counter,times[0]);
       print_timer("masking and packing  ", times[9],counter,times[0]);
       print_timer("ffm_dgemm            ", times[15],counter,times[0]);
       print_timer("fmf_dgemm            ", times[16],counter,times[0]);
       print_timer("m_diagonalize        ", times[17],counter,times[0]);
       print_timer("mmm_multiply         ", times[18],counter,times[0]);
       print_timer("queue fft            ", times[30],counter,times[0]);
       print_timer("queue fft serial     ", times[31],counter,times[0]);
       print_timer("queue fft parallel   ", times[32],counter,times[0]);

       print_timer("projector generate   ", times[70],counter,times[0]);
       print_timer("<P|psi> overlap/mpi  ", times[71],counter,times[0]);
       print_timer("sw1/sw2 generation   ", times[72],counter,times[0]);
       print_timer("psi^t*sw1            ", times[73],counter,times[0]);
   }

};
}

#endif
