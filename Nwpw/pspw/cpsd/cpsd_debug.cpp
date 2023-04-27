
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;

#include "Parallel.hpp"
//#include	"control.hpp"
#include "Control2.hpp"
#include "Coulomb.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "Lattice.hpp"
#include "PGrid.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "exchange_correlation.hpp"
#include "inner_loop.hpp"
#include "psi.hpp"
#include "util_date.hpp"
//#include	"rtdb.hpp"
#include "mpi.h"

#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "psp_file_check.hpp"
#include "psp_library.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {
using namespace pwdft;

/******************************************
 *                                        *
 *                cpsd_debug              *
 *                                        *
 ******************************************/
int cpsd_debug(MPI_Comm comm_world0, string &rtdbstring) {
  // Parallel myparallel(argc,argv);
  Parallel myparallel(comm_world0);
  // RTDB myrtdb(&myparallel, "eric.db", "old");

  int version, nfft[3], ne[2], ispin;
  int i, ii, ia, nn, ngrid[3], matype, nelem, icount, done;
  char date[26];
  double sum1, sum2, ev, zv;
  double cpu1, cpu2, cpu3, cpu4;
  double E[50], deltae, deltac, deltar, viral, unita[9];
  double *psi1, *psi2, *Hpsi, *psi_r;
  double *dn;
  double *hml, *lmbda, *eig;

  for (ii = 0; ii < 50; ++ii)
    E[ii] = 0.0;

  if (myparallel.is_master()) {
    seconds(&cpu1);
    ios_base::sync_with_stdio();
    cout << "          *****************************************************\n";
    cout << "          *                                                   *\n";
    cout << "          *        Car-Parrinello DEBUG calculation           *\n";
    cout << "          *                                                   *\n";
    cout << "          *     [     steepest descent minimization   ]       *\n";
    cout << "          *     [          C++ implementation         ]       *\n";
    cout << "          *                                                   *\n";
    cout << "          *            version #x.xx   11/11/21               *\n";
    cout << "          *                                                   *\n";
    cout << "          *    This code was developed by Eric J. Bylaska     *\n";
    cout << "          *                                                   *\n";
    cout << "          *****************************************************\n";
    cout << "          >>> job started at       " << util_date() << " <<<\n";
  }
  cout << "HELLO A np = " << myparallel.np() << endl;

  Control2 control(myparallel.np(), rtdbstring);

  cout << "HELLO B np = " << myparallel.np() << endl;

  /* initialize processor grid structure */
  myparallel.init2d(control.np_orbital(), control.pfft3_qsize());
  MPI_Barrier(comm_world0);

  /* initialize lattice */
  Lattice mylattice(control);

  /* read in ion structure */
  // Ion myion(myrtdb);
  Ion myion(rtdbstring, control);
  MPI_Barrier(comm_world0);

  /* Check for and generate psp files                       */
  /* - this routine also sets the valence charges in myion, */
  /*   and total_ion_charge and ne in control               */
  psp_file_check(&myparallel, &myion, control, cout);
  cout << "HELLO F\n";

  MPI_Barrier(comm_world0);

  /* fetch ispin and ne psi information from control */
  ispin = control.ispin();
  ne[0] = control.ne(0);
  ne[1] = control.ne(1);
  nfft[0] = control.ngrid(0);
  nfft[1] = control.ngrid(1);
  nfft[2] = control.ngrid(2);
  unita[0] = mylattice.unita1d(0);
  unita[1] = mylattice.unita1d(1);
  unita[2] = mylattice.unita1d(2);
  unita[3] = mylattice.unita1d(3);
  unita[4] = mylattice.unita1d(4);
  unita[5] = mylattice.unita1d(5);
  unita[6] = mylattice.unita1d(6);
  unita[7] = mylattice.unita1d(7);
  unita[8] = mylattice.unita1d(8);
  version = 3;
  MPI_Barrier(comm_world0);

  cout << "HELLO G ne=" << control.ne_ptr()[0] << " " << control.ne_ptr()[1]
       << std::endl;

  // initialize parallel grid structure
  // MPI_Barrier(comm_world0);
  Pneb mygrid(&myparallel, &mylattice, control, control.ispin(),
              control.ne_ptr());
  MPI_Barrier(comm_world0);

  /*

     cout << "HELLO H\n";

     // initialize psi1 and psi2
     psi1  = mygrid.g_allocate(1);
     psi2  = mygrid.g_allocate(1);
     Hpsi  = mygrid.g_allocate(1);
     psi_r = mygrid.h_allocate();
     dn    = mygrid.r_nalloc(ispin);
     hml   = mygrid.m_allocate(-1,1);
     lmbda = mygrid.m_allocate(-1,1);
     eig   = new double[ne[0]+ne[1]];
     gdevice_psi_alloc(mygrid.npack(1),mygrid.neq[0]+mygrid.neq[1],1);
     MPI_Barrier(comm_world0);

     cout << "HELLO I\n";



     cout << "HELLO J\n";


     // deallocate memory
     MPI_Barrier(comm_world0);
     mygrid.g_deallocate(psi1);
     mygrid.g_deallocate(psi2);
     mygrid.g_deallocate(Hpsi);
     mygrid.h_deallocate(psi_r);
     mygrid.r_dealloc(dn);
     mygrid.m_deallocate(hml);
     mygrid.m_deallocate(lmbda);
     delete [] eig;
     MPI_Barrier(comm_world0);

     gdevice_psi_dealloc();
   */
  MPI_Barrier(comm_world0);

  cout << "HELLO K\n";

  return 0;
}

} // namespace pwdft
