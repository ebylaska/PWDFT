
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Parallel.hpp"
#include "iofmt.hpp"
//#include	"control.hpp"
#include "Control2.hpp"
#include "Coulomb12.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "Lattice.hpp"
#include "PGrid.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "exchange_correlation.hpp"
#include "inner_loop_md.hpp"
#include "nwpw_Nose_Hoover.hpp"
#include "nwpw_aimd_running_data.hpp"
#include "psi.hpp"
#include "util_date.hpp"
//#include	"rtdb.hpp"
#include "mpi.h"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

static Parallel *myparallel;
static Lattice *mylattice;
static Control2 *control;
static Pneb *mygrid;
static Ion *myion;
static Strfac *mystrfac;
static Kinetic_Operator *mykin;
static Coulomb12_Operator *mycoulomb12;
static XC_Operator *myxc;
static Pseudopotential *mypsp;
static Ewald *myewald;
static nwpw_Nose_Hoover *mynose;

static bool verlet = false;
static bool SA;
static int version, nfft[3], ne[2], ispin, current_iteration;
static double sa_alpha[2], sa_decay[2];
static double unita[9];
static double cpu1, cpu2, cpu3, cpu4;
static double eave, evar, have, hvar, qave, qvar, eke;

static double *psi0;
static double *psi1;
static double *psi2;
static double *Hpsi;
static double *psi_r;
static double *dn;
static double *hml;
static double *lmbda;
static double *eig;
static double *E;

/******************************************
 *                                        *
 *           ctask_cpmd_start             *
 *                                        *
 ******************************************/
int ctask_cpmd_start(MPI_Comm comm_world0, std::string &rtdbstring,
                     double *rion, double *uion, double *qion, double *fion,
                     double *Etot, double *Eapc, std::ostream &coutput)

{
  myparallel = new Parallel(comm_world0);

  int i, ii, ia, nn, ngrid[3], matype, nelem, icount, done;
  char date[26];
  double sum1, sum2, ev;
  // double *psi0,*psi1,*psi2,*Hpsi,*psi_r;
  // double *dn;
  // double *hml,*lmbda,*eig;
  // double sa_alpha[2],sa_decay[2],Te_init,Tr_init,Te_new,Tr_new;
  double Te_init, Tr_init, Te_new, Tr_new;
  double kb = 3.16679e-6;

  E = new double[60];

  control = new Control2(myparallel->np(), rtdbstring);

  bool hprint = (myparallel->is_master() && control->print_level("high"));
  bool oprint = (myparallel->is_master() && control->print_level("medium"));
  bool lprint = (myparallel->is_master() && control->print_level("low"));

  /* reset Parallel base_stdio_print = lprint */
  myparallel->base_stdio_print = lprint;

  for (ii = 0; ii < 60; ++ii)
    E[ii] = 0.0;

  if (myparallel->is_master())
    seconds(&cpu1);
  if (oprint) {
    std::ios_base::sync_with_stdio();
    coutput
        << "          *****************************************************\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *     Car-Parrinello calculation for molecules,     *\n";
    coutput
        << "          *       microclusters, liquids, and materials       *\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *     [     extended Lagrangian molecular   ]       *\n";
    coutput
        << "          *     [        dynamics simulation          ]       *\n";
    coutput
        << "          *     [          C++ implementation         ]       *\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *            version #7.00   03/20/20               *\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *    This code was developed by Eric J. Bylaska     *\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *****************************************************\n";
    coutput << "          >>> job started at       " << util_date() << " <<<\n";
  }

  mylattice = new Lattice(*control);
  myparallel->init2d(control->np_orbital(), control->pfft3_qsize());

  /* initialize lattice, parallel grid structure */
  psi_get_header(myparallel, &version, nfft, unita, &ispin, ne,
                 control->input_movecs_filename());
  mygrid = new Pneb(myparallel, mylattice, *control, ispin, ne);

  /* initialize psi0, psi1, and psi2 */
  psi0 = mygrid->g_allocate(1);
  psi1 = mygrid->g_allocate(1);
  psi2 = mygrid->g_allocate(1);
  Hpsi = mygrid->g_allocate(1);
  psi_r = mygrid->h_allocate();
  dn = mygrid->r_nalloc(ispin);
  hml = mygrid->m_allocate(-1, 1);
  lmbda = mygrid->m_allocate(-1, 1);
  eig = new double[ne[0] + ne[1]];

  /* read wavefunction */
  psi_read0(mygrid, &version, nfft, unita, &ispin, ne, psi2,
            control->input_movecs_filename());

  /* ortho check */
  sum2 = mygrid->gg_traceall(psi2, psi2);
  sum1 = ne[0] + ne[1];
  if (ispin == 1)
    sum1 *= 2;
  if (std::fabs(sum2 - sum1) > 1.0e-10) {
    if (oprint)
      coutput << "Warning: Gram-Schmidt Being performed on psi2" << std::endl;
  }

  /* read wavefunction velocities */
  mygrid->g_zero(psi1);
  if (psi_filefind(mygrid, control->input_v_movecs_filename()))
    psi_read0(mygrid, &version, nfft, unita, &ispin, ne, psi1,
              control->input_v_movecs_filename());

  /* read in ion structure */
  // Ion myion(myrtdb);
  myion = new Ion(rtdbstring, *control);

  /* setup structure factor */
  mystrfac = new Strfac(myion, mygrid);
  mystrfac->phafac();

  /* initialize operators */
  mykin = new Kinetic_Operator(mygrid);
  mycoulomb12 = new Coulomb12_Operator(mygrid, *control);
  myxc = new XC_Operator(mygrid, *control);

  mypsp = new Pseudopotential(myion, mygrid, mystrfac, *control, coutput);

  /* setup ewald */
  myewald = new Ewald(myparallel, myion, mylattice, *control, mypsp->zv);
  myewald->phafac();

  /* scaling psi velocity */
  mygrid->gg_copy(psi1, psi0);
  mygrid->g_Scale(control->elc_scaling(), psi1);
  double eke0 = control->fake_mass() * mygrid->gg_traceall(psi0, psi0);
  double eke1 = control->fake_mass() * mygrid->gg_traceall(psi1, psi1);

  /* initialize thermostats */
  double w = mykin->ke_ave(psi2);
  mynose = new nwpw_Nose_Hoover(*myion, (mygrid->ne[0] + mygrid->ne[1]), w,
                                *control);

  /* initialize simulated annealing */
  SA = false;
  Te_init = 0.0;
  Tr_init = 0.0;
  sa_alpha[0] = 1.0;
  sa_alpha[1] = 1.0;
  if (control->SA()) {
    sa_decay[0] = control->SA_decay(0);
    sa_decay[1] = control->SA_decay(1);
    if (mynose->on()) {
      SA = true;
      Te_init = mynose->Te;
      Tr_init = mynose->Tr;
    } else {
      double dt = control->time_step();
      SA = false;
      sa_alpha[0] = exp(-(dt / sa_decay[0]));
      sa_alpha[1] = exp(-(dt / sa_decay[1]));
    }
  }

  /* Initialize simulated annealing */

  /* initialize two-electron Gaussian integrals */
  /* initialize paw ncmp*Vloc */

  /* initialize metadynamics and tamd */

  /* initialize dplot */

  /* initialize SIC and HFX */

  /* initialize DFT+U */

  /* initialize META GGA */

  /* initialize vdw */

  /* initialize pressure */

  // input rion, uion
  std::memcpy(myion->rion1, rion, 3 * myion->nion);
  std::memcpy(mypsp->myapc->uion, uion, myion->nion);

  //                 |**************************|
  // *****************   summary of input data  **********************
  //                 |**************************|

  if (oprint) {
    coutput << "\n\n";
    coutput
        << "          ==============  summary of input  ==================\n";
    coutput << "\n input psi filename:  " << control->input_movecs_filename()
            << "\n";
    coutput << " input vpsi filename: " << control->input_v_movecs_filename()
            << "\n";
    coutput << "\n";
    coutput << " number of processors used: " << myparallel->np() << "\n";
    coutput << " processor grid           : " << myparallel->np_i() << " x"
            << myparallel->np_j() << "\n";
    if (mygrid->maptype == 1)
      coutput << " parallel mapping         : 1d-slab"
              << "\n";
    if (mygrid->maptype == 2)
      coutput << " parallel mapping         : 2d-hilbert"
              << "\n";
    if (mygrid->maptype == 3)
      coutput << " parallel mapping         : 2d-hcurve"
              << "\n";
    if (mygrid->isbalanced())
      coutput << " parallel mapping         : balanced"
              << "\n";
    else
      coutput << " parallel mapping         : not balanced"
              << "\n";

    coutput << "\n options:\n";
    coutput << "   ion motion           = ";
    coutput << "yes\n";
    coutput << "   boundary conditions  = ";
    if (control->version == 3)
      coutput << "periodic\n";
    if (control->version == 4)
      coutput << "aperiodic\n";

    coutput << "   electron spin        = ";
    if (ispin == 1)
      coutput << "restricted\n";
    else
      coutput << "unrestricted\n";
    coutput << *myxc;

    coutput << mypsp->print_pspall();

    coutput << "\n total charge =" << Ffmt(8, 3) << control->total_charge()
            << std::endl;

    coutput << "\n atom composition:" << std::endl;
    for (ia = 0; ia < myion->nkatm; ++ia)
      coutput << "   " << myion->atom(ia) << " : " << myion->natm[ia];
    coutput << "\n\n initial ion positions of ions (au):" << std::endl;
    for (ii = 0; ii < myion->nion; ++ii)
      coutput << Ifmt(4) << ii + 1 << " " << myion->symbol(ii) << "\t( "
              << Ffmt(10, 5) << myion->rion(0, ii) << " " << Ffmt(10, 5)
              << myion->rion(1, ii) << " " << Ffmt(10, 5) << myion->rion(2, ii)
              << " ) - atomic mass = " << Ffmt(6, 3) << myion->amu(ii)
              << std::endl;
    coutput << "   G.C.\t( " << Ffmt(10, 5) << myion->gc(0) << " "
            << Ffmt(10, 5) << myion->gc(1) << " " << Ffmt(10, 5) << myion->gc(2)
            << " )" << std::endl;
    coutput << " C.O.M.\t( " << Ffmt(10, 5) << myion->com(0) << " "
            << Ffmt(10, 5) << myion->com(1) << " " << Ffmt(10, 5)
            << myion->com(2) << " )" << std::endl;

    coutput << "\n\n initial velocity of ions (au):" << std::endl;
    for (ii = 0; ii < myion->nion; ++ii)
      coutput << Ifmt(4) << ii + 1 << " " << myion->symbol(ii) << "\t( "
              << Ffmt(10, 5) << myion->vion(0, ii) << " " << Ffmt(10, 5)
              << myion->vion(1, ii) << " " << Ffmt(10, 5) << myion->vion(2, ii)
              << " )" << std::endl;
    coutput << "   G.C.\t( " << Ffmt(10, 5) << myion->vgc(0) << " "
            << Ffmt(10, 5) << myion->vgc(1) << " " << Ffmt(10, 5)
            << myion->vgc(2) << " )" << std::endl;
    coutput << " C.O.M.\t( " << Ffmt(10, 5) << myion->vcom(0) << " "
            << Ffmt(10, 5) << myion->vcom(1) << " " << Ffmt(10, 5)
            << myion->vcom(2) << " )" << std::endl;
    coutput << " number of constraints = " << Ifmt(5) << 0
            << " ( DOF = " << Ifmt(6) << myion->ndof() << " )" << std::endl;
    coutput << std::endl;

    coutput << " number of electrons: spin up =" << Ifmt(6) << mygrid->ne[0]
            << " (" << Ifmt(4) << mygrid->neq[0]
            << " per task) down =" << Ifmt(6) << mygrid->ne[ispin - 1] << " ("
            << Ifmt(4) << mygrid->neq[ispin - 1] << " per task)" << std::endl;

    coutput << std::endl;
    coutput << " supercell:" << std::endl;
    coutput << "      volume = " << Ffmt(10, 2) << mylattice->omega()
            << std::endl;
    coutput << "      lattice:    a1 = < " << Ffmt(8, 3)
            << mylattice->unita(0, 0) << " " << Ffmt(8, 3)
            << mylattice->unita(1, 0) << " " << Ffmt(8, 3)
            << mylattice->unita(2, 0) << " >\n";
    coutput << "                  a2 = < " << Ffmt(8, 3)
            << mylattice->unita(0, 1) << " " << Ffmt(8, 3)
            << mylattice->unita(1, 1) << " " << Ffmt(8, 3)
            << mylattice->unita(2, 1) << " >\n";
    coutput << "                  a3 = < " << Ffmt(8, 3)
            << mylattice->unita(0, 2) << " " << Ffmt(8, 3)
            << mylattice->unita(1, 2) << " " << Ffmt(8, 3)
            << mylattice->unita(2, 2) << " >\n";
    coutput << "      reciprocal: b1 = < " << Ffmt(8, 3)
            << mylattice->unitg(0, 0) << " " << Ffmt(8, 3)
            << mylattice->unitg(1, 0) << " " << Ffmt(8, 3)
            << mylattice->unitg(2, 0) << " >\n";
    coutput << "                  b2 = < " << Ffmt(8, 3)
            << mylattice->unitg(0, 1) << " " << Ffmt(8, 3)
            << mylattice->unitg(1, 1) << " " << Ffmt(8, 3)
            << mylattice->unitg(2, 1) << " >\n";
    coutput << "                  b3 = < " << Ffmt(8, 3)
            << mylattice->unitg(0, 2) << " " << Ffmt(8, 3)
            << mylattice->unitg(1, 2) << " " << Ffmt(8, 3)
            << mylattice->unitg(2, 2) << " >\n";

    {
      double aa1, bb1, cc1, alpha1, beta1, gamma1;
      mylattice->abc_abg(&aa1, &bb1, &cc1, &alpha1, &beta1, &gamma1);
      coutput << "      lattice:    a =    " << Ffmt(8, 3) << aa1
              << " b =   " << Ffmt(8, 3) << bb1 << " c =    " << Ffmt(8, 3)
              << cc1 << std::endl;
      coutput << "                  alpha =" << Ffmt(8, 3) << alpha1
              << " beta =" << Ffmt(8, 3) << beta1 << " gamma =" << Ffmt(8, 3)
              << gamma1 << std::endl;
    }
    coutput << "      density cutoff =" << Ffmt(7, 3) << mylattice->ecut()
            << " fft =" << Ifmt(4) << mygrid->nx << " x " << Ifmt(4)
            << mygrid->ny << " x " << Ifmt(4) << mygrid->nz << "  (" << Ifmt(8)
            << mygrid->npack_all(0) << " waves " << Ifmt(8) << mygrid->npack(0)
            << " per task)" << std::endl;
    coutput << "      wavefnc cutoff =" << Ffmt(7, 3) << mylattice->wcut()
            << " fft =" << Ifmt(4) << mygrid->nx << " x " << Ifmt(4)
            << mygrid->ny << " x " << Ifmt(4) << mygrid->nz << "  (" << Ifmt(8)
            << mygrid->npack_all(1) << " waves " << Ifmt(8) << mygrid->npack(1)
            << " per task)" << std::endl;
    coutput << "\n";
    coutput << " Ewald parameters:\n";
    coutput << "      energy cutoff = " << Ffmt(7, 3) << myewald->ecut()
            << " fft =" << Ifmt(4) << myewald->nx() << " x " << Ifmt(4)
            << myewald->ny() << " x " << Ifmt(4) << myewald->nz() << "  ("
            << Ifmt(8) << myewald->npack_all() << " waves " << Ifmt(8)
            << myewald->npack() << " per task)" << std::endl;
    coutput << "      Ewald summation: cut radius = " << Ffmt(7, 3)
            << myewald->rcut() << " and " << Ifmt(3) << myewald->ncut()
            << std::endl;
    coutput << "                       Mandelung Wigner-Seitz =" << Ffmt(12, 8)
            << myewald->mandelung() << " (alpha=" << Ffmt(12, 8)
            << myewald->rsalpha() << " rs =" << Ffmt(12, 8) << myewald->rs()
            << ")" << std::endl;

    coutput << std::endl;
    coutput << " technical parameters:" << std::endl;
    if (myion->fix_translation)
      coutput << "      translation constrained" << std::endl;
    if (myion->fix_rotation)
      coutput << "      rotation constrained" << std::endl;
    coutput << "      time step =" << Ffmt(11, 2) << control->time_step()
            << " ficticious mass =" << Ffmt(11, 2) << control->fake_mass()
            << std::endl;

    coutput << "      max iterations = " << Ifmt(10)
            << control->loop(0) * control->loop(1) << " (" << Ifmt(5)
            << control->loop(0) << " inner " << Ifmt(5) << control->loop(1)
            << " outer)" << std::endl;
    coutput << std::endl;
    coutput << " velocity scaling: " << std::endl;
    coutput << "      cooling/heating rates  =" << Efmt(12, 5)
            << control->elc_scaling() << " (psi) " << Efmt(12, 5)
            << control->ion_scaling() << " (ion)" << std::endl;
    coutput << "      initial kinetic energy =" << Efmt(12, 5) << eke0
            << " (psi) " << Efmt(12, 5) << myion->eki0 << " (ion)" << std::endl;
    coutput << "                                                 "
            << Efmt(12, 5) << myion->ekg << " (C.O.M.)" << std::endl;
    coutput << "      after scaling          =" << Efmt(12, 5) << eke1
            << " (psi) " << Efmt(12, 5) << myion->eki1 << " (ion)" << std::endl;
    coutput << "      increased energy       =" << Efmt(12, 5) << eke1 - eke0
            << " (psi) " << Efmt(12, 5) << myion->eki1 - myion->eki0 << " (ion)"
            << std::endl;
    coutput << std::endl;

    if (mynose->on())
      coutput << mynose->inputprint();
    else
      coutput << " constant energy simulation" << std::endl;

    if (SA)
      coutput << "      SA decay rate = " << Efmt(10, 3) << sa_decay[0]
              << "  (elc) " << Efmt(10, 3) << sa_decay[1] << " (ion)"
              << std::endl;

    coutput << std::endl << std::endl;
  }

  //                 |**************************|
  // *****************     start iterations     **********************
  //                 |**************************|

  current_iteration = 1;

  if (myparallel->is_master())
    seconds(&cpu2);
  if (oprint) {
    coutput << "         ================ Car-Parrinello iteration "
               "================\n";
    coutput << "     >>> iteration started at " << util_date() << " <<<\n";
    coutput << "     iter.          KE+Energy             Energy        KE_psi "
               "       KE_Ion   Temperature\n";
    coutput << "     "
               "---------------------------------------------------------------"
               "----------------------\n";
  }

  // Newton step - first step using velocity
  bool verlet = false;
  inner_loop_md(verlet, sa_alpha, *control, mygrid, myion, mynose, mykin,
                mycoulomb12, myxc, mypsp, mystrfac, myewald, psi0, psi1, psi2,
                Hpsi, psi_r, dn, hml, lmbda, 1, E);

  if (oprint) {
    coutput << Ifmt(10) << current_iteration << Efmt(19, 10) << E[0]
            << Efmt(19, 10) << E[1] << Efmt(14, 5) << E[2] << Efmt(14, 5)
            << E[3] << Ffmt(14, 2) << myion->Temperature() << std::endl;
  }

  // output fion, qion, Etot, and Eapc
  std::memcpy(fion, myion->fion1, 3 * myion->nion);
  std::memcpy(qion, mypsp->myapc->qion, myion->nion);

  *Etot = E[1] + E[2];
  *Eapc = mypsp->myapc->Eapc;

  return 0;
}

/*********************************************
 *                                           *
 *            ctask_cpmd_run                 *
 *                                           *
 *********************************************/
int ctask_cpmd_run(MPI_Comm comm_world0, double *rion, double *uion,
                   double *qion, double *fion, double *Etot, double *Eapc,
                   std::ostream &coutput) {
  ++current_iteration;
  // input rion, uion
  std::memcpy(myion->rion1, rion, 3 * myion->nion);
  std::memcpy(mypsp->myapc->uion, uion, myion->nion);

  bool verlet = true;
  inner_loop_md(verlet, sa_alpha, *control, mygrid, myion, mynose, mykin,
                mycoulomb12, myxc, mypsp, mystrfac, myewald, psi0, psi1, psi2,
                Hpsi, psi_r, dn, hml, lmbda, 1, E);

  bool oprint = (myparallel->is_master() && control->print_level("medium"));
  if (oprint) {
    coutput << Ifmt(10) << current_iteration << Efmt(19, 10) << E[0]
            << Efmt(19, 10) << E[1] << Efmt(14, 5) << E[2] << Efmt(14, 5)
            << E[3] << Ffmt(14, 2) << myion->Temperature() << std::endl;
  }

  // output fion, qion, Etot, and Eapc
  std::memcpy(fion, myion->fion1, 3 * myion->nion);
  std::memcpy(qion, mypsp->myapc->qion, myion->nion);

  *Etot = E[1] + E[2];
  *Eapc = mypsp->myapc->Eapc;

  return 0;
}

/*********************************************
 *                                           *
 *            ctask_cpmd_stop                *
 *                                           *
 *********************************************/

/* deallocate memory */
int ctask_cpmd_stop(MPI_Comm comm_world0, std::ostream &coutput) {
  /* write wavefunction and velocity wavefunction */
  psi_write(mygrid, &version, nfft, unita, &ispin, ne, psi2,
            control->output_movecs_filename(), coutput);
  psi_write(mygrid, &version, nfft, unita, &ispin, ne, psi0,
            control->output_v_movecs_filename(), coutput);

  /* deallocate memory */
  mygrid->g_deallocate(psi0);
  mygrid->g_deallocate(psi1);
  mygrid->g_deallocate(psi2);
  mygrid->g_deallocate(Hpsi);
  mygrid->h_deallocate(psi_r);
  mygrid->r_dealloc(dn);
  mygrid->m_deallocate(hml);
  mygrid->m_deallocate(lmbda);
  delete[] eig;

  delete mystrfac;
  delete mykin;
  delete mycoulomb12;
  delete myxc;
  delete mypsp;
  delete myewald;
  delete mynose;
  delete myion;

  delete control;
  delete mygrid;
  delete mylattice;
  delete myparallel;

  MPI_Barrier(comm_world0);
  return 0;
}

} // namespace pwdft
