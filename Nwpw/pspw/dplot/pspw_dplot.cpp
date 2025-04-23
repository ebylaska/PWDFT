
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
//
#include "Parallel.hpp"
#include "iofmt.hpp"
#include "parsestring.hpp"
#include "util_linesearch.hpp"
//#include	"control.hpp"
#include "Control2.hpp"
#include "Electron.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "Lattice.hpp"
#include "PGrid.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "exchange_correlation.hpp"
#include "psi.hpp"
#include "util_date.hpp"
//#include	"rtdb.hpp"
#include "mpi.h"

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "psp_file_check.hpp"
#include "psp_library.hpp"

#include "nwpw_dplot.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

/******************************************
 *                                        *
 *              pspw_dplot                *
 *                                        *
 ******************************************/
int pspw_dplot(MPI_Comm comm_world0, std::string &rtdbstring,
               std::ostream &coutput) {
  Parallel myparallel(comm_world0);

  int ne[2],nextra[2],ispin;
  char date[26];
  double *psi1, *psi_r, *dn, *rho;

  Control2 control(myparallel.np(), rtdbstring);
  int flag = control.task();

  bool hprint = (myparallel.is_master() && control.print_level("high"));
  bool oprint = (myparallel.is_master() && control.print_level("medium"));
  bool lprint = (myparallel.is_master() && control.print_level("low"));

  /* reset Parallel base_stdio_print = lprint */
  myparallel.base_stdio_print = lprint;

  if (oprint) {
    std::ios_base::sync_with_stdio();
    coutput
        << "          *****************************************************\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *                 PWDFT PSPW dplot                  *\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *  [              C++ implementation             ]  *\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *              version #1.00   03/15/23             *\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *    This code was developed by Eric J. Bylaska,    *\n";
    coutput
        << "          *    Abhishek Bagusetty, David H. Bross, ...        *\n";
    coutput
        << "          *                                                   *\n";
    coutput
        << "          *****************************************************\n";
    coutput << "          >>> job started at       " << util_date() << " <<<\n";
  }

  /* initialize processor grid structure */
  myparallel.init2d(control.np_orbital(), control.pfft3_qsize());

  /* initialize lattice */
  Lattice mylattice(control);

  /* read in ion structure */
  // Ion myion(myrtdb);
  Ion myion(rtdbstring, control);

  /* Check for and generate psp files                       */
  /* - this routine also sets the valence charges in myion, */
  /*   and total_ion_charge and ne in control               */
  psp_file_check(&myparallel, &myion, control, coutput);
  MPI_Barrier(comm_world0);

   bool fractional = control.fractional();
   if (fractional)
   {
      nextra[0] = control.fractional_orbitals(0);
      if (control.ispin()==2)
         nextra[1] = control.fractional_orbitals(1);
      else
         nextra[1] = 0;
   }
   else
   {
      nextra[0] = 0;
      nextra[1] = 0;
   }
   ispin = control.ispin();
   ne[0] = control.ne(0) + nextra[0];
   ne[1] = control.ne(1) + nextra[1];


  /* initialize parallel grid structure */
  Pneb mygrid(&myparallel, &mylattice, control, ispin,ne);

  /* initialize psi1 */
  int neall = mygrid.neq[0] + mygrid.neq[1];
  int shift1 = 2 * (mygrid.npack(1));
  int shift2 = (mygrid.n2ft3d);
  int n2ft3d = (mygrid.n2ft3d);

  double omega = mylattice.omega();
  double scal1 = 1.0 / ((double)((mygrid.nx) * (mygrid.ny) * (mygrid.nz)));
  double scal2 = 1.0 / omega;
  double dv = omega * scal1;
  //ne[0] = control.ne(0);
  //ne[1] = control.ne(1);
  psi1 = mygrid.g_allocate(1);
  psi_r = mygrid.h_allocate();
  dn = mygrid.r_nalloc(ispin);
  rho = mygrid.r_alloc();
  mygrid.d3db::mygdevice.psi_alloc(mygrid.npack(1),mygrid.neq[0]+mygrid.neq[1],
                                   control.tile_factor());

  bool newpsi = psi_read(&mygrid, control.input_movecs_filename(),
                         control.input_movecs_initialize(), psi1, coutput);
  MPI_Barrier(comm_world0);

  /* initialize nwpw_dplot */
  nwpw_dplot mydplot(&myion, &mygrid, control);

  //                 |**************************|
  // *****************   summary of input data  **********************
  //                 |**************************|
  if (oprint) {
    coutput << "\n";
    coutput
        << "          ==============  summary of input  ==================\n";
    coutput << "\n input psi filename: " << control.input_movecs_filename()
            << "\n";
    coutput << "\n";
    coutput << " number of processors used: " << myparallel.np() << "\n";
    coutput << " processor grid           : " << myparallel.np_i() << " x"
            << myparallel.np_j() << "\n";
    if (mygrid.maptype == 1)
      coutput << " parallel mapping         : 1d-slab"
              << "\n";
    if (mygrid.maptype == 2)
      coutput << " parallel mapping         : 2d-hilbert"
              << "\n";
    if (mygrid.maptype == 3)
      coutput << " parallel mapping         : 2d-hcurve"
              << "\n";
    if (mygrid.isbalanced())
      coutput << " parallel mapping         : balanced"
              << "\n";
    else
      coutput << " parallel mapping         : not balanced"
              << "\n";

    coutput << "\n options:\n";
    coutput << "   boundary conditions  = ";
    if (control.version == 3)
      coutput << "periodic\n";
    if (control.version == 4)
      coutput << "aperiodic\n";

    coutput << std::endl;
    coutput << " number of orbitals: spin up =" << Ifmt(6) << mygrid.ne[0]
            << " (" << Ifmt(4) << mygrid.neq[0]
            << " per task) down =" << Ifmt(6) << mygrid.ne[ispin - 1] << " ("
            << Ifmt(4) << mygrid.neq[ispin - 1] << " per task)" << std::endl;

    coutput << std::endl;
    coutput << " ncell              = " << Ifmt(2) << mydplot.ncell[0]
            << Ifmt(2) << mydplot.ncell[1] << Ifmt(2) << mydplot.ncell[2]
            << std::endl;
    coutput << " position tolerance = " << Efmt(12, 6)
            << mydplot.position_tolerance << std::endl;
    coutput << "             origin = " << Ffmt(8, 3) << mydplot.origin[0]
            << Ffmt(8, 3) << mydplot.origin[1] << Ffmt(8, 3)
            << mydplot.origin[2] << std::endl;

    coutput << std::endl;
    coutput << " supercell:" << std::endl;
    coutput << "      volume = " << Ffmt(10, 2) << mylattice.omega()
            << std::endl;
    coutput << "      lattice:    a1 = < " << Ffmt(8, 3)
            << mylattice.unita(0, 0) << " " << Ffmt(8, 3)
            << mylattice.unita(1, 0) << " " << Ffmt(8, 3)
            << mylattice.unita(2, 0) << " >\n";
    coutput << "                  a2 = < " << Ffmt(8, 3)
            << mylattice.unita(0, 1) << " " << Ffmt(8, 3)
            << mylattice.unita(1, 1) << " " << Ffmt(8, 3)
            << mylattice.unita(2, 1) << " >\n";
    coutput << "                  a3 = < " << Ffmt(8, 3)
            << mylattice.unita(0, 2) << " " << Ffmt(8, 3)
            << mylattice.unita(1, 2) << " " << Ffmt(8, 3)
            << mylattice.unita(2, 2) << " >\n";
    coutput << "      reciprocal: b1 = < " << Ffmt(8, 3)
            << mylattice.unitg(0, 0) << " " << Ffmt(8, 3)
            << mylattice.unitg(1, 0) << " " << Ffmt(8, 3)
            << mylattice.unitg(2, 0) << " >\n";
    coutput << "                  b2 = < " << Ffmt(8, 3)
            << mylattice.unitg(0, 1) << " " << Ffmt(8, 3)
            << mylattice.unitg(1, 1) << " " << Ffmt(8, 3)
            << mylattice.unitg(2, 1) << " >\n";
    coutput << "                  b3 = < " << Ffmt(8, 3)
            << mylattice.unitg(0, 2) << " " << Ffmt(8, 3)
            << mylattice.unitg(1, 2) << " " << Ffmt(8, 3)
            << mylattice.unitg(2, 2) << " >\n";

    {
      double aa1, bb1, cc1, alpha1, beta1, gamma1;
      mylattice.abc_abg(&aa1, &bb1, &cc1, &alpha1, &beta1, &gamma1);
      coutput << "      lattice:    a =    " << Ffmt(8, 3) << aa1
              << " b =   " << Ffmt(8, 3) << bb1 << " c =    " << Ffmt(8, 3)
              << cc1 << std::endl;
      coutput << "                  alpha =" << Ffmt(8, 3) << alpha1
              << " beta =" << Ffmt(8, 3) << beta1 << " gamma =" << Ffmt(8, 3)
              << gamma1 << std::endl;
    }
    coutput << "      density cutoff =" << Ffmt(7, 3) << mylattice.ecut()
            << " fft =" << Ifmt(4) << mygrid.nx << " x " << Ifmt(4) << mygrid.ny
            << " x " << Ifmt(4) << mygrid.nz << "  (" << Ifmt(8)
            << mygrid.npack_all(0) << " waves " << Ifmt(8) << mygrid.npack(0)
            << " per task)" << std::endl;
    coutput << "      wavefnc cutoff =" << Ffmt(7, 3) << mylattice.wcut()
            << " fft =" << Ifmt(4) << mygrid.nx << " x " << Ifmt(4) << mygrid.ny
            << " x " << Ifmt(4) << mygrid.nz << "  (" << Ifmt(8)
            << mygrid.npack_all(1) << " waves " << Ifmt(8) << mygrid.npack(1)
            << " per task)" << std::endl;
    coutput << std::endl;
  }

  /* translate system if origin is not zero */
  if ((mydplot.origin[0] != 0.0) && (mydplot.origin[1] != 0.0) &&
      (mydplot.origin[2] != 0.0)) {
    if (oprint)
      coutput << "...translating origin..." << std::endl
              << std::endl
              << std::endl;
    /*
    auto trans[3] = { -mydplot.origin[0],
                       -mydplot.origin[1],
                       -mydplot.origin[2]};
    myion.translate(trans);

    Strfac mystrfac(&myion, &mygrid);
    mystrfac.phafac();
    call psi_translate(origin, npack1,(ne(1)+ne(2)), dcpl_mb(psi2(1)))
    */
  }

  /* convert psi(G) to psi(r) - Expensive */
  int indx1 = 0;
  int indx2 = 0;
  for (auto i = 0; i < neall; ++i) {
    mygrid.cc_pack_copy(1, psi1 + indx1, psi_r + indx2);
    mygrid.c_unpack(1, psi_r + indx2);
    mygrid.cr_fft3d(psi_r + indx2);

    indx1 += shift1;
    indx2 += shift2;
  }

  /* generate dn */
  mygrid.hr_aSumSqr(scal2, psi_r, dn);

  // std::cout << "number of cubefiles = " << control.number_cubefiles() <<
  // std::endl;
  if (oprint)
    coutput << std::endl;
  for (auto i = 0; i < control.number_cubefiles(); ++i) {
    int cubetype = control.cubetype_cubefiles(i);
    if (cubetype != 0) {
      std::string cubename = control.cubename_cubefiles(i);
      std::string cubekey = control.cubekey_cubefiles(i);
      std::string cube_comment;
      // orbital
      if (cubetype > 0) {
        if (oprint)
          coutput << " generating cubefile - orbital     " << Ifmt(5)
                  << cubetype << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;
        int ishift = (cubetype - 1) * n2ft3d;
        ;
        std::memcpy(rho, psi_r + ishift, n2ft3d * sizeof(double));
        cube_comment = "SCF Molecular Orbitals";
      }
      // orbital2
      else if (cubetype == -99) {
        int i = std::stoi(
            mystring_split(mystring_split(cubekey, "orbital2-")[1], "-")[0]);
        int j = std::stoi(
            mystring_split(mystring_split(cubekey, "orbital2-")[1], "-")[1]);
        if (oprint)
          coutput << " generating cubefile - orbital2    " << Ifmt(5) << i
                  << Ifmt(5) << j << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;
        mygrid.rrr_Mul(psi_r + (i - 1) * n2ft3d, psi_r + (j - 1) * n2ft3d, rho);
        cube_comment = "SCF Molecular Orbitals squared";
      }
      // total density
      else if (cubetype == -1) {
        if (oprint)
          coutput << " generating cubefile - total density    "
                  << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;
        mygrid.rrr_Sum(dn, dn + (ispin - 1) * n2ft3d, rho);
        cube_comment = "SCF Total Density";
      }
      // difference density
      else if (cubetype == -2) {
        if (oprint)
          coutput << " generating cubefile - diff density     "
                  << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;
        std::memcpy(rho, dn, n2ft3d * sizeof(double));
        mygrid.rrr_Minus(dn, dn + (ispin - 1) * n2ft3d, rho);
        cube_comment = "SCF Spin Density";

      }
      // alpha density
      else if (cubetype == -3) {
        if (oprint)
          coutput << " generating cubefile - alpha density    "
                  << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;
        std::memcpy(rho, dn, n2ft3d * sizeof(double));
        cube_comment = "SCF Alpha Density";
      }
      // beta density
      else if (cubetype == -4) {
        if (oprint)
          coutput << " generating cubefile - beta density     "
                  << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;
        std::memcpy(rho, dn + (ispin - 1) * n2ft3d, n2ft3d * sizeof(double));
        cube_comment = "SCF Beta Density";
      }
      // laplacian density
      else if (cubetype == -5) {
        if (oprint)
          coutput << " generating cubefile - laplacian density"
                  << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;
        mygrid.rrr_Sum(dn, dn + (ispin - 1) * n2ft3d, rho);
        mygrid.r_SMul(scal1, rho);
        mygrid.r_zero_ends(rho);
        mygrid.rc_fft3d(rho);
        mygrid.c_pack(0, rho);
        double *Gx = mygrid.Gpackxyz(0, 0);
        double *Gy = mygrid.Gpackxyz(0, 1);
        double *Gz = mygrid.Gpackxyz(0, 2);
        int kk = 0;
        for (auto k = 0; k < mygrid.npack(0); ++k) {
          rho[kk] *= -(Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k]);
          rho[kk + 1] *= -(Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k]);
          kk += 2;
        }
        mygrid.c_unpack(0, rho);
        mygrid.cr_fft3d(rho);
        mygrid.r_zero_ends(rho);
        cube_comment = "SCF Laplacian Density";
      }
      // electrostatic potential density
      else if (cubetype == -6) {
        if (oprint)
          coutput << " generating cubefile - potential density"
                  << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;

        Coulomb12_Operator mycoulomb12(&mygrid, control);

        /* generate coulomb potential */
        if (control.version == 3) {
          double *vc = mygrid.c_pack_allocate(0);
          double *tmp = mygrid.r_alloc();

          /* generate dng */
          mygrid.rrr_Sum(dn, dn + (ispin - 1) * n2ft3d, tmp);
          mygrid.r_zero_ends(tmp);
          mygrid.rc_fft3d(tmp);
          mygrid.c_pack(0, tmp);
          mygrid.c_pack_SMul(0, scal1, tmp);

          mygrid.c_pack_SMul(0, scal2, tmp);

          mycoulomb12.mycoulomb1->vcoulomb(tmp, vc);

          mygrid.cc_pack_copy(0, vc, tmp);
          mygrid.c_unpack(0, tmp);
          mygrid.cr_fft3d(tmp);
          mygrid.r_zero_ends(tmp);

          mygrid.rr_copy(tmp, rho);

          mygrid.r_dealloc(tmp);
          mygrid.c_pack_deallocate(vc);

        } else if (control.version == 4) {
          double *tmp = mygrid.r_alloc();

          mygrid.rrr_Sum(dn, dn + (ispin - 1) * n2ft3d, tmp);
          mycoulomb12.mycoulomb2->vcoulomb(tmp, rho);

          mygrid.r_dealloc(tmp);
        }
        cube_comment = "SCF Electrostatic Potential";
      }
      // full electrostatic potential density
      else if (cubetype == -7) {
        if (oprint)
          coutput << " generating cubefile - full potential density"
                  << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;

        /* setup Coulomb operator and structure factor */
        Strfac mystrfac(&myion, &mygrid);
        mystrfac.phafac();
        Coulomb12_Operator mycoulomb12(&mygrid, control);

        /* define rcut */
        double rcut = 1.0;

        /* generate coulomb potential */
        if (control.version == 3) {
          double *vc = mygrid.c_pack_allocate(0);
          double *vg = mygrid.t_pack_allocate(0);
          double *exi = mygrid.c_pack_allocate(0);
          double *tmp = mygrid.r_alloc();
          double *Gx = mygrid.Gpackxyz(0, 0);
          double *Gy = mygrid.Gpackxyz(0, 1);
          double *Gz = mygrid.Gpackxyz(0, 2);
          double w = 0.25 * rcut * rcut;

          /* generate dng */
          mygrid.rrr_Sum(dn, dn + (ispin - 1) * n2ft3d, tmp);
          mygrid.r_zero_ends(tmp);
          mygrid.rc_fft3d(tmp);
          mygrid.c_pack(0, tmp);
          mygrid.c_pack_SMul(0, scal1, tmp);

          // add negative ion charges
          for (auto ii = 0; ii < myion.nion; ++ii) {
            auto ia = myion.katm[ii];
            auto scal = -myion.zv_psp[ia] / omega;
            for (auto k = 0; k < mygrid.npack(0); ++k) {
              auto gg = (Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k]);
              vg[k] = scal * std::exp(-w * gg);
            }
            mystrfac.strfac_pack(0, ii, exi);
            mygrid.tc_pack_Mul(0, vg, exi);
            mygrid.cc_pack_Sum2(0, exi, tmp);
          }
          mygrid.c_pack_SMul(0, scal2, tmp);

          mycoulomb12.mycoulomb1->vcoulomb(tmp, vc);

          mygrid.cc_pack_copy(0, vc, tmp);
          mygrid.c_unpack(0, tmp);
          mygrid.cr_fft3d(tmp);
          mygrid.r_zero_ends(tmp);

          mygrid.rr_copy(tmp, rho);

          mygrid.r_dealloc(tmp);
          mygrid.t_pack_deallocate(vg);
          mygrid.c_pack_deallocate(exi);
          mygrid.c_pack_deallocate(vc);

        } else if (control.version == 4) {
          double sqrt_pi = std::sqrt(4.0 * std::atan(1.0));
          double c = 1.00 / rcut;

          mygrid.initialize_r_grid();
          double *r_grid = mygrid.r_grid;
          double *tmp = mygrid.r_alloc();

          mygrid.rrr_Sum(dn, dn + (ispin - 1) * n2ft3d, tmp);
          mycoulomb12.mycoulomb2->vcoulomb(tmp, rho);

          for (auto ii = 0; ii < myion.nion; ++ii) {
            auto ia = myion.katm[ii];
            auto q = -myion.zv_psp[ia];
            auto xii = myion.rion(0, ii);
            auto yii = myion.rion(1, ii);
            auto zii = myion.rion(2, ii);
            for (auto i = 0; i < n2ft3d; ++i) {
              auto x = (r_grid[3 * i] - xii);
              auto y = (r_grid[3 * i + 1] - yii);
              auto z = (r_grid[3 * i + 2] - zii);
              auto r = std::sqrt(x * x + y * y + z * z);

              rho[i] += ((r > 1.0e-15) ? (q / r) * std::erf(r * c)
                                       : 2.0 * q * c / sqrt_pi);
            }
          }
          mygrid.r_dealloc(tmp);
        }
        cube_comment = "SCF Full Electrostatic Potential";
      }

      // ELF density
      else if (cubetype == -8) {
        if (oprint)
          coutput << " generating cubefile - ELF density"
                  << " - cubefilename = " << cubename
                  << " (json key=" << cubekey << ")" << std::endl;
        cube_comment = "SCF ELF";
      }

      mydplot.gcube_write(cubename, cubetype, cube_comment, rho);
    }
  }

  //                 |**************************|
  // *****************   report consumed time   **********************
  //                 |**************************|

  /* deallocate memory */
  mygrid.g_deallocate(psi1);
  mygrid.h_deallocate(psi_r);
  mygrid.r_dealloc(dn);
  mygrid.r_dealloc(rho);
  mygrid.d3db::mygdevice.psi_dealloc();

  MPI_Barrier(comm_world0);

  return 0;
}

} // namespace pwdft
