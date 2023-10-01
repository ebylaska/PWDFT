/* HFX.cpp -
   Author - Eric Bylaska
*/

/*
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>
#include	<string>
*/
#include <cmath>
#include "Coulomb2.hpp"
#include "HFX.hpp"

#include "parsestring.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *       HFX_Operator::HFX_Operator        *
 *                                         *
 *******************************************/
HFX_Operator::HFX_Operator(Pneb *mygrid, bool has_coulomb2_in, Coulomb2_Operator *mycoulomb2_in, Control2 &control) 
{
   mypneb = mygrid;

   has_coulomb2 = has_coulomb2_in;
   mycoulomb2   = mycoulomb2_in;

   // parse for exchange-correlation functional that contain exact exchange
   std::string xc_name = control.xc_name();

   if (mystring_contains(mystring_lowercase(xc_name), "bnl")) { hfx_on=true;hfx_parameter=1.00;screening_type= 2;}
   if (mystring_contains(mystring_lowercase(xc_name), "hse")) { hfx_on=true;hfx_parameter=0.25;screening_type= 4;rcut=1.0/0.207; }
   if (mystring_contains(mystring_lowercase(xc_name), "pbe0"))         { hfx_on = true; hfx_parameter = 0.25; }
   if (mystring_contains(mystring_lowercase(xc_name), "blyp0"))        { hfx_on = true; hfx_parameter = 0.25; }
   if (mystring_contains(mystring_lowercase(xc_name), "revpbe0"))      { hfx_on = true; hfx_parameter = 0.25; }
   if (mystring_contains(mystring_lowercase(xc_name), "b3lyp"))        { hfx_on = true; hfx_parameter = 0.20; }
   if (mystring_contains(mystring_lowercase(xc_name), "hartree-fock")) { hfx_on = true; hfx_parameter = 1.00; }
   if (mystring_contains(mystring_lowercase(xc_name), "hf"))           { hfx_on = true; hfx_parameter = 1.00; }

   if (mystring_contains(mystring_lowercase(xc_name), "-unrelaxed")) { relaxed = false; }
   if (mystring_contains(mystring_lowercase(xc_name), "-pert"))      { relaxed = false; }
   if (mystring_contains(mystring_lowercase(xc_name), "-screening_radius"))
   {
      rcut = std::stod(mystring_split(mystring_lowercase(xc_name), "-screening_radius")[1]);
   }
   if (mystring_contains(mystring_lowercase(xc_name), "-screening_power"))
   {
      pp = std::stod(mystring_split(mystring_lowercase(xc_name), "-screening_power")[1]);
   }
   if (mystring_contains(mystring_lowercase(xc_name), "-screening_type"))
   {
      screening_type = std::stoi(mystring_split(mystring_lowercase(xc_name), "-screening_type")[1]);
   }
   if (mystring_contains(mystring_lowercase(xc_name), "-scaling_parameter"))
   {
      hfx_parameter = std::stod(mystring_split(mystring_lowercase(xc_name), "-scaling_parameter")[1]);
   }

   // set the number of hfx orbitals
   ispin = mygrid->ispin;
   for (auto ms=0; ms<ispin; ++ms)
   {
      norbs[ms] = mygrid->ne[ms];
      orbital_list[ms] = new (std::nothrow) int[norbs[ms]]();
      ehfx_orb[ms]     = new (std::nothrow) double[norbs[ms]]();
      for (auto n=0; n<norbs[ms]; ++n)
         orbital_list[ms][n] = n + ms*norbs[0];
   }

   // Set up HFX
   if (hfx_on)
   {
      // default to periodic solver
      if (control.version==3) solver_type = 0; // periodic solver
      if (control.version==4) solver_type = 1; // aperiodic solver
      if (mystring_contains(mystring_lowercase(xc_name), "-periodic"))  { solver_type = 0; }
      if (mystring_contains(mystring_lowercase(xc_name), "-aperiodic")) { solver_type = 1; }

      int k, pzero, zero, taskid;
      double gg;
      double *Gx = mygrid->Gxyz(0);
      double *Gy = mygrid->Gxyz(1);
      double *Gz = mygrid->Gxyz(2);
      vg = new double[mygrid->npack(0)];
      double *tmp = new double[mygrid->nfft3d];
      double fourpi = 16.0 * atan(1.0);

      
      new_coulomb2 = false;

      // periodic solver
      if (solver_type==0)
      {
         //call coulomb_screened_init(flag,rcut,pp)
      }
      // aperiodic solver 
      else if (solver_type==1)
      {
          if (!has_coulomb2)
          {
             new_coulomb2 = true;
             mycoulomb2   = new (std::nothrow) Coulomb2_Operator(mygrid, control);
          }
      }
 

//*        **** initialize coulomb_screened ****
//*        **** initialize free-space coulomb if necessary ****
//*     ***** flag == 1 use periodic definitions of kernel ******
//*     ***** flag==3 use periodic definitions of cutoff-kernel ******
//*     ***** flag==4 use erfc definitions of kernel            ******
//*     ***** flag==5 Filtered exchange definitions of kernel   ******
 
 /*
      taskid = mypneb->d3db::parall->taskid_i();
      pzero = mypneb->ijktop(0, 0, 0);
      zero = mypneb->ijktoindex(0, 0, 0);
 
      for (k = 0; k < (mypneb->nfft3d); ++k) {
        gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
        if ((pzero == taskid) && (k == zero))
          tmp[k] = 0.0;
        else
          tmp[k] = fourpi / gg;
      }
      mypneb->t_pack(0, tmp);
      mypneb->tt_pack_copy(0, tmp, vg);
      */
  
      delete[] tmp;
   }
}

/*******************************************
 *                                         *
 *        HFX_Operator::v_exchange         *
 *                                         *
 *******************************************/
void HFX_Operator::v_exchange(const double *psi_r, double *Hpsi_r) 
{
    /*
   int ksize = (mypneb->npack(0));
   int k1 = 0;
   for (auto k=0; k<ksize; ++k) 
   {
      vcout[k1]   = vg[k] * dng[k1];
      vcout[k1+1] = vg[k] * dng[k1+1];
      k1 += 2;
   }
   */
}


/*******************************************
 *                                         *
 *        HFX_Operator::e_exchange         *
 *                                         *
 *******************************************/
void HFX_Operator::e_exchange(const double *psi_r,double& ehfx_out, double& phfx_out) 
{
    /*
   int ksize1 = (mypneb->nzero(0));
   int ksize2 = (mypneb->npack(0));
   double ave = 0.0;

   int k1 = 0;
   for (auto k=0; k<ksize1; ++k) 
   {
      ave += vg[k] * (dng[k1]*dng[k1] + dng[k1+1]*dng[k1+1]);
      k1 += 2;
   }
   for (auto k=ksize1; k<ksize2; ++k) 
   {
      ave += 2.0 * vg[k] * (dng[k1]*dng[k1] + dng[k1+1]*dng[k1+1]);
      k1  += 2;
   
   ave = mypneb->d3db::parall->SumAll(1,ave);
   ave *= 0.5*(mypneb->lattice->omega());
 
   return ave;
   */
}

} // namespace pwdft
