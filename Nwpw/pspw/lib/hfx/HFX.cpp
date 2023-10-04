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
#include "filon_filter.hpp"
#include "HFX.hpp"

#include "parsestring.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *         coulomb_screened_kernel         *
 *                                         *
 *******************************************/
//extern void coulomb_filter(Pneb *, double *, const std::string);
static void coulomb_screened_kernel(Pneb *mygrid,
                                   const int screening_type, const double rcut, const double pp, 
                                   const double attenuation,
                                   const std::string kernel_filter_filename,
                                   double *vg)
{
   double epsilon = 1.0;
   if (screening_type==2) epsilon=attenuation;

   int nxh = (mygrid->nx/2);
   int nyh = (mygrid->ny/2);
   int nzh = (mygrid->nz/2);

   double pp2 = pp + 2.0;
   double omega = mygrid->lattice->omega();
   double scal1 = 1.0/((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
   double dv = omega * scal1;
   double pi = 4.0*atan(1.0);
   double fourpi = 4*pi;
   double sqrt_pi = std::sqrt(pi);
   double eps = 1.0e-12;
   double *Gx = mygrid->Gpackxyz(0,0);
   double *Gy = mygrid->Gpackxyz(0,1);
   double *Gz = mygrid->Gpackxyz(0,2);
   double gg;

   mygrid->initialize_r_grid();
   double *r_grid = mygrid->r_grid;

   int taskid = mygrid->d3db::parall->taskid_i();
   int pzero  = mygrid->ijktop(0, 0, 0);
   int zero   = mygrid->ijktoindex(0, 0, 0);

   int npack0 = mygrid->npack(0);
   int n2ft3d = mygrid->n2ft3d;
   int nfft3d = mygrid->nfft3d;
   double *glr = new double[n2ft3d];
   double *gk  = new double[2*nfft3d];
   double *tmp = new double[nfft3d];

   // Use aperiodic definitions of kernel 
   if ((screening_type==0) || (screening_type==2))
   {
      double *gk = new double[2*nfft3d];
      double *gr = new double[n2ft3d];
      double epsilon2 = epsilon*epsilon;
 
      // short-range part of Greens function set only for short-range
      std::memset(gk,0,2*nfft3d*sizeof(double));
      for (auto k=0; k<nfft3d; ++k) 
      {
         gg = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
         if ((pzero == taskid) && (k == zero))
            gk[2*k] = pi/epsilon2;
         else
            gk[2*k] = (fourpi / gg)*(1.0 - std::exp(-gg/(4.0*epsilon2)));
      }
      mygrid->cr_fft3d(gk);
      mygrid->r_SMul(scal1,gk);

      // long-range part of Greens function
      std::memset(glr,0,n2ft3d*sizeof(double));
      for (auto k=0; k<n2ft3d; ++k)
      {
          double x = r_grid[3*k];
          double y = r_grid[3*k+1];
          double z = r_grid[3*k+2];
          double temp = std::sqrt(x*x + y*y + z*z);
          if (temp>1.0e-10) 
              temp = std::erf((epsilon*temp)/temp);
           else
              temp = 2.0*epsilon/sqrt_pi;
         glr[k] = temp*dv;
      }
      mygrid->rr_Sum(gk,glr);

      // multiply by the screening function ****
      for (auto k=0; k<n2ft3d; ++k)
      {
         double x = r_grid[3*k];
         double y = r_grid[3*k+1];
         double z = r_grid[3*k+2];
         double temp = std::sqrt(x*x + y*y + z*z);
         glr[k] = glr[k]* (1.0 - std::pow((1.0 - std::exp(-std::pow((temp/rcut),pp2))), pp));
      }
      mygrid->r_zero_ends(glr);
      mygrid->rc_fft3d(glr);

      for (auto k=0; k<nfft3d; ++k)
          tmp[k] = glr[2*k];
      mygrid->t_pack(0, tmp);
      mygrid->tt_pack_copy(0, tmp, vg);

   }

   // screening_type == 1 use periodic definitions of kernel
   else if (screening_type==1)
   {
      std::memset(vg,0,npack0*sizeof(double));
      for (auto k=0; k<npack0; ++k) 
      {
         gg = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
         if (gg < eps)
            vg[k] = 5.633714987781071* std::pow(omega,(2.0/3.0)) /pi;
         else
            vg[k] = (fourpi / gg);
      }
   }

   // screening_type == 3 use periodic definitions of cutoff-kernel
   else if (screening_type==3)
   {
      std::memset(vg,0,npack0*sizeof(double));
      for (auto k=0; k<npack0; ++k) 
      {
         gg = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
         if (gg < eps)
            vg[k] = util_kiril_coulomb_transform0(3,rcut,pp);
         else
            vg[k] = util_kiril_coulomb_transform(3,gg,rcut,pp);
      }
   }

   // screeing_type == 4 use erfc definitions of kernel 
   else if (screening_type==4)
   {
      double invrcut  = 1.0/rcut;
      double invrcut2 = invrcut*invrcut;
      std::memset(vg,0,npack0*sizeof(double));
      for (auto k=0; k<npack0; ++k)
      {
         gg = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
         if (gg < eps )
            vg[k] = pi/invrcut2;
         else
            vg[k] = (fourpi/gg) * (1.0 - exp(-gg/(4.0*invrcut2)));
      }
   }


   // screeing_type == 5 Filtered exchange definitions of kernel   ******
   else if (screening_type==5)
   {
      std::memset(vg,0,npack0*sizeof(double));
      //coulomb filter here
      coulomb_filter(mygrid,vg,kernel_filter_filename);
   }
   delete[] tmp;
   delete[] glr;
   delete[] gk;
}

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

   kernel_filter_filename = control.permanent_dir_str  + "/filtered_kernel.dat";

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
   if (mystring_contains(mystring_lowercase(xc_name), "-attenuation"))
   {
      attenuation = std::stod(mystring_split(mystring_lowercase(xc_name), "-attenuation")[1]);
   }
   if (mystring_contains(mystring_lowercase(xc_name), "-filter_filename"))
   {
      kernel_filter_filename = std::stod(mystring_split(mystring_lowercase(xc_name), "-filter_filename")[1]);
      if (!mystring_contains(kernel_filter_filename,"/"))
         kernel_filter_filename = control.permanent_dir_str  + kernel_filter_filename;
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


      vg = new double[mygrid->npack(0)];
      
      new_coulomb2 = false;

      // periodic solver
      if (solver_type==0)
      {
         coulomb_screened_kernel(mygrid,screening_type,rcut,pp,attenuation,kernel_filter_filename,vg);
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
