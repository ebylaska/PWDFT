/* Paw_compcharge.cpp - 
   Author - Eric Bylaska
*/

#include        <iostream>
#include        <cstdio>
#include        <cstdlib>
#include        <cstring>
#include	<complex>
#include        <cmath>

using namespace std;


#include	"blas.h"
#include	"util.hpp"
#include        "util_gaunt.hpp"
#include        "util_wgaussian.hpp"
#include        "Paw_compcharge.hpp"

namespace pwdft {
using namespace pwdft;


/***********************************
 *                                 *
 *  Paw_compcharge_gen_smoothpsp   *
 *                                 *
 ***********************************/
static void Paw_compcharge_gen_smoothpsp(const int nray, const double Gray[], const double vlray[],
                                         const int npack0, const double Gx[], const double Gy[], const double Gz[],
                                         double vlsmooth[])
{
   double dG = Gray[2]-Gray[1];
   for (auto k=0; k<npack0; ++k)
   {
      double Q = sqrt(Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k]);
      if (Q>1.0e-9)
      {
         double nxray = (Q/dG) + 1.0;
         double P = util_splint(&Gray[1],&vlray[1],&vlray[1+nray],nray-1,nxray-1,Q);
         vlsmooth[k]=P;
      }
      else
         vlsmooth[k]=vlray[0];
   }

   std::cout << "out smoothpsp" << std::endl;
}


/***********************************
 *                                 *
 *    Paw_compcharge_gen_vlray     *
 *                                 *
 ***********************************/
static void Paw_compcharge_gen_vlray(const double sigma_smooth,
                                     const double log_amesh, const int nrho, const double rho[], const double f[],
                                     const int nray, const double Gray[], double vlray[], double tmpray[])
{
   double twopi  = 8.0*atan(1.0);
   double fourpi = 2.0*twopi;
   //double ecut   = control_ecut()
   //double  rlocal = 1.0d0

   memset(vlray,0,2*nray*sizeof(double));
   for (auto k=1; k<nray; ++k)
   {
      double q = Gray[k];
      double x = sigma_smooth*q;
      vlray[k] = -(fourpi/(q*q))*exp(-0.25*x*x);
   }

   vlray[0] = 0.50*twopi*sigma_smooth*sigma_smooth;

   //if (control_kbpp_filter())  nwpw_kbpp_filter_ray(nray,Gray,ecut,vlray)

   double dG = Gray[2]-Gray[1];

   //**** five point formula ***
   double yp1 = ( -50.0*vlray[1]
                 + 96.0*vlray[2]
                 - 72.0*vlray[3]
                 + 32.0*vlray[4]
                 -  6.0*vlray[5])/(24.0*dG);

   std::cout << "into util_spline vlray" << std::endl;

   util_spline(&Gray[1],&vlray[1],nray-1,yp1,0.0,&vlray[1+nray],tmpray);
   std::cout << "out util_spline vlray" << std::endl;
}



/*************************************************
 *                                               *
 *          Paw_compcharge_gen_vk_smooth         *
 *                                               *
 *************************************************/
static void Paw_compcharge_gen_vk_smooth(const int nray, const double *Gray,
                                         const double sigma_smooth,
                                         const int npack0, 
                                         const double Gx[], const double Gy[], const double Gz[], 
                                         double vk[])
{

   double vlray[2*nray],tmpray[nray];

   double bmesh     = 1.0050;
   double log_bmesh = log(bmesh);
   int    nrho      = ((int) (log(25.0/0.00025)/log_bmesh)) + 1;

   //*** make sure loggrid is odd ***
   if ((nrho%2)==0) ++nrho;

   double f[nrho],rho[nrho];

   rho[0] = 0.00025;
   for (auto i=1; i<nrho; ++i)
      rho[i] = bmesh*rho[i-1];

   std::cout << "Into gen_vlray, nray=" << nray << std::endl;
      
   Paw_compcharge_gen_vlray(sigma_smooth,log_bmesh,nrho,rho,f,nray,Gray,vlray,tmpray);

   std::cout << "Into gen_smoothpsp" << std::endl;
   Paw_compcharge_gen_smoothpsp(nray,Gray,vlray,npack0,Gx,Gy,Gz,vk);
}



#define matindx2(i,j,l,nbasis)          (i + j*nbasis + l*nbasis*nbasis)
#define matindx4(i1,j1,i2,j2,l,nbasis)  (i1 + j1*nbasis + i2*nbasis*nbasis + j2*nbasis*nbasis*nbasis + l*nbasis*nbasis*nbasis*nbasis)

/*******************************************************
 *                                                     *
 *            Paw_compcharge::Paw_compcharge           *
 *                                                     *
 *******************************************************/

Paw_compcharge::Paw_compcharge(Ion *myion0, Pneb *mypneb0,  Control2& control,
                               const int nprj[], const int nbasis[], const int psp_type[], const int lmax0[],
                               const double sigma[],
                               const int nprj_max, int* l_prj[], int* m_prj[], int* b_prj[],
                               double* comp_charge_matrix[], double* hartree_matrix[])
{
   int ia,iia,ii,iii;

   myion    = myion0;
   mypneb   = mypneb0;


   int nion  = myion->nion;
   int nkatm = myion->nkatm;

   //**** set isgamma and use_grid_cmp routines ****
   isgamma = true;
   use_grid_cmp = control.use_grid_cmp();


   double pi   = 4.0*atan(1.0);
   double rcut = control.ewald_rcut();
   if ((control.version==4) && (rcut<=0.0)) rcut = 1.0;
   if (rcut<=0.0)
   {
      double rs,w;
      rs = std::pow(mypneb->lattice->unita(0,0),2) 
         + std::pow(mypneb->lattice->unita(1,0),2) 
         + std::pow(mypneb->lattice->unita(2,0),2);
      rs = sqrt(rs);
      rcut=rs/pi;

      rs = std::pow(mypneb->lattice->unita(0,1),2) 
         + std::pow(mypneb->lattice->unita(1,1),2) 
         + std::pow(mypneb->lattice->unita(2,1),2);
      rs = sqrt(rs);
      w=rs/pi;
      if (w<rcut) rcut = w;

      rs = std::pow(mypneb->lattice->unita(0,2),2) 
         + std::pow(mypneb->lattice->unita(1,2),2) 
         + std::pow(mypneb->lattice->unita(2,2),2);
      rs = sqrt(rs);
      w=rs/pi;
      if (w<rcut) rcut = w;
   }
   sigma_smooth = rcut;


   //**** determine nion_paw, nkatm_paw, katm_paw, ****
   //****           katm_pawtoion, ion_pawtoion,  ****
   //****           katm_iontopaw, ion_iontopaw     ****
   nion_paw = 0;
   for (ii=0; ii<nion; ++ii)
   {
      ia = myion->katm[ii];
      if (psp_type[ia]==4) ++nion_paw;
   }

   nkatm_paw = 0;
   for (ia=0; ia<nkatm; ++ia)
   {
      if (psp_type[ia]==4) ++nkatm_paw;
   }

   ion_iontopaw  = new int[nion];
   ion_pawtoion  = new int[nion_paw];
   katm_iontopaw = new int[nkatm];
   katm_pawtoion = new int[nkatm_paw];
   katm_paw      = new int[nion_paw];
   sigma_paw     = new double[nkatm_paw];


   iia = 0;
   for (ia=0; ia<nkatm; ++ia)
   {
      if (psp_type[ia]==4) 
      {
         katm_pawtoion[iia] = ia;
         katm_iontopaw[ia]  = iia;
         ++iia;
      }
      else
      {
         katm_iontopaw[ia] = -1;
      }
   }

   iii = 0;
   for (ii=0; ii<nion; ++ii)
   {
      int ia  = myion->katm[ii];
      int iia = katm_iontopaw[ia];
      if (psp_type[ia]==4)
      {
         katm_paw[iii]     = iia;
         ion_pawtoion[iii] = ii;
         ion_iontopaw[ii]  = iii;
         ++iii;
      }
      else
      {
         ion_iontopaw[ii] = -1;
      }
   }

   for (iia=0; iia<nkatm_paw; ++iia)
   {
      ia = katm_pawtoion[iia];
      sigma_paw[iia] = sigma[ia];
   }

   //**** allocate mult_l and lm_size ****
   mult_l  = new int[nkatm_paw];
   lm_size = new int[nkatm_paw];


   mult_l_max = 0;
   for (iia=0; iia<nkatm_paw; ++iia)
   {
      ia = katm_pawtoion[iia];
      mult_l[iia]  = 2*lmax0[ia];
      lm_size[iia] = std::pow((2*lmax0[ia]+1),2);
      if (mult_l_max<(2*lmax0[ia])) mult_l_max = 2*lmax0[ia];
   }

   //*** initialize gaunt
   util_gaunt_init(false,2*mult_l_max);

   //*** allocate gk_smooth, gk,and glm ***
   int npack0  = mypneb->npack(0);
   lm_size_max = (mult_l_max+1)*(mult_l_max+1);

   vk_smooth = new (std::nothrow) double [npack0]();
   gk_smooth = new (std::nothrow) double [npack0]();

   gk  = new (std::nothrow) double [nkatm_paw*  npack0]();
   glm = new (std::nothrow) double [lm_size_max*npack0]();

   Qlm         = new (std::nothrow) double [2*lm_size_max*nion_paw]();
   Qlmx        = new (std::nothrow) double [2*lm_size_max*nion_paw]();
   Qlmy        = new (std::nothrow) double [2*lm_size_max*nion_paw]();
   Qlmz        = new (std::nothrow) double [2*lm_size_max*nion_paw]();
   dEmult_Qlm  = new (std::nothrow) double [2*lm_size_max*nion_paw]();
   dElocal_Qlm = new (std::nothrow) double [2*lm_size_max*nion_paw]();
   dE_Qlm      = new (std::nothrow) double [2*lm_size_max*nion_paw]();

   memset(Qlm,0,2*lm_size_max*nion_paw*sizeof(double));
   memset(Qlmx,0,2*lm_size_max*nion_paw*sizeof(double));
   memset(Qlmy,0,2*lm_size_max*nion_paw*sizeof(double));
   memset(Qlmz,0,2*lm_size_max*nion_paw*sizeof(double));
   memset(dEmult_Qlm,0,2*lm_size_max*nion_paw*sizeof(double));
   memset(dElocal_Qlm,0,2*lm_size_max*nion_paw*sizeof(double));
   memset(dE_Qlm,0,2*lm_size_max*nion_paw*sizeof(double));

   double scal;
   double fourpioveromega = 16.0*atan(1.0)/mypneb->lattice->omega();
   double* Gall[3];
   Gall[0] = mypneb->Gpackxyz(0,0);
   Gall[1] = mypneb->Gpackxyz(0,1);
   Gall[2] = mypneb->Gpackxyz(0,2);

   int zero = 0;
   int one  = 1;
   int lm = 0;
   for (auto l=0; l<=mult_l_max; ++l)
   {
      double phase = 1.0;
      if ((l%4)==0) 
         phase = 1.0;
      else if ((l%4)==1)
         phase = -1.00;
      else if ((l%4)==2)
         phase = -1.0;
      else if ((l%4)==3) 
         phase = 1.0;

      //**** define  |k|**l / (2*l+1)!! ****
      scal = 1.0/util_double_factorial(2*l+1);
      if (l>0)
      {
         for (auto k=0; k<npack0; ++k)
         {
            double gg = Gall[0][k]*Gall[0][k]
                      + Gall[1][k]*Gall[1][k]
                      + Gall[2][k]*Gall[2][k];
            gk[k] = scal*std::pow(sqrt(gg),l);
         }
      }
      else
         DCOPY_PWDFT(npack0,&scal,zero,gk,one);

      //**** define glm = (-i)**l *  |k|**l * Tlm(k)/ (2*l+1)!! ****
      for (auto m=-l; m<=l; ++m)
      {
         util_Tesseral3_vector_lm(l,m,npack0,Gall[0],Gall[1],Gall[2],&glm[lm*npack0]);

         for (auto k=0; k<npack0; ++k)
            glm[lm*npack0+k] = phase*glm[lm*npack0+k]*gk[k];

         ++lm;
      }
   }


   //**** define vk_smooth(k) ****
   int n_ray  = mypneb->n_ray();
   double *G_ray = mypneb->generate_G_ray();

   Paw_compcharge_gen_vk_smooth(n_ray,G_ray,sigma_smooth,npack0,Gall[0],Gall[1],Gall[2],vk_smooth);

   delete [] G_ray;


   //**** define gk_smooth(k)  = 4*pi * Exp[-k*k*sigma_smooth**2 / 4] ****
   scal = 0.25*sigma_smooth*sigma_smooth;
   for (auto k=0; k<npack0; ++k)
   {
      double gg = Gall[0][k]*Gall[0][k]
                + Gall[1][k]*Gall[1][k]
                + Gall[2][k]*Gall[2][k];
      gk_smooth[k] = fourpioveromega*exp(-gg*scal);
   }

   std::cout << "out define gk_smooth" << std::endl;



   //**** define gk(k,iia)  = 4*pi * Exp[-k*k*sigma(iia**2 / 4] ****
   for (auto iia=0; iia<=nkatm_paw; ++iia)
   {
      for (auto k=0; k<npack0; ++k)
      {
         scal = 0.25*sigma_paw[iia]*sigma_paw[iia];
         double gg = Gall[0][k]*Gall[0][k]
                   + Gall[1][k]*Gall[1][k]
                   + Gall[2][k]*Gall[2][k];
         gk[iia*npack0+k] = fourpioveromega*exp(-gg*scal);
      }
   }

   //*******************************************************
   //*****  define indexing for compcharge evalulations ****
   //*******************************************************
   nindx_Tndiff = new (std::nothrow) int [nkatm_paw]();
   shift_Tndiff = new (std::nothrow) int [nkatm_paw]();

   int indx = 0;
   for (auto iia=0; iia<nkatm_paw; ++iia)
   {
      ia = katm_pawtoion[iia];
      shift_Tndiff[iia] = indx;


      for (auto l=0; l<=mult_l[iia]; ++l)
      {
         for (auto m=-l; m<=l; ++m)
         {
            for (auto j=0; j<nprj[ia]; ++j)
            {
               int lj = l_prj[ia][j];
               int mj = m_prj[ia][j];
               int bj = b_prj[ia][j];

               for (auto i=0; i<nprj[ia]; ++i)
               {
                  int li = l_prj[ia][i];
                  int mi = m_prj[ia][i];
                  int bi = b_prj[ia][i];
                  double taunt = util_gaunt(false,l,m,li,mi,lj,mj)*comp_charge_matrix[ia][matindx2(bi,bj,l,nbasis[ia])];
                  if (abs(taunt)>1.0e-15) ++indx;
               }
            }
         }
      }

      nindx_Tndiff[iia] = indx - shift_Tndiff[iia];
   }
   std::cout << "out define taunt" << std::endl;

   lm_Tndiff   = new (std::nothrow) int [indx]();
   iprj_Tndiff = new (std::nothrow) int [indx]();
   jprj_Tndiff = new (std::nothrow) int [indx]();
   coeff_Tndiff = new (std::nothrow) double [indx]();

   indx = 0;
   for (auto iia=0; iia<nkatm_paw; ++iia)
   {
      ia = katm_pawtoion[iia];

      int lm = 0;
      for (auto l=0; l<=mult_l[iia]; ++l)
      {
         for (auto m=-l; m<=l; ++m)
         {
            for (auto j=0; j<nprj[ia]; ++j)
            {
               int lj = l_prj[ia][j];
               int mj = m_prj[ia][j];
               int bj = b_prj[ia][j];

               for (auto i=0; i<nprj[ia]; ++i)
               {
                  int li = l_prj[ia][i];
                  int mi = m_prj[ia][i];
                  int bi = b_prj[ia][i];
                  double taunt = util_gaunt(false,l,m,li,mi,lj,mj)*comp_charge_matrix[ia][matindx2(bi,bj,l,nbasis[ia])];

                  if (abs(taunt)>1.0e-15) 
                  {
                     lm_Tndiff[indx]   = lm;
                     iprj_Tndiff[indx] = i;
                     jprj_Tndiff[indx] = j;
                     coeff_Tndiff[indx] = taunt;
                     ++indx;
                  }
               }
            }
            ++lm;
         }
      }
   }

   //************************************************************
   //*****  define indexing for hartree matrix evalulations  ****
   //************************************************************
   nindx_hartree = new (std::nothrow) int [nkatm_paw]();
   shift_hartree = new (std::nothrow) int [nkatm_paw]();

   indx = 0;
   for (auto iia=0; iia<nkatm_paw; ++iia)
   {
      ia = katm_pawtoion[iia];
      shift_hartree[iia] = indx;

      for (auto j=0; j<nprj[ia]; ++j)
      {
         int lj = l_prj[ia][j];
         int mj = m_prj[ia][j];
         int bj = b_prj[ia][j];
         for (auto i=0; i<nprj[ia]; ++i)
         {
            int li = l_prj[ia][i];
            int mi = m_prj[ia][i];
            int bi = b_prj[ia][i];
            for (auto j1=0; j1<nprj[ia]; ++j1)
            {
               int lj1 = l_prj[ia][j1];
               int mj1 = m_prj[ia][j1];
               int bj1 = b_prj[ia][j1];
               for (auto i1=0; i1<nprj[ia]; ++i1)
               {
                  int li1 = l_prj[ia][i1];
                  int mi1 = m_prj[ia][i1];
                  int bi1 = b_prj[ia][i1];
                  for (auto l=0; l<mult_l[iia]; ++l)
                  {
                     for (auto m=-l; m<=l; ++m)
                     {
                        double taunt = util_gaunt(false,l,m,li,mi,lj,mj)
                                      *util_gaunt(false,l,m,li1,mi1,lj1,mj1)
                                      *hartree_matrix[ia][matindx4(bi,bj,bi1,bj1,l,nbasis[ia])];
                        if (abs(taunt)>1.0e-15) ++indx; 
                     }
                  }

               }
            }

         }
      }
      nindx_hartree[iia] = indx - shift_hartree[iia];
   }

   iprj_hartree = new (std::nothrow) int [indx]();
   jprj_hartree = new (std::nothrow) int [indx]();
   iprj1_hartree = new (std::nothrow) int [indx]();
   jprj1_hartree = new (std::nothrow) int [indx]();
   coeff_hartree = new (std::nothrow) double [indx]();

   indx = 0;
   for (auto iia=0; iia<nkatm_paw; ++iia)
   {
      ia = katm_pawtoion[iia];

      for (auto j=0; j<nprj[ia]; ++j)
      {
         int lj = l_prj[ia][j];
         int mj = m_prj[ia][j];
         int bj = b_prj[ia][j];
         for (auto i=0; i<nprj[ia]; ++i)
         {
            int li = l_prj[ia][i];
            int mi = m_prj[ia][i];
            int bi = b_prj[ia][i];
            for (auto j1=0; j1<nprj[ia]; ++j1)
            {
               int lj1 = l_prj[ia][j1];
               int mj1 = m_prj[ia][j1];
               int bj1 = b_prj[ia][j1];
               for (auto i1=0; i1<nprj[ia]; ++i1)
               {
                  int li1 = l_prj[ia][i1];
                  int mi1 = m_prj[ia][i1];
                  int bi1 = b_prj[ia][i1];
                  for (auto l=0; l<=mult_l[iia]; ++l)
                  {
                     for (auto m=-l; m<=l; ++m)
                     {
                        double taunt = util_gaunt(false,l,m,li,mi,lj,mj)
                                      *util_gaunt(false,l,m,li1,mi1,lj1,mj1)
                                      *hartree_matrix[ia][matindx4(bi,bj,bi1,bj1,l,nbasis[ia])];
                        if (abs(taunt)>1.0e-15)
                        {
                           iprj_hartree[indx]  = i;
                           jprj_hartree[indx]  = j;
                           iprj1_hartree[indx] = i1;
                           jprj1_hartree[indx] = j1;
                           coeff_hartree[indx] = taunt;
                           ++indx;
                        }
                     }
                  }

               }
            }

         }
      }
   }


   //**** initialize the gaussian integrals ****

}





}
