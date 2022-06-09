/* Paw_gintegrals.cpp - 
   Author - Eric Bylaska
*/

#include        <iostream>
#include        <cstdio>
#include        <cstdlib>
#include        <cstring>
#include	<complex>
#include        <cmath>




#include        "util_wgaussian.hpp"
#include        "Paw_gintegrals.hpp"

namespace pwdft {


/*************************************************
 *                                               *
 *          Paw_gintegrals_set_gcount            *
 *                                               *
 *************************************************/

static void Paw_gintegrals_set_gcount(const int nshl3d,
                                      const int taskid, const int np, const int nthr,
                                      const int nion_paw, const int katm_paw[], const int mult_l[],                                    
                                      int& ngauss, int& ngauss_max,
                                      int tgauss[],int tgauss_shift[])
{
   std::memset(tgauss,0,nthr*sizeof(int));

   int nn;
   int pcount = 0;
   int gcount = 0;
   int tcount = 0;
   for (auto iii=0; iii<nion_paw; ++iii)
   {
      int iia = katm_paw[iii];

      //**** calculate on-site integrals ****
      for (auto l1=0; l1<=mult_l[iia]; ++l1)
      for (auto m1=0; m1<=l1; ++m1)
      {
         if (m1==0) nn = 1;
         if (m1>0)  nn = 2;
         if ((pcount%np)==taskid)
         {
            int tid = (gcount%nthr);
            tgauss[tid] += nn;
            tcount += nn;
            gcount += 1;
         }
         ++pcount;

         if (nshl3d>1)
         {
            for (auto l2=0; l2<=mult_l[iia]; ++l2)
            for (auto m2=0; m2<=l2; ++m2)
            {
               if ((m1==0) && (m2==0)) nn = 1;
               if ((m1==0) && (m2>0))  nn = 2;
               if ((m1>0)  && (m2==0)) nn = 2;
               if ((m1>0)  && (m2>0))  nn = 4;
               if ((pcount%np)==taskid) 
               {
                  int tid = (gcount%nthr);
                  tgauss[tid] += nn;
                  tcount += nn;
                  gcount += 1;
               }
               ++pcount;
            }
         }
      }
      
      //**** calculate IJ integrals ****
      for (auto jjj=iii+1; jjj<nion_paw; ++jjj)
      {
         int jja = katm_paw[jjj];

         for (auto l1=0; l1<=mult_l[iia]; ++l1)
         for (auto m1=0; m1<=l1; ++m1)
            for (auto l2=0; l2<=mult_l[jja]; ++l2)
            for (auto m2=0; m2<=l2; ++m2)
            {
            
               if ((m1==0) && (m2==0)) nn = 1;
               if ((m1==0) && (m2>0))  nn = 2;
               if ((m1>0)  && (m2==0)) nn = 2;
               if ((m1>0)  && (m2>0))  nn = 4;
               if ((pcount%np)==taskid)
               {
                  int tid = (gcount%nthr);
                  tgauss[tid] += nn;
                  tcount += nn;
                  gcount += 1;
               }
               ++pcount;
            } 
      }
   }
   ngauss_max = tcount;
   ngauss     = tcount;

   tcount = 0;
   for (auto l1=0; l1<nthr; ++l1)
   {
      tgauss_shift[l1] = tcount;
      tcount += tgauss[l1];
   }
}


/*************************************************
 *                                               *
 *             Paw_gintegrals_stripper           *
 *                                               *
 *************************************************/
/*
   This routine is used to remove unecessary integrals
*/
static void Paw_gintegral_stripper(const int ng_in,
                                   const int iii1_in[], const int iii2_in[],
                                   const int lm1_in[],  const int lm2_in[],
                                   const double e_in[], const double f_in[],
                                   int& ng_out,
                                   int iii1_out[], int iii2_out[],
                                   int lm1_out[],  int lm2_out[],
                                   double e_out[], double f_out[])
{
   double  tole=1.0e-25;

   ng_out = 0;
   for (auto i=0; i<ng_in; ++i)
   {
      if (std::abs(e_in[i])>tole) 
      {
         iii1_out[ng_out] = iii1_in[i];
         iii2_out[ng_out] = iii2_in[i];
         lm1_out[ng_out]  = lm1_in[i];
         lm2_out[ng_out]  = lm2_in[i];
         e_out[ng_out]    = e_in[i];
         f_out[3*ng_out]   = f_in[3*i];
         f_out[3*ng_out+1] = f_in[3*i+1];
         f_out[3*ng_out+2] = f_in[3*i+2];

         ++ng_out;
      }
   }
}

/*******************************************************
 *                                                     *
 *              Paw_WGaussian_block                    *
 *                                                     *
 *******************************************************/

#define	m1pow(n)	((double) (1 - (((n)&1) << 1)))

static void Paw_WGaussian_block(const int l1, const int mod_m1, const double sigma1,
                                const int l2, const int mod_m2, const double sigma2,
                                const double  sigma_smooth, const double Rab[],
                                int& n, int lm1[], int lm2[], double W[])
{
   std::complex<double> CW4,CW4p,CW4m,CW4pp,CW4pm,CW4mp,CW4mm;

   if ((mod_m1==0) && (mod_m2==0))
   {
      n = 1;
      lm1[0] = l1*(l1+1);
      lm2[0] = l2*(l2+1);
 
      CW4 = util_CWGaussian3(l1,mod_m1,sigma1,l2,mod_m2,sigma2,sigma_smooth,Rab);
      W[0] = CW4.real();
   }
   else if (mod_m1==0)
   {
      n = 2;
      lm1[0] = l1*(l1+1);
      lm2[0] = l2*(l2+1) - mod_m2;

      lm1[1] = l1*(l1+1);
      lm2[1] = l2*(l2+1) + mod_m2;

      CW4p = util_CWGaussian3(l1,mod_m1,sigma1,l2,mod_m2,sigma2,sigma_smooth,Rab);
      CW4m = util_CWGaussian3(l1,mod_m1,sigma1,l2,-mod_m2,sigma2,sigma_smooth,Rab);

      //** m2<0 **
      CW4 = (CW4m - m1pow(mod_m2) * CW4p)*std::complex<double>(0.0,1.0/sqrt(2.0));
      W[0] = CW4.real();

      //** m2>0 **
      CW4 = (CW4m + m1pow(mod_m2) * CW4p)/sqrt(2.0);
      W[1] = CW4.real();
   }
   else if (mod_m2==0)
   {
      n = 2;
      lm1[0] = l1*(l1+1) - mod_m1;
      lm2[0] = l2*(l2+1);

      lm1[1] = l1*(l1+1) + mod_m1;
      lm2[1] = l2*(l2+1);

      CW4p = util_CWGaussian3(l1,mod_m1,sigma1,l2,mod_m2,sigma2,sigma_smooth,Rab);
      CW4m = util_CWGaussian3(l1,-mod_m1,sigma1,l2,mod_m2,sigma2,sigma_smooth,Rab);

      //** m1<0 **
      CW4 = (CW4m - m1pow(mod_m1) * CW4p) * std::complex<double>(0.0,1.0/sqrt(2.0));
      W[0] = CW4.real();

      //** m1>0 **
      CW4 = (CW4m + m1pow(mod_m1) * CW4p)/sqrt(2.0);
      W[1] = CW4.real();
   }
   else
   {
      n = 4;
      lm1[0] = l1*(l1+1) - mod_m1;
      lm2[0] = l2*(l2+1) - mod_m2;

      lm1[1] = l1*(l1+1) - mod_m1;
      lm2[1] = l2*(l2+1) + mod_m2;

      lm1[2] = l1*(l1+1) + mod_m1;
      lm2[2] = l2*(l2+1) - mod_m2;

      lm1[3] = l1*(l1+1) + mod_m1;
      lm2[3] = l2*(l2+1) + mod_m2;

      CW4pp = util_CWGaussian3(l1, mod_m1,sigma1,l2, mod_m2,sigma2,sigma_smooth,Rab);
      CW4pm = util_CWGaussian3(l1, mod_m1,sigma1,l2,-mod_m2,sigma2,sigma_smooth,Rab);
      CW4mp = util_CWGaussian3(l1,-mod_m1,sigma1,l2, mod_m2,sigma2,sigma_smooth,Rab);
      CW4mm = util_CWGaussian3(l1,-mod_m1,sigma1,l2,-mod_m2,sigma2,sigma_smooth,Rab);

      //** m1<0 and m2<0 **
      CW4 = -(CW4mm
             + m1pow((mod_m1+mod_m2)) * CW4pp
             - m1pow(mod_m1)     * CW4pm
             - m1pow(mod_m2)     * CW4mp)/2.0;
      W[0] = CW4.real();

      //** m1<0 and m2>0 **
      CW4 = (CW4mm
            - m1pow((mod_m1+mod_m2)) * CW4pp
            - m1pow(mod_m1)          * CW4pm
            + m1pow(mod_m2)          * CW4mp)
            * std::complex<double>(0.0,1.0/2.0);
      W[1] = CW4.real();

      //** m1>0 and m2<0 **
      CW4 = (CW4mm
            - m1pow((mod_m1+mod_m2)) * CW4pp
            + m1pow(mod_m1)          * CW4pm
            - m1pow(mod_m2)          * CW4mp)
           * std::complex<double>(0.0,1.0/2.0);
      W[2] = CW4.real();

      //** m1>0 and m2>0 **
      CW4 =  (CW4mm
             + m1pow((mod_m1+mod_m2)) * CW4pp
             + m1pow(mod_m1)          * CW4pm
             + m1pow(mod_m2)          * CW4mp)/2.0;
      W[3] = CW4.real();
   }
}


/*******************************************************
 *                                                     *
 *              Paw_WGaussian2_block                   *
 *                                                     *
 *******************************************************/

static void Paw_WGaussian2_block(const int l1, const int mod_m1, const double sigma1, 
                                 const int l2, const int mod_m2, 
                                 const double sigma_smooth, const double Rab[], 
                                 int& n, int lm1[], int lm2[], double W[])
{
   std::complex<double> CW4,CW4p,CW4m,CW4pp,CW4pm,CW4mp,CW4mm;

   if ((mod_m1==0) && (mod_m2==0))
   {
      n = 1;
      lm1[0] = l1*(l1+1);
      lm2[0] = l2*(l2+1);
 
      CW4 = util_CWGaussian2(l1,mod_m1,sigma1,l2,mod_m2,sigma_smooth,Rab);
      W[0] = CW4.real();
   }
   else if (mod_m1==0)
   {
      n = 2;
      lm1[0] = l1*(l1+1);
      lm2[0] = l2*(l2+1) - mod_m2;

      lm1[1] = l1*(l1+1);
      lm2[1] = l2*(l2+1) + mod_m2;

      CW4p = util_CWGaussian2(l1,mod_m1,sigma1,l2, mod_m2,sigma_smooth,Rab);
      CW4m = util_CWGaussian2(l1,mod_m1,sigma1,l2,-mod_m2,sigma_smooth,Rab);

      // ** m2<0 **
      CW4 = (CW4m - m1pow(mod_m2) * CW4p)*std::complex<double>(0.0,1.0/sqrt(2.0));
      W[0] = CW4.real();

      // ** m2>0 **
      CW4 = (CW4m + m1pow(mod_m2) * CW4p)/sqrt(2.0);
      W[1] = CW4.real();
   }
   else if (mod_m2==0)
   {
      n = 2;
      lm1[0] = l1*(l1+1) - mod_m1;
      lm2[0] = l2*(l2+1);

      lm1[1] = l1*(l1+1) + mod_m1;
      lm2[1] = l2*(l2+1);

      CW4p = util_CWGaussian2(l1, mod_m1,sigma1,l2,mod_m2,sigma_smooth,Rab);
      CW4m = util_CWGaussian2(l1,-mod_m1,sigma1,l2,mod_m2,sigma_smooth,Rab);

      // ** m1<0 **
      CW4 = (CW4m - m1pow(mod_m1) * CW4p) * std::complex<double>(0.0,1.0/sqrt(2.0));
      W[0] = CW4.real();

      // ** m1>0 **
      CW4 = (CW4m + m1pow(mod_m1) * CW4p)/sqrt(2.0);
      W[1] = CW4.real();
   }
   else
   {
      n = 4;
      lm1[0] = l1*(l1+1) - mod_m1;
      lm2[0] = l2*(l2+1) - mod_m2;

      lm1[1] = l1*(l1+1) - mod_m1;
      lm2[1] = l2*(l2+1) + mod_m2;

      lm1[2] = l1*(l1+1) + mod_m1;
      lm2[2] = l2*(l2+1) - mod_m2;

      lm1[3] = l1*(l1+1) + mod_m1;
      lm2[3] = l2*(l2+1) + mod_m2;

      CW4pp = util_CWGaussian2(l1, mod_m1,sigma1,l2, mod_m2,sigma_smooth,Rab);
      CW4pm = util_CWGaussian2(l1, mod_m1,sigma1,l2,-mod_m2,sigma_smooth,Rab);
      CW4mp = util_CWGaussian2(l1,-mod_m1,sigma1,l2, mod_m2,sigma_smooth,Rab);
      CW4mm = util_CWGaussian2(l1,-mod_m1,sigma1,l2,-mod_m2,sigma_smooth,Rab);

      //** m1<0 and m2<0 **
      CW4 = -(CW4mm
             + m1pow((mod_m1+mod_m2)) * CW4pp
             - m1pow(mod_m1)          * CW4pm
             - m1pow(mod_m2)          * CW4mp)/2.0;
      W[0] = CW4.real();

      //** m1<0 and m2>0 **
      CW4 = (CW4mm
            - m1pow((mod_m1+mod_m2)) * CW4pp
            - m1pow(mod_m1)          * CW4pm
            + m1pow(mod_m2)          * CW4mp)
            *std::complex<double>(0.0,1.0/2.0);
      W[1] = CW4.real();

      //** m1>0 and m2<0 **
      CW4 = (CW4mm
            - m1pow((mod_m1+mod_m2)) * CW4pp
            + m1pow(mod_m1)          * CW4pm
            - m1pow(mod_m2)          * CW4mp)
           *std::complex<double>(0.0,1.0/2.0);
      W[2] = CW4.real();

      //** m1>0 and m2>0 **
      CW4 =  (CW4mm
             + m1pow((mod_m1+mod_m2)) * CW4pp
             + m1pow(mod_m1)          * CW4pm
             + m1pow(mod_m2)          * CW4mp)/2.0;
      W[3] = CW4.real();
   }
}



/*******************************************************
 *                                                     *
 *              Paw_dWGaussian_block                   *
 *                                                     *
 *******************************************************/

static void Paw_dWGaussian_block(const int l1, const int mod_m1, const double sigma1,
                                 const int l2, const int mod_m2, const double sigma2,
                                 const double sigma_smooth, const double Rab[],
                                 int& n, int lm1[], int lm2[], double W[], double dW[])
{
   std::complex<double> CW4,CW4p,CW4m,CW4pp,CW4pm,CW4mp,CW4mm;
   std::complex<double> dCW4[3],dCW4p[3],dCW4m[3];
   std::complex<double> dCW4pp[3],dCW4pm[3];
   std::complex<double> dCW4mp[3],dCW4mm[3];

   if ((mod_m1==0) && (mod_m2==0)) 
   {
      n = 1;
      lm1[0] = l1*(l1+1);
      lm2[0] = l2*(l2+1);
 
      util_dCWGaussian3(l1,mod_m1,sigma1,l2,mod_m2,sigma2,sigma_smooth,Rab,CW4,dCW4);

      W[0]  = CW4.real();
      dW[0] = dCW4[0].real();
      dW[1] = dCW4[1].real();
      dW[2] = dCW4[2].real();
   }
   else if (mod_m1==0)
   {
      n = 2;
      lm1[0] = l1*(l1+1);
      lm2[0] = l2*(l2+1) - mod_m2;

      lm1[1] = l1*(l1+1);
      lm2[1] = l2*(l2+1) + mod_m2;

      util_dCWGaussian3(l1,mod_m1,sigma1,l2, mod_m2,sigma2,sigma_smooth,Rab,CW4p,dCW4p);
      util_dCWGaussian3(l1,mod_m1,sigma1,l2,-mod_m2,sigma2,sigma_smooth,Rab,CW4m,dCW4m);

      //** m2<0 **
      CW4 = (CW4m - m1pow(mod_m2) * CW4p) * std::complex<double>(0.0,1.0/sqrt(2.0));
      W[0] = CW4.real();
      for (auto i=0; i<3; ++i)
      {
         CW4 = (dCW4m[i] - m1pow(mod_m2) * dCW4p[i]) * std::complex<double>(0.0,1.0/sqrt(2.0));
         dW[i] = CW4.real();
      }

      //** m2>0 **
      CW4 = (CW4m + m1pow(mod_m2) * CW4p)/sqrt(2.0);
      W[1] = CW4.real();
      for (auto i=0; i<3; ++i)
      {
         CW4 = (dCW4m[i] + m1pow(mod_m2) * dCW4p[i])/sqrt(2.0);
         dW[i+3] = CW4.real();
      }
   }
   else if (mod_m2==0) 
   {
      n = 2;
      lm1[0] = l1*(l1+1) - mod_m1;
      lm2[0] = l2*(l2+1);

      lm1[1] = l1*(l1+1) + mod_m1;
      lm2[1] = l2*(l2+1);

      util_dCWGaussian3(l1, mod_m1,sigma1,l2,mod_m2,sigma2,sigma_smooth,Rab,CW4p,dCW4p);
      util_dCWGaussian3(l1,-mod_m1,sigma1,l2,mod_m2,sigma2,sigma_smooth,Rab,CW4m,dCW4m);

      //** m1<0 **
      CW4 = (CW4m - m1pow(mod_m1) * CW4p) * std::complex<double>(0.0,1.0/sqrt(2.0));
      W[0] = CW4.real();
      for (auto i=0; i<3; ++i)
      {
         CW4 = (dCW4m[i] - m1pow(mod_m1) * dCW4p[i]) *std::complex<double>(0.0,1.0/sqrt(2.0));
         dW[i] = CW4.real();
      }

      //** m1>0 **
      CW4 = (CW4m + m1pow(mod_m1) * CW4p)/sqrt(2.0);
      W[1] = CW4.real();
      for (auto i=0; i<3; ++i)
      {
         CW4 = (dCW4m[i] + m1pow(mod_m1) * dCW4p[i])/sqrt(2.0);
         dW[i+3] = CW4.real();
      }

   }
   else
   {
      n = 4;
      lm1[0] = l1*(l1+1) - mod_m1;
      lm2[0] = l2*(l2+1) - mod_m2;

      lm1[1] = l1*(l1+1) - mod_m1;
      lm2[1] = l2*(l2+1) + mod_m2;

      lm1[2] = l1*(l1+1) + mod_m1;
      lm2[2] = l2*(l2+1) - mod_m2;

      lm1[3] = l1*(l1+1) + mod_m1;
      lm2[3] = l2*(l2+1) + mod_m2;

      util_dCWGaussian3(l1, mod_m1,sigma1,l2, mod_m2,sigma2,sigma_smooth,Rab,CW4pp,dCW4pp);
      util_dCWGaussian3(l1, mod_m1,sigma1,l2,-mod_m2,sigma2,sigma_smooth,Rab,CW4pm,dCW4pm);
      util_dCWGaussian3(l1,-mod_m1,sigma1,l2, mod_m2,sigma2,sigma_smooth,Rab,CW4mp,dCW4mp);
      util_dCWGaussian3(l1,-mod_m1,sigma1,l2,-mod_m2,sigma2,sigma_smooth,Rab,CW4mm,dCW4mm);

      //** m1<0 and m2<0 **
      CW4 = -(CW4mm
             + m1pow((mod_m1+mod_m2)) * CW4pp
             - m1pow(mod_m1)          * CW4pm
             - m1pow(mod_m2)          * CW4mp)/2.0;
      W[0] = CW4.real();
      for (auto i=0; i<3; ++i)
      {
         CW4 = -(dCW4mm[i]
             + m1pow((mod_m1+mod_m2)) * dCW4pp[i]
             - m1pow(mod_m1)          * dCW4pm[i]
             - m1pow(mod_m2)          * dCW4mp[i])/2.0;
         dW[i] = CW4.real();
      }

      //** m1<0 and m2>0 **
      CW4 = (CW4mm
            - m1pow((mod_m1+mod_m2)) * CW4pp
            - m1pow(mod_m1)          * CW4pm
            + m1pow(mod_m2)          * CW4mp)*std::complex<double>(0.0,1.0/2.0);
      W[1] = CW4.real();
      for (auto i=0; i<3; ++i)
      {
         CW4 = (dCW4mm[i]
            - m1pow((mod_m1+mod_m2)) * dCW4pp[i]
            - m1pow(mod_m1)          * dCW4pm[i]
            + m1pow(mod_m2)          * dCW4mp[i])*std::complex<double>(0.0,1.0/2.0);
         dW[i+3] = CW4.real();
      }

      //** m1>0 and m2<0 **
      CW4 = (CW4mm
            - m1pow((mod_m1+mod_m2)) * CW4pp
            + m1pow(mod_m1)          * CW4pm
            - m1pow(mod_m2)          * CW4mp)*std::complex<double>(0.0,1.0/2.0);
      W[2] = CW4.real();
      for (auto i=0; i<3; ++i)
      {
         CW4 = (dCW4mm[i]
            - m1pow((mod_m1+mod_m2)) * dCW4pp[i]
            + m1pow(mod_m1)          * dCW4pm[i]
            - m1pow(mod_m2)          * dCW4mp[i])*std::complex<double>(0.0,1.0/2.0);
         dW[i+6] = CW4.real();
      }

      //** m1>0 and m2>0 **
      CW4 =  (CW4mm
             + m1pow((mod_m1+mod_m2)) * CW4pp
             + m1pow(mod_m1)          * CW4pm
             + m1pow(mod_m2)          * CW4mp)/2.0;
      W[3] = CW4.real();
      for (auto i=0; i<3; ++i)
      {
         CW4 =  (dCW4mm[i]
             + m1pow((mod_m1+mod_m2)) * dCW4pp[i]
             + m1pow(mod_m1)          * dCW4pm[i]
             + m1pow(mod_m2)          * dCW4mp[i])/2.0;
         dW[i+9] = CW4.real();
      }
   }
}


/*******************************************************
 *                                                     *
 *            Paw_gintegrals::Paw_gintegrals           *
 *                                                     *
 *******************************************************/

Paw_gintegrals::Paw_gintegrals(Parallel *myparall0, Ion *myion0, Ewald *myewald0, const bool periodic0, const int lm_size_max0,
                               const int nion_paw0, int katm_paw0[], int mult_l0[], int ion_pawtoion0[], double sigma_paw0[], const double sigma_smooth0)
{
   myparall = myparall0;
   myion    = myion0;
   myewald  = myewald0;

   periodic = periodic0;
   lm_size_max = lm_size_max0;
   nion_paw = nion_paw0;
   katm_paw = katm_paw0;
   mult_l   = mult_l0;
   ion_pawtoion = ion_pawtoion0;
   sigma_paw = sigma_paw0;
   sigma_smooth = sigma_smooth0;

   int nshl3d=1;
   if (periodic) nshl3d = myewald->nshl3d();

   tgauss       = new int[myparall->maxthreads()];
   tgauss_shift = new int[myparall->maxthreads()];

   Paw_gintegrals_set_gcount(nshl3d,myparall->taskid(),myparall->np(),myparall->maxthreads(),
                             nion_paw,katm_paw,mult_l,
                             ngauss,ngauss_max,
                             tgauss,tgauss_shift);

   lm1_gauss = new int[ngauss_max];
   lm2_gauss = new int[ngauss_max];
   iii1_gauss = new int[ngauss_max];
   iii2_gauss = new int[ngauss_max];
   e_gauss = new double[ngauss_max];
   f_gauss = new double[3*ngauss_max];  
}


/*******************************************************
 *                                                     *
 *               Paw_gintegrals::set                   *
 *                                                     *
 *******************************************************/
/*
   The logic of this routine needs to be completely reworked for threading.
   It's well designed for MPI parallelism, so one option is to expand all the
   data structures over tasks and threads instead of just tasks.
   Another option is to define thread shifts for indx.... However the threshold
   check with tole would have to be eliminated.
*/
void Paw_gintegrals::set(const bool move)
{
   double tole = 1.0e-25;

   int taskid = myparall->taskid();
   int np     = myparall->np();
   int tid    = myparall->threadid();
   int nthr   = myparall->nthreads();
   int shft   = tgauss_shift[tid];

   int    lm1[4],lm2[4];
   double R1[3],R12[3],Rab[3],Rba[3];
   double W1,W2,W3,W4,dW1[3],dW2[3],dW3(3),dW4[3];
   double W[4],dW[12];
   double e1[4],de1[12];

   int lm1_tauss[ngauss_max],  lm2_tauss[ngauss_max];
   int iii1_tauss[ngauss_max], iii2_tauss[ngauss_max];
   double e_tauss[ngauss_max], f_tauss[3*ngauss_max];

   std::memset(e_tauss,0,ngauss_max*sizeof(double));
   std::memset(f_tauss,0,3*ngauss_max*sizeof(double));

   int nshl3d=1;
   double *rcell;

   if (periodic) 
   {
      rcell  = myewald->rcell;
      nshl3d = myewald->nshl3d();
   }
   else
   {
      nshl3d=1;
      rcell = new double[3];
      rcell[0] = 0.0;
      rcell[1] = 0.0;
      rcell[2] = 0.0;
   }

   int nn;
   int pcount = 0;
   int gcount = 0;
   int indx   = 0;
   for (auto iii=0; iii<nion_paw; ++iii)
   {
      int iia = katm_paw[iii];
      double s1 = sigma_paw[iia];

      //**** calculate on-site integrals ****
      for (auto l1=0; l1<=mult_l[iia]; ++l1)
      for (auto m1=0; m1<=l1; ++m1)
      { 
         if (m1==0) nn = 1;
         if (m1==0) nn = 2;
         if ((pcount%np)==taskid)
         {
            if ((gcount%nthr)==tid)
            {
               W1 = util_UGaussian(l1,m1,s1,l1,m1,s1);
               W2 = util_UGaussian(l1,m1,s1,l1,m1,sigma_smooth);
               W4 = util_UGaussian(l1,m1,sigma_smooth,l1,m1,sigma_smooth);
               e1[0] = 0.5*W1 + 0.5*W4 - W2;
               lm1[0] = l1*(l1+1) + m1;
               if (nn>1)
               {
                  W1 = util_UGaussian(l1,-m1,s1,l1,-m1,s1);
                  W2 = util_UGaussian(l1,-m1,s1,l1,-m1,sigma_smooth);
                  W4 = util_UGaussian(l1,-m1,sigma_smooth,l1,-m1,sigma_smooth);
                  e1[1] = 0.5*W1 + 0.5*W4 - W2;
                  lm1[1] = l1*(l1+1) - m1;
               }

               for (auto i=0; i<nn; ++i)
               {
                  int inds = indx + shft;
                  e_tauss[inds]    = e1[i];
                  lm1_tauss[inds]  = (iii)*2*lm_size_max+lm1[i];
                  lm2_tauss[inds]  = (iii)*2*lm_size_max+lm1[i];
                  iii1_tauss[inds] = iii;
                  iii2_tauss[inds] = iii;
                  ++indx;
               }
            }
            ++gcount;
         }
         ++pcount;

         if (nshl3d>1)
         {
            for (auto l2=0; l2<=mult_l[iia]; ++l2)
            for (auto m2=0; m2<=l2; ++m2)
            {
               if ((m1==0) && (m2==0)) nn = 1;
               if ((m1==0) && (m2>0))  nn = 2;
               if ((m1>0)  && (m2==0)) nn = 2;
               if ((m1>0)  && (m2>0))  nn = 4;
               if ((pcount%np)==taskid)
               {
                  if ((gcount%nthr)==tid)
                  {
                     std::memset(e1,0,4*sizeof(double));
                     for (auto l=1; l<nshl3d; ++l)
                     {
                        Rab[0] = rcell[l];
                        Rab[1] = rcell[l+nshl3d];
                        Rab[2] = rcell[l+2*nshl3d];
                        double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
                        if (R<(4*sigma_smooth)) 
                        {
                           int n;
                           Paw_WGaussian2_block(l1,m1,s1,l2,m2,sigma_smooth,Rab,n,lm1,lm2,W);
                           for (auto i=0; i<n; ++i)
                              e1[i] = e1[i] + 0.5*W[i];
                        }
                     }

                     for (auto i=0; i<nn; ++i)
                     {
                        int inds = indx + shft;
                        e_tauss[inds]    = e1[i];
                        lm1_tauss[inds]  = (iii)*2*lm_size_max+lm1[i];
                        lm2_tauss[inds]  = (iii)*2*lm_size_max+lm2[i];
                        iii1_tauss[inds] = iii;
                        iii2_tauss[inds] = iii;

                        ++indx;
                     }
                  }
                  ++gcount;
               }
               ++pcount;
            }
         }
      }
      

      //**** calculate IJ integrals ****
      int ii = ion_pawtoion[iii];
      R1[0] = myion->rion(0,ii);
      R1[1] = myion->rion(1,ii);
      R1[2] = myion->rion(2,ii);
      for (auto jjj=iii+1; jjj<nion_paw; ++jjj) 
      {
         int jja = katm_paw[jjj];
         double s2 = sigma_paw[jja];

         int jj = ion_pawtoion[jjj];
         R12[0] = R1[0] - myion->rion(0,jj);
         R12[1] = R1[1] - myion->rion(1,jj);
         R12[2] = R1[2] - myion->rion(2,jj);

         for (auto l1=0; l1<=mult_l[iia]; ++l1)
         for (auto m1=0; m1<=l1; ++m1)
         {
            for (auto l2=0; l2<=mult_l[jja]; ++l2)
            for (auto m2=0; m2<=l2; ++m2)
            {
               if ((m1==0) && (m2==0)) nn = 1;
               if ((m1==0) && (m2>0))  nn = 2;
               if ((m1>0)  && (m2==0)) nn = 2;
               if ((m1>0)  && (m2>0))  nn = 4;
               if ((pcount%np)==taskid)
               {
                  if ((gcount%nthr)==tid)
                  {
                     std::memset(e1, 0,  4*sizeof(double));
                     std::memset(de1,0,3*4*sizeof(double));
                     
                     for (auto l=0; l<nshl3d; ++l)
                     {
                        Rab[0] = R12[0] + rcell[l];
                        Rab[1] = R12[1] + rcell[l+nshl3d];
                        Rab[2] = R12[2] + rcell[l+2*nshl3d];
                        double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
                        if (R<(4*sigma_smooth))
                        {
                           if (move) 
                           {
                              int n;
                              Paw_dWGaussian_block(l1,m1,s1,l2,m2,s2,sigma_smooth,Rab,n,lm1,lm2,W,dW);
                              for (auto i=0; i<n; ++i)
                              {
                                 e1[i]  +=  W[i];
                                 de1[0+3*i] +=  dW[0+3*i];
                                 de1[1+3*i] +=  dW[1+3*i];
                                 de1[2+3*i] +=  dW[2+3*i];
                              }
                           }
                           else
                           {
                              int n;
                              Paw_WGaussian_block(l1,m1,s1,l2,m2,s2,sigma_smooth,Rab,n,lm1,lm2,W);
                              for (auto i=0; i<n; ++i)
                              {
                                 e1[i] +=  W[i];
                              }
                           }
                        }
                     }
                  
                     for (auto i=0; i<nn; ++i)
                     {
                        int inds = indx + shft;
                        e_tauss[inds] = e1[i];
                        if (move)
                        {
                           f_tauss[3*inds]   = de1[0+3*i];
                           f_tauss[3*inds+1] = de1[1+3*i];
                           f_tauss[3*inds+2] = de1[2+3*i];
                        }
                        lm1_tauss[inds] = (iii)*2*lm_size_max+lm1[i];
                        lm2_tauss[inds] = (jjj)*2*lm_size_max+lm2[i];
                        iii1_tauss[inds] = iii;
                        iii2_tauss[inds] = jjj;

                        ++indx;
                     }
                  }
                  ++gcount;
               }
               ++pcount;
            }

         }

      }

   }

   Paw_gintegral_stripper(ngauss_max,iii1_tauss,iii2_tauss,lm1_tauss,lm2_tauss,e_tauss,f_tauss,
                          ngauss,    iii1_gauss,iii2_gauss,lm1_gauss,lm2_gauss,e_gauss,f_gauss);


   //Need to have barrier before deallocation below.  This barrier could be removed
   //if the tauss variables were on the heap and not deallocated rather than stack.


   //**** deallocate rcell memory ****
   if (!periodic) delete [] rcell;
}




}
