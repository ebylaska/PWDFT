/* Psp1d_pawppv1.cpp - 
   Author - Eric Bylaska
*/

#include	<iostream>
#include	<cstdio>
#include	<cstdlib>
#include	<cstring>
#include	<cmath>
//#include        "blas.h"

using namespace std;


#include	"Parallel.hpp"
#include	"util.hpp"
#include	"Psp1d_pawppv1.hpp"

namespace pwdft {
using namespace pwdft;


#define FMT1    "%lf"
#define FMT2    " %lf %lf"
#define FMT10   "%10.3lf %10.3lf %10.3lf"



/*************************************************
 *                                               *
 *        pawppv1_derivative_ngrid               *
 *                                               *
 *************************************************/

/* This routine computes the seven point derivative of f.
   where f and df are stored on a logarithmic grid. The
   dimensions of f and df are, f(1:ng), and df(1:ng)     */

static void pawppv1_derivative_ngrid(const int ng, const double log_amesh,
                                      const double r[], const double f[], double df[])
{
   double one_over_60 = 1.0/60.0;
   double aa = one_over_60/log_amesh;

   df[0] = aa*(-147.00*f[0]
              + 360.00*f[1]
              - 450.00*f[2]
              + 400.00*f[3]
              - 225.00*f[4]
              +  72.00*f[5]
              -  10.00*f[6])/r[0];
   df[1] = aa*( -10.00*f[0]
              -  77.00*f[1]
              + 150.00*f[2]
              - 100.00*f[3]
              +  50.00*f[4]
              -  15.00*f[5]
              +   2.00*f[6])/r[1];
   df[2] = aa*(  +2.00*f[0]
              -  24.00*f[1]
              -  35.00*f[2]
              +  80.00*f[3]
              -  30.00*f[4]
              +   8.00*f[5]
              -   1.00*f[6])/r[2];
   df[ng-3] = aa*(  -2.00*f[ng-1]
                 +  24.00*f[ng-2]
                 +  35.00*f[ng-3]
                 -  80.00*f[ng-4]
                 +  30.00*f[ng-5]
                 -   8.00*f[ng-6]
                 +   1.00*f[ng-7])/r[ng-3];
   df[ng-2] = aa*( +10.00*f[ng-1]
                 +  77.00*f[ng-2]
                 - 150.00*f[ng-3]
                 + 100.00*f[ng-4]
                 -  50.00*f[ng-5]
                 +  15.00*f[ng-6]
                 -   2.00*f[ng-7])/r[ng-2];
   df[ng-1] = aa*(+147.00*f[ng-1]
                 - 360.00*f[ng-2]
                 + 450.00*f[ng-3]
                 - 400.00*f[ng-4]
                 + 225.00*f[ng-5]
                 -  72.00*f[ng-6]
                 +  10.00*f[ng-7])/r[ng-1];
   for (auto i=3; i<(ng-3); ++i)
      df[i] = aa*(  -1.00*f[i-3]
                 +   9.00*f[i-2]
                 -  45.00*f[i-1]
                 +  45.00*f[i+1]
                 -   9.00*f[i+2]
                 +   1.00*f[i+3])/r[i];
}





/* Constructors */

/*******************************************
 *                                         *
 *     Psp1d_pawppv1::Psp1d_pawppv1        *
 *                                         *
 *******************************************/
Psp1d_pawppv1::Psp1d_pawppv1(Parallel *myparall, const char *psp_name)
{
   double xx;
   FILE *fp;
   int nn;


   if (myparall->is_master())
   {
      fp = std::fopen(psp_name,"r");

      std::fscanf(fp,"%d",&psp_type);
      std::fscanf(fp,"%s",atom);
      std::fscanf(fp,FMT1,&zv);
      std::fscanf(fp,FMT1,&r1);
      std::fscanf(fp,FMT1,&rmax);
      std::fscanf(fp,"%d",&n1dgrid);
      std::fscanf(fp,"%d",&nbasis);
      for (int l=0; l<nbasis; ++l)
         std::fscanf(fp,FMT1,&rc[l]);
      std::fscanf(fp,"%d",&icut);
      std::fscanf(fp,"%s",comment);
      std::fscanf(fp,FMT1,&core_kin_energy);
   }

   myparall->Brdcst_iValue(0,0,&psp_type);
   myparall->Brdcst_Values(0,0,1,&zv);
   myparall->Brdcst_Values(0,0,1,&r1);
   myparall->Brdcst_Values(0,0,1,&rmax);
   myparall->Brdcst_iValue(0,0,&n1dgrid);
   myparall->Brdcst_iValue(0,0,&nbasis);
   myparall->Brdcst_Values(0,0,nbasis,rc);
   myparall->Brdcst_iValue(0,0,&icut);

   myparall->Brdcst_cValues(0,0,80,comment);
   myparall->Brdcst_cValues(0,0,2,atom);

   /* define rgrid */
   log_amesh = log(rmax/r1)/((double) n1dgrid-1);
   amesh     = exp(log_amesh);
   rgrid = (double *) new double[n1dgrid];
   rgrid[0] = r1;
   for (auto i=1; i<n1dgrid; ++i)
      rgrid[i] = rgrid[i-1]*amesh;

   /* allocate rest of grid data */
   nae = (int *) new int[nbasis];
   nps = (int *) new int[nbasis];
   lps = (int *) new int[nbasis];
   eig           = (double *) new double[nbasis];
   phi_ae        = (double *) new double[nbasis*n1dgrid];
   dphi_ae       = (double *) new double[nbasis*n1dgrid];
   phi_ps        = (double *) new double[nbasis*n1dgrid];
   dphi_ps       = (double *) new double[nbasis*n1dgrid];
   core_ae       = (double *) new double[n1dgrid];
   core_ps       = (double *) new double[n1dgrid];
   core_ae_prime = (double *) new double[n1dgrid];
   core_ps_prime = (double *) new double[n1dgrid];
   v_ps          = (double *) new double[n1dgrid];
   prj_ps        = (double *) new double[nbasis*n1dgrid];
   prj_ps0       = (double *) new double[nbasis*n1dgrid];
   

   if (myparall->is_master())
   {
      for (auto j=0; j<nbasis; ++j)
      {
         std::fscanf(fp,"%d",&nae[j]);
         std::fscanf(fp,FMT1,&eig[j]);
         std::fscanf(fp,"%d",&nps[j]);
         std::fscanf(fp,"%d",&lps[j]);
      }

      for (auto j=0; j<nbasis;  ++j)
      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&phi_ae[i+j*n1dgrid]);

      for (auto j=0; j<nbasis;  ++j)
      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&dphi_ae[i+j*n1dgrid]);

      for (auto j=0; j<nbasis;  ++j)
      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&phi_ps[i+j*n1dgrid]);

      for (auto j=0; j<nbasis;  ++j)
      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&dphi_ps[i+j*n1dgrid]);

      for (auto j=0; j<nbasis;  ++j)
      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&prj_ps[i+j*n1dgrid]);

      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&core_ae[i]);

      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&core_ps[i]);

      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&v_ps[i]);

      std::fscanf(fp,FMT1,&sigma);
      std::fscanf(fp,FMT1,&zion);

      for (auto j=0; j<nbasis;  ++j)
      for (auto i=0; i<n1dgrid; ++i)
         std::fscanf(fp,FMT1,&prj_ps0[i+j*n1dgrid]);

      std::fclose(fp);
   }

   myparall->Brdcst_iValues(0,0,nbasis,nae);
   myparall->Brdcst_iValues(0,0,nbasis,nps);
   myparall->Brdcst_iValues(0,0,nbasis,lps);
   myparall->Brdcst_Values(0,0,nbasis,eig);

   myparall->Brdcst_Values(0,0,nbasis*n1dgrid,phi_ae);
   myparall->Brdcst_Values(0,0,nbasis*n1dgrid,dphi_ae);
   myparall->Brdcst_Values(0,0,nbasis*n1dgrid,phi_ps);
   myparall->Brdcst_Values(0,0,nbasis*n1dgrid,dphi_ps);
   myparall->Brdcst_Values(0,0,nbasis*n1dgrid,prj_ps);
   myparall->Brdcst_Values(0,0,nbasis*n1dgrid,prj_ps0);

   myparall->Brdcst_Values(0,0,n1dgrid,core_ae);
   myparall->Brdcst_Values(0,0,n1dgrid,core_ps);
   myparall->Brdcst_Values(0,0,n1dgrid,v_ps);

   myparall->Brdcst_Values(0,0,1,&sigma);
   myparall->Brdcst_Values(0,0,1,&zion);


   // define nprj and lmax 
   locp = -1;
   lmax = -1;
   nprj = 0;
   for (auto ii=0; ii<nbasis; ++ii)
   {
      auto l = lps[ii];
      nprj += 2*l+1;
      if (l>lmax) lmax = l;
   }

   // define nmax
   nmax = 0;
   {
      int nmaxl[lmax+1];
      for (auto i=0; i<(lmax+1); ++i) nmaxl[i] = 0;
      for (auto ii=1; ii<nbasis; ++ii)
         nmaxl[lps[ii]] += 1;
      for (auto l=0; l<=lmax; ++l)
         if (nmaxl[l]>nmax) 
            nmax = nmaxl[l];
   }


   // allocate Gijl,Sijl,Tijl,vcore,Vpseuo 
   //Gijl               = new double[nmax*nmax*(lmax+1)];
   //hartree_matrix     = new double[nbasis*nbasis*nbasis*nbasis*(2*lmax+1)];
   //comp_charge_matrix = new double[nbasis*nbasis*(2*lmax+1)];
   //comp_pot_matrix    = new double[nbasis*nbasis*(2*lmax+1)];


   // define n_prj,l_prj, and m_prj
   n_prj = new int[nprj];
   l_prj = new int[nprj];
   m_prj = new int[nprj];
   b_prj = new int[nprj];

   int lcount = 0;
   for (auto i=0; i<nbasis; ++i)
   {
      int la=lps[i];
      int na=nps[i] - la;

      // f-wave
      if (la==3)
      {
         n_prj[lcount] = na;
         l_prj[lcount] = 3;
         m_prj[lcount] = -3;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 3;
         m_prj[lcount] = -2;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 3;
         m_prj[lcount] = -1;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 3;
         m_prj[lcount] = 0;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 3;
         m_prj[lcount] = 1;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 3;
         m_prj[lcount] = 2;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 3;
         m_prj[lcount] = 3;
         b_prj[lcount] = i;
         ++lcount;
      }

      // d-wave
      if (la==2)
      {
         n_prj[lcount] = na;
         l_prj[lcount] = 2;
         m_prj[lcount] = -2;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 2;
         m_prj[lcount] = -1;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 2;
         m_prj[lcount] = 0;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 2;
         m_prj[lcount] = 1;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 2;
         m_prj[lcount] = 2;
         b_prj[lcount] = i;
         ++lcount;
      }

      // p-wave
      if (la==1)
      {
         n_prj[lcount] = na;
         l_prj[lcount] = 1;
         m_prj[lcount] = -1;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 1;
         m_prj[lcount] = 0;
         b_prj[lcount] = i;
         ++lcount;

         n_prj[lcount] = na;
         l_prj[lcount] = 1;
         m_prj[lcount] = 1;
         b_prj[lcount] = i;
         ++lcount;
      }

      // s-wave
      if (la==0)
      {
         n_prj[lcount] = na;
         l_prj[lcount] = 0;
         m_prj[lcount] = 0;
         b_prj[lcount] = i;
         //vnl(1,1,1,lcount)=vnl_ray(1,i,1) !*** j0!=0 at G=0***
         ++lcount;
      }
   }





}


/****************************************************
 *                                                  *
 *     Psp1d_pawppv1::vpp_generate_paw_matrices     *
 *                                                  *
 ****************************************************/
void Psp1d_pawppv1::vpp_generate_paw_matrices(Parallel *myparall, double *Gijl, double *comp_charge_matrix, double *comp_pot_matrix, double *hartree_matrix)
{

   double twopi = 8.0*atan(1.0);
   int nb1 = nbasis;
   int nb2 = nbasis*nbasis;
   int nb3 = nbasis*nbasis*nbasis;
   int nb4 = nbasis*nbasis*nbasis*nbasis;

   double *f1 = new double[n1dgrid];
   double *f2 = new double[n1dgrid];
   double *f3 = new double[n1dgrid];
   double *f4 = new double[n1dgrid];

   // **** comp_charge_matrix(nbasis,nbasis,0:2*lmax) ****
   for (auto i=0; i<nbasis; ++i)
      for (auto j=0; j<=i; ++j)
      {
         for (auto k=0; k<icut; ++k)
            f1[k]=phi_ae[k+i*n1dgrid]*phi_ae[k+j*n1dgrid]-phi_ps[k+i*n1dgrid]*phi_ps[k+j*n1dgrid];

         for (auto l=0; l<(2*lmax+1); ++l)
         {
            double d = util_log_integrate_def(2*l+2,f1,l,rgrid,log_amesh,icut);
            comp_charge_matrix[i+j*nb1+l*nb2] = d;
            comp_charge_matrix[j+i*nb1*l*nb2] = d;
         }
      }

   // **** comp_pot_matrix(nbasis,nbasis,0:2*lmax) ****
   for (auto l=0; l<(2*lmax+1); ++l)
   {
      int k1 = l+2;
      int k2 = 2*l+2;

      util_compcharge_gen_rgaussian(l,sigma,n1dgrid,rgrid,f1);
      for (auto k=0; k<icut; ++k)
         f1[k] *= rgrid[k]*rgrid[k];
 
      for (auto i=0; i<nbasis; ++i)
      {
         for (auto j=0; j<=i; ++j)
         {
            for (auto k=0; k<icut; ++k)
               f2[k] = phi_ps[k+i*n1dgrid]*phi_ps[k+j*n1dgrid];

            double d = util_log_multipole_energy(l,icut,rgrid,k1,f1,k2,f2,log_amesh);
            comp_pot_matrix[i+j*nb1+l*nb2] = d;
            comp_pot_matrix[j+i*nb1+l*nb2] = d;
         }
      }
   }

   // **** hartree_matrix(nbasis,nbasis,nbasis,nbasis,0:2*lmax)        ****
   // **** Note - This is the effective hartree matrix which includes  ****
   // **** comp_charge_matrix and comp_pot_matrix terms in it.         ****

   memset(hartree_matrix,0,nbasis*nbasis*nbasis*nbasis*(2*lmax+1)*sizeof(double));

   for (auto i=0; i<nbasis; ++i)
      for (auto j=0; j<=i; ++j)
      {
         int sum_l = lps[i]+lps[j];
         int dif_l = abs(lps[i]-lps[j]);
         int k1    = sum_l + 2;

         for (auto k=0; k<icut; ++k)
         {
            f1[k] = phi_ae[k+i*n1dgrid]*phi_ae[k+j*n1dgrid];
            f3[k] = phi_ps[k+i*n1dgrid]*phi_ps[k+j*n1dgrid];
         }

         for (auto i2=0; i2<nbasis; ++i2)
             for (auto j2=0; j2<=i2; ++j2)
             {
                int sum_l2 = lps[i2]+lps[j2];
                int dif_l2 = abs(lps[i2]-lps[j2]);
                int k2    = sum_l2 + 2;

                for (auto k=0; k<icut; ++k)
                {
                   f2[k] = phi_ae[k+i2*n1dgrid]*phi_ae[k+j2*n1dgrid];
                   f4[k] = phi_ps[k+i2*n1dgrid]*phi_ps[k+j2*n1dgrid];
                }

                for (auto l=0; l<(2*lmax+1); ++l)
                {
                   double d   = ((double) (2*l+1)) * util_double_factorial(2*l+1)*pow(sigma,(2*l+1));
                   double vgl = 4.00*sqrt(twopi)/d;

                   if ((l<=sum_l) && (l>=dif_l) && (l<=sum_l2) && (l>=dif_l2)) 
                   {
                       d = util_log_multipole_energy(l,icut,rgrid,k1,f1,k2,f2,log_amesh)
                         - util_log_multipole_energy(l,icut,rgrid,k1,f3,k2,f4,log_amesh);

                       hartree_matrix[i + j*nb1 + i2*nb2 + j2*nb3 + l*nb4] = d
                               - 2.00*    comp_pot_matrix[i  + j *nb1 + l*nb2]
                                      *comp_charge_matrix[i2 + j2*nb1 + l*nb2]
                               - vgl  *comp_charge_matrix[i  + j *nb1 + l*nb2]
                                      *comp_charge_matrix[i2 + j2*nb1 + l*nb2];
                       hartree_matrix[j + i*nb1 + i2*nb2 + j2*nb3 + l*nb4] = d
                               - 2.00*    comp_pot_matrix[j  + i *nb1 + l*nb2]
                                      *comp_charge_matrix[i2 + j2*nb1 + l*nb2]
                               - vgl  *comp_charge_matrix[j  + i *nb1 + l*nb2]
                                      *comp_charge_matrix[i2 + j2*nb1 + l*nb2];
                       hartree_matrix[i + j*nb1 + j2*nb2 + i2*nb3 + l*nb4] = d
                               - 2.00*    comp_pot_matrix[i  + j *nb1 + l*nb2]
                                      *comp_charge_matrix[j2 + i2*nb1 + l*nb2]
                               - vgl  *comp_charge_matrix[i  + j *nb1 + l*nb2]
                                      *comp_charge_matrix[j2 + i2*nb1 + l*nb2];

                       hartree_matrix[j + i*nb1 + j2*nb2 + i2*nb3 + l*nb4] = d
                               - 2.00*    comp_pot_matrix[j  + i *nb1 + l*nb2]
                                      *comp_charge_matrix[j2 + i2*nb1 + l*nb2]
                               - vgl  *comp_charge_matrix[j  + i *nb1 + l*nb2]
                                      *comp_charge_matrix[j2 + i2*nb1 + l*nb2];
                   }
                }

             }
      }

}


/*******************************************
 *                                         *
 *     Psp1d_pawppv1::vpp_generate_ray     *
 *                                         *
 *******************************************/
void Psp1d_pawppv1::vpp_generate_ray(Parallel *myparall, int nray, double *G_ray, double *vl_ray, double *vlpaw_ray, double *vnl_ray)
{

   double pi    = 4.00*atan(1.0);
   double twopi = 2.0*pi;
   double forpi = 4.0*pi;

   double p0 = sqrt(forpi);
   double p1 = sqrt(3.0*forpi);
   double p2 = sqrt(15.0*forpi);
   double p3 = sqrt(105.0*forpi);

   double zero = 0.0;
   int    izero = 0;
   int    ione  = 1;
   int    nray2 = 2*nray;
   

   double q;
   double *cs = new double[n1dgrid];
   double *sn = new double[n1dgrid];
   double *f  = new double[n1dgrid];
   double a,xx;

   memset(vl_ray,0,nray*sizeof(double));
   memset(vlpaw_ray,0,nray*sizeof(double));
   memset(vnl_ray,0,nbasis*nray*sizeof(double));

   for (auto k1=(1+myparall->taskid()); k1<nray; k1+=myparall->np())
   {
      q=G_ray[k1];
      for (auto i=0; i<n1dgrid; ++i)
      {
         cs[i]=cos(q*rgrid[i]);
         sn[i]=sin(q*rgrid[i]);
      }
      for (auto ib=0; ib<nbasis; ++ib)
      {
         int la=lps[ib];

         /* h projectors */
         /* g projectors */
         // ::::::::::::::::::::::::::::::  f-wave  ::::::::::::::::::::::::::::::
         if (la==3)
         {
            int power_f=8;
            f[0]=0.0;
            for (auto i=1; i<n1dgrid; ++i)
            {
               a=sn[i]/(q*rgrid[i]);
               a=15.00*(a-cs[i])/pow((q*rgrid[i]),2.0) - 6.0*a + cs[i];
               f[i]=a*prj_ps[i+ib*n1dgrid];
             }
             vnl_ray[k1+ib*nray]= p3*util_log_integrate_def(power_f,f,0,rgrid,log_amesh,n1dgrid)/q;
         }
         // ::::::::::::::::::::::::::::::  d-wave  ::::::::::::::::::::::::::::::
         if (la==2) 
         {
             int power_f=6;
             f[0]=0.0;
             for (auto i=0; i<n1dgrid; ++i)
             {
               a=3.00*(sn[i]/(q*rgrid[i])-cs[i])/(q*rgrid[i])-sn[i];
               f[i]=a*prj_ps[i+ib*n1dgrid];
             }
             vnl_ray[k1+ib*nray]=p2*util_log_integrate_def(power_f,f,0,rgrid,log_amesh,n1dgrid)/q;
         }
         // ::::::::::::::::::::::::::::::  p-wave  ::::::::::::::::::::::::::::::
         if (la==1) 
         {
             int power_f=4;
             f[0]=0.0;
             for (auto i=0; i<n1dgrid; ++i)
             {
               f[i]=(sn[i]/(q*rgrid[i])-cs[i])*prj_ps[i+ib*n1dgrid];
             }
             vnl_ray[k1+ib*nray]=p1*util_log_integrate_def(power_f,f,0,rgrid,log_amesh,n1dgrid)/q;
         }
         // ::::::::::::::::::::::::::::::  s-wave  ::::::::::::::::::::::::::::::
         if (la==0) 
         {
            int power_f=2;
            for (auto i=0; i<n1dgrid; ++i)
            {
              f[i]=sn[i]*prj_ps[i+ib*n1dgrid];
            }
            vnl_ray[k1+ib*nray]=p0*util_log_integrate_def(power_f,f,0,rgrid,log_amesh,n1dgrid)/q;
         }
         
      }

      // **** local ****
      if (version==3)
      {
         vl_ray[k1]=-zv*(forpi/(q*q))*exp(-0.25*(q*sigma)*(q*sigma));
      }
      else if (version==4)
      {
           vl_ray[k1]=-zv*(forpi/(q*q))*(exp(-0.25*(q*sigma)*(q*sigma))
                                        -exp(-0.25*(q*rlocal)*(q*rlocal)));
      }

      // **** vlpaw potential ****
      for (auto i=0; i<n1dgrid; ++i)
          f[i]=(rgrid[i]*v_ps[i])*sn[i];
      vlpaw_ray[k1]=forpi*util_log_integrate_def(0,f,0,rgrid,log_amesh,n1dgrid)/q;
   }

   myparall->Vector_SumAll(0,nray,vl_ray);
   myparall->Vector_SumAll(0,nray,vlpaw_ray);
   myparall->Vector_SumAll(0,nbasis*nray,vnl_ray);

   // :::::::::::::::::::::::::::::::  G=0  ::::::::::::::::::::::::::::::::
   if (version==3)
   {
      vl_ray[0] = zv*pi*sigma*sigma;
   }
   else if (version==4) 
   {
      vl_ray[0] = zv*pi*(sigma*sigma-rlocal*rlocal);
   }

   for (auto i=0; i<n1dgrid; ++i)
      f[i]=v_ps[i]*rgrid[i]*rgrid[i];
   vlpaw_ray[0]=forpi*util_log_integrate_def(0,f,0,rgrid,log_amesh,n1dgrid);

   // *** only j0 is non-zero at zero ****
   int power_f=2;
   for (auto ib=0; ib<nbasis; ++ib)
   {
      int la=lps[ib];
      vnl_ray[ib*nray]=0.0;
      if (la==0)
      {
         for (auto i=0; i<n1dgrid; ++i)
            f[i]=rgrid[i]*prj_ps[i+ib*n1dgrid];
         vnl_ray[ib*nray]=p0*util_log_integrate_def(power_f,f,0,rgrid,log_amesh,n1dgrid);
      }
   }

   delete [] f;
   delete [] sn;
   delete [] cs;
}



/*******************************************
 *                                         *
 *   Psp1d_pawppv1::vpp_generate_spline    *
 *                                         *
 *******************************************/
void Psp1d_pawppv1::vpp_generate_spline(PGrid *mygrid, int nray, double *G_ray, double *vl_ray, double *vlpaw_ray, double *vnl_ray,
                                        double *vl, double *vlpaw, double *vnl)
{
   double pi    = 4.00*atan(1.0);

   /* allocate spline grids */
   double *vl_splineray    = new double [nray];
   double *vlpaw_splineray = new double [nray];
   double *vnl_splineray   = new double [nbasis*nray];
   double *tmp_splineray   = new double [nray];


   /* setup cubic bsplines */
   double dG = G_ray[2] - G_ray[1];

   /* five point formula */
   double yp1 = ( -50.00*vl_ray[1]
                 + 96.00*vl_ray[2]
                 - 72.00*vl_ray[3]
                 + 32.00*vl_ray[4]
                 -  6.00*vl_ray[5])/(24.00*dG);
   util_spline(&(G_ray[1]),&(vl_ray[1]),nray-1,yp1,0.00,&(vl_splineray[1]),tmp_splineray);
   yp1 = ( -50.00*vlpaw_ray[1]
          + 96.00*vlpaw_ray[2]
          - 72.00*vlpaw_ray[3]
          + 32.00*vlpaw_ray[4]
          -  6.00*vlpaw_ray[5])/(24.00*dG);
   util_spline(&(G_ray[1]),&(vlpaw_ray[1]),nray-1,yp1,0.00,&(vlpaw_splineray[1]),tmp_splineray);

   for (auto i=0; i<nbasis; ++i)
      util_spline(G_ray,&(vnl_ray[i*nray]),nray,0.00,0.00,&(vnl_splineray[i*nray]),tmp_splineray);


   double q, qx, qy, qz, xx;
   double *gx, *gy, *gz;
   int npack0 = mygrid->npack(0);
   int npack1 = mygrid->npack(1);
   int nx,lcount;
   mygrid->t_pack_nzero(0,1,vl);
   mygrid->t_pack_nzero(1,nbasis,vnl);

   // **** generate vl and vlpaw ****
   gx = mygrid->Gpackxyz(0,0);
   gy = mygrid->Gpackxyz(0,1);
   gz = mygrid->Gpackxyz(0,2);
   for (auto k=0; k<npack0; ++k)
   {
      qx = gx[k]; qy = gy[k]; qz = gz[k];
      q = sqrt(qx*qx + qy*qy + qz*qz);
      nx = (int) floor(q/dG);

      if (q>1.0e-9)
      {
         qx /= q; qy /= q; qz /= q;
         vl[k] = util_splint(&(G_ray[1]),&(vl_ray[1]),&(vl_splineray[1]),nray-1,nx,q);
         vlpaw[k] = util_splint(&(G_ray[1]),&(vl_ray[1]),&(vl_splineray[1]),nray-1,nx,q);
      }
      else
      {
         vl[k] = vl_ray[0];
         vlpaw[k] = vl_ray[0];
      }
   }

   /* generate vnl */
   gx = mygrid->Gpackxyz(1,0);
   gy = mygrid->Gpackxyz(1,1);
   gz = mygrid->Gpackxyz(1,2);

   for (auto k=0; k<npack1; ++k)
   {
      qx = gx[k]; qy = gy[k]; qz = gz[k];
      q = sqrt(qx*qx + qy*qy + qz*qz);
      nx = (int) floor(q/dG);

      if (q>1.0e-9)
      {
         qx /= q; qy /= q; qz /= q;
         lcount=0;
         for (auto i=0; i<nbasis; ++i)
         {
            int la=lps[i];

            // **** f projectors ****
            if (la==3)
            {
               xx = util_splint(G_ray,&(vnl_ray[i*nray]),&(vnl_splineray[i*nray]),nray,nx,q);
               vnl[k+lcount*npack1]=  xx * qy*(3.00*(1.00-qz*qz)-4.00*qy*qy)/sqrt(24.00);
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qx*qy*qz;
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qy*(5.00*qz*qz-1.00)/sqrt(40.00);
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qz*(5.00*qz*qz-3.00)/sqrt(60.00);
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qx*(5.00*qz*qz-1.00)/sqrt(40.00);
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qz*(qx*qx - qy*qy)/2.00;
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qx*(4.00*qx*qx-3.00*(1.00-qz*qz))/sqrt(24.00);
               ++lcount;
            }

            // **** d projectors ****
            if (la==2)
            {
               xx = util_splint(G_ray,&(vnl_ray[i*nray]),&(vnl_splineray[i*nray]),nray,nx,q);
               vnl[k+lcount*npack1]=  xx * qx*qy;
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qy*qz;
               ++lcount;
               vnl[k+lcount*npack1]=  xx * (3.00*qz*qz-1.00)/(2.00*sqrt(3.00));
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qz*qx;
               ++lcount;
               vnl[k+lcount*npack1]=  xx * (qx*qx-qy*qy)/(2.00);
               ++lcount;
            }

            // **** p projectors ****
            if (la==1)
            {
               xx = util_splint(G_ray,&(vnl_ray[i*nray]),&(vnl_splineray[i*nray]),nray,nx,q);
               vnl[k+lcount*npack1]=  xx * qy;
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qz;
               ++lcount;
               vnl[k+lcount*npack1]=  xx * qx;
               ++lcount;
            }

            /* s projectors */
            if (la==0)
            {
               xx = util_splint(G_ray,&(vnl_ray[i*nray]),&(vnl_splineray[i*nray]),nray,nx,q);
               vnl[k+lcount*npack1]=  xx;
               ++lcount;
            }
         }
      }
   }

   delete [] tmp_splineray;
   delete [] vnl_splineray;
   delete [] vlpaw_splineray;
   delete [] vl_splineray;

}

}


