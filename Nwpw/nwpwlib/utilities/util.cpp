
//extern "C" {
//#include        "compressed_io.h"
//}
#include        "compressed_io.hpp"

#include	<iostream>
#include	<cstdlib>
#include	<cmath>
#include        "Parallel.hpp"
#include	"util.hpp"

namespace pwdft {
using namespace pwdft;

void c_aindexcopy(const int n, const int *indx, double *A, double *B)
{
   int ii,jj;
   ii = 0;
   for (int i=0; i<n; ++i)
   {
      jj      = 2*indx[i];
      B[ii]   = A[jj];
      B[ii+1] = A[jj+1];
      ii += 2;
   }
}

void c_bindexcopy(const int n, const int *indx, double *A, double *B)
{
   int ii,jj;
   ii = 0;
   for (int i=0; i<n; ++i)
   {
      jj      = 2*indx[i];
      B[jj]   = A[ii];
      B[jj+1] = A[ii+1];
      ii += 2;
   }
}
void c_bindexcopy_conjg(const int n, const int *indx, double *A, double *B)
{
   int ii,jj;
   ii = 0;
   for (int i=0; i<n; ++i)
   {
      jj      = 2*indx[i];
      B[jj]   = A[ii];
      B[jj+1] = -A[ii+1];
      ii += 2;
   }
}


void t_aindexcopy(const int n, const int *indx, double *A, double *B)
{
   for (int i=0; i<n; ++i)
      B[i] = A[indx[i]];
}

void t_bindexcopy(const int n, const int *indx, double *A, double *B)
{
   for (int i=0; i<n; ++i)
      B[indx[i]] = A[i];
}


void i_aindexcopy(const int n, const int *indx, int *A, int *B)
{
   for (int i=0; i<n; ++i)
      B[i] = A[indx[i]];
}




/**************************************
 *                                    *
 *            eigsrt                  *
 *                                    *
 **************************************/

void eigsrt(double *D, double *V, int n)
{
   int i,j,k;
   double p;

   for (i=0; i<(n-1); ++i)
   {
      k = i;
      p = D[i];
      for(j=i+1; j<n; ++j)
         if (D[j]>=p)
         {
            k = j;
            p = D[j];
         }

      if (k!=i)
      {
         D[k] = D[i];
         D[i] = p;
         for (j=0; j<n; ++j)
         {
            p = V[j+i*n];
            V[j+i*n] = V[j+k*n];
            V[j+k*n] = p;
         }
      }
   }
}



/**************************************
 *                                    *
 *         util_getfilling            *
 *                                    *
 **************************************/
void util_getfilling(int f, int nfft[], int *filling, double zvalue[])
{
   int h = (f%2);
   int f2 = (f-h)/2;
   int k=0;
   while ((k+1)*(k+2)*(k+3) <= (6*f2))
      ++k;

   f2 -= k*(k+1)*(k+2)/6;
   int j=0;
   while ((j+1)*(j+2) <= (2*f2))
      ++j;

   int i = f2 - j*(j+1)/2;

   filling[0] = i;
   if (i==0) 
   {
      filling[1] = j;
      filling[2] = k;
   }
   else
   {

      if (((j+1)%2)==0)
         filling[1] = (j+1)/2;
      else 
         filling[1] = -(j+1)/2;

      if (((k+1)%2)==0)
         filling[2] = (k+1)/2;
      else 
         filling[2] = -(k+1)/2;
   }

   if ((i==0) && (j==0) && (k==0))
   {
      filling[3] = 0;
      zvalue[0] = 1.0;
      zvalue[1] = 0.0;
   }
   else if (h==0)
   {
      filling[3] = 2;
      zvalue[0] = 1.0/sqrt(2.0);
      zvalue[1] = 0.0;
   }
   else
   {
      filling[3] = -2;
      zvalue[0] = 0.0;
      zvalue[1] = 1.0/sqrt(2.0);
   }

   /* modularize the filling */
   int inc2c =( nfft[0]/2 + 1);
   filling[0] = (filling[0]+inc2c)%inc2c;
   filling[1] = (filling[1]+nfft[1])%nfft[1];
   filling[2] = (filling[2]+nfft[2])%nfft[2];
}


/**************************************
 *                                    *
 *           util_random              *
 *                                    *
 **************************************/

/* returns a random number between 0 and 1

   Entry - seed - if zero set a seed with srand
   Exit - returns a random number between 0 and 1
   Uses - rand and srand stdlib functions
*/
double util_random(const int seed)
{
   if (seed>0) std::srand(((double) seed));
   return ( (double) std::rand()/RAND_MAX);
}



/**************************************
 *                                    *
 *        util_double_factorial       *
 *                                    *
 **************************************/
double util_double_factorial(const int n)
{
   int n11[18] = {1,1,1,2,3,8,15,48,105,384,945,3840,10395,46080,135135,645120,2027025,10321920};

   if ((n>=-1) && (n<=16))
   {
      return ((double) n11[n+1]);
   }
   else
   {
      std::cout << "too big parameter in nwpw_double_factorial" << std::endl;
      return  -999999.0;
   }
}

/*************************************************
 *                                               *
 *         util_compcharge_gen_rgaussian         *
 *                                               *
 *************************************************/
void util_compcharge_gen_rgaussian(const int l, const double sigma, const int nr, const double r[], double gl[])
{
   double pi = 4.00*atan(1.00);
   double c  = pow(2.00,(l+2))/(sqrt(pi)*(util_double_factorial(2*l+1))*pow(sigma,(2*l+3)));

   // *** this fixes possible underflow error ***
   for (auto i=0; i<nr; ++i) gl[i] = 0.00;

   for (auto i=0; i<nr; ++i)
      if (abs(r[i]) < (8.00*sigma))
         gl[i] = c*pow(r[i],l)*exp(-(r[i]/sigma)*(r[i]/sigma));
}


/*******************************************
 *                                         *
 *         util_log_integrate_def          *
 *                                         *
 *******************************************/
double util_log_integrate_def(const int power_f, const double f[],
                              const int power_r, const double r[],
                              const double log_amesh, const int  nrange)
{
   double integrand[nrange];
   //double *integrand = new double[nrange];

   for (auto k=0; k<nrange; ++k)
       integrand[k] = f[k]*pow(r[k],power_r+1);

   // *** integrate from the origin to the first point ***
   double sum_f = integrand[0]/((double) (power_r+power_f+1));

   // *** the rest via trapesoidal rule ***
   double tmp_sum = 0.0;
   for (auto k=0; k<nrange; ++k)
      tmp_sum += integrand[k];

   // *** the rest via trapesoidal rule ***
   sum_f +=  log_amesh*tmp_sum - 0.50*log_amesh*(integrand[0]+integrand[nrange-1]);

   //delete [] integrand;

   return sum_f;
}


/************************************************
 *                                              *
 *            log_integrate_def0                *
 *                                              *
 ************************************************/
/*     
     Computes the following integral
      
      
          /               /
          | f(r)/r dr =   | f(i) * log(a) di
          /               /
 
     where r(i) = r0*a**i   and f(r-->0) = r**power_f
*/

double util_log_integrate_def0(const int power_f,
                               const double f[],
                               const double r[],
                               const double log_amesh,
                               const int nrange)
{
   // *** integrate from the origin to the first point ***
   double sum_f = f[0]/((double) power_f);

   // *** the rest via trapesoidal rule ***
   double tmp_sum = 0.0;
   for (auto i=0; i<nrange; ++i)
      tmp_sum += f[i];

   // *** the rest via trapesoidal rule ***
   sum_f += log_amesh*tmp_sum - 0.50*log_amesh*(f[0] + f[nrange-1]);

   return sum_f;
}



/************************************************
 *                                              *
 *            log_integrate_indef               *
 *                                              *
 ************************************************/
void util_log_integrate_indef(const int power_f, const double f[],
                              const int power_r, const double r[],
                              const double  log_amesh,
                              const int nrange, double sum_f[])
{
   double integrand[nrange];
   //double *integrand = new double[nrange];

   for (auto k=0; k<nrange; ++k)
      integrand[k] = f[k]*pow(r[k],(power_r+1));

   if (nrange<=5)
   {
      for (auto k=0; k<nrange; ++k)
         sum_f[k] = integrand[k]/((double) (power_r+power_f+1));
   }
   else
   {
      for (auto k=0; k<5; ++k)
         sum_f[k] = integrand[k]/((double) (power_r+power_f+1));

      for (auto k=5; k<nrange; ++k)
          sum_f[k] = sum_f[k-1] + log_amesh*0.50*(integrand[k-1] + integrand[k]);
   }
   //delete [] integrand;
}


/*******************************************
 *                                         *
 *         util_log_multipole_energy       *
 *                                         *
 *******************************************/
double util_log_multipole_energy(const int l, const int nrange, const double g_r[],
                                 const int power_q1, const double q1[],
                                 const int power_q2, const double q2[],
                                 const double log_amesh)
{
   double fourpi  = 16.0*atan(1.0);
   double power_f = power_q1 + power_q2;
   double q1_l[nrange],q2_l[nrange], ftmp[nrange];

   util_log_integrate_indef(power_q1,q1,l,g_r,log_amesh,nrange,q1_l);
   util_log_integrate_indef(power_q1,q2,l,g_r,log_amesh,nrange,q2_l);

   for (auto k=0; k<nrange; ++k)
      ftmp[k] = (q1[k]*q2_l[k]+q1_l[k]*q2[k])/pow(g_r[k],(l+1));

   return (util_log_integrate_def(power_f,ftmp,0,g_r,log_amesh,nrange) * fourpi/((double) (2*l+1)));
}


/************************************************
 *                                              *
 *           util_log_r2integrate_eric          *
 *                                              *
 ************************************************/
double util_log_r2integrate_eric(const int ngrid, const double log_amesh, const double rgrid[], const double f[])
{
   double mysum = (   9.00*f[0]*pow(rgrid[0],3)
                   + 28.00*f[1]*pow(rgrid[1],3)
                   + 23.00*f[2]*pow(rgrid[2],3)
                  )/24.00;

   for (auto i=3; i<ngrid; ++i)
      mysum += f[i]*pow(rgrid[i],3);

   return (log_amesh*mysum + f[0]*pow(rgrid[0],3)/3.00);
}

/************************************************
 *                                              *
 *            log_corrector_iF                  *
 *                                              *
 ************************************************/
/*
     Computes y(i) <-- y(i+1) + F(f(i),f(i+1),f(i+2),f(i+3))
    where F is a 5 point corrector                           */

double util_log_corrector_iF(const int i, const double f[])
{
   double oneo24 = 1.0/24.0;

   return ( -oneo24 * (   9.00*f[i]
                       + 19.00*f[i+1]
                       -  5.00*f[i+2]
                       +  1.00*f[i+3]));
}


/************************************************
 *                                              *
 *            util_log_coulomb0_energy          *
 *                                              *
 ************************************************/
double util_log_coulomb0_energy(const double rho[],     const double charge,
                                const double rgrid[],   const int ngrid,
                                const double log_amesh, const double zion)
{
   double fourpi = 16.00*atan(1.0);
   double tmp[ngrid];

   for (auto i=0; i<ngrid; ++i)
      tmp[i] = fourpi*rho[i]/rgrid[i];

   double E  = util_log_r2integrate_eric(ngrid,log_amesh,rgrid,tmp);

   E  = -E*zion;

   return E;
}


/************************************************
 *                                              *
 *            util_log_coulomb_energy           *
 *                                              *
 ************************************************/
double util_log_coulomb_energy(const double rho[],   const double charge,
                               const double rgrid[], const int ngrid, 
                               const double log_amesh)
{
   double fourpi = 16.0*atan(1.0);
   double tmp[ngrid],vh[ngrid];

   for (auto i=0; i<ngrid; ++i)      
      tmp[i] = fourpi*rho[i]*log_amesh*pow(rgrid[i],3);

   vh[ngrid-1] = charge;
   vh[ngrid-2] = charge;
   vh[ngrid-3] = charge;

   for (auto i=(ngrid-4); i>=0; --i)
      vh[i] = vh[i+1] + util_log_corrector_iF(i,tmp);

   for (auto i=0; i<ngrid; ++i)
      tmp[i] = fourpi*rho[i]*log_amesh*pow(rgrid[i],2);
    
   double tt = 0.0;

   for (auto i=(ngrid-4); i>=0; --i)
   {
      tt += util_log_corrector_iF(i,tmp);
      vh[i] -= rgrid[i]*tt;
   }

   for (auto i=0; i<ngrid; ++i)
      vh[i] /= rgrid[i];

   for (auto i=0; i<ngrid; ++i)
      tmp[i] = fourpi*rho[i]*vh[i];

   double E = util_log_r2integrate_eric(ngrid,log_amesh,rgrid,tmp);

   return (E*0.50);
}



/**************************************
 *                                    *
 *           util_filefind            *
 *                                    *
 **************************************/
bool util_filefind(Parallel *myparall, char *fname)
{
   int ifound;

   if (myparall->is_master())
      ifound = cfileexists(fname);

   myparall->Brdcst_iValue(0,0,&ifound);

   return (ifound>0);
}


/**************************************
 *                                    *
 *           util_spline              *
 *                                    *
 **************************************/
void   util_spline(double *x, double *y, int n, double yp1, double ypn, double *y2, double *utmp)
{
   double sig,p,qn,un;

   if (yp1>0.99e30)
   {
      y2[0]   = 0.0;
      utmp[0] = 0.0;
   }
   else
   {
      y2[0]   = -0.5;
      utmp[0] = 3.0 / (x[1]-x[0]) * ( (y[1]-y[0]) / (x[1]-x[0]) - yp1);
   }
   for (auto i=1; i<(n-1); ++i)
   {
      sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
      p = sig*y2[i-1] + 2.00;
      y2[i] = (sig-1.00) / p;
      utmp[i] = (
                 6.00 *
                 (
                   (y[i+1]-y[i  ]) / (x[i+1]-x[i  ])
                 - (y[i  ]-y[i-1]) / (x[i  ]-x[i-1])
                 )
                 / (x[i+1]-x[i-1]) - sig*utmp[i-1]
               )
               / p;
   }
   
   if (ypn>0.99e30)
   {
      qn = 0.0;
      un = 0.0;
   }
   else
   {
      qn = 0.5;
      un = 3.00 / (x[n-1]-x[n-2])
              * ( ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]) );
   
   }

   y2[n-1] = (un-qn*utmp[n-2]) / (qn*y2[n-2]+1.00);
   for (auto k=n-2; k>=0; --k)
      y2[k] = y2[k]*y2[k+1] + utmp[k];

}



/**************************************
 *                                    *
 *           util_splint              *
 *                                    *
 **************************************/
double util_splint(double *xa, double *ya, double *y2a, int n, int nx, double x)
{
   int khi = nx;
   int klo = nx-1;

   while ((xa[klo]>x) || (xa[khi]<x))
   {
      if (xa[klo]>x)
      {
        klo = klo - 1;
        khi = khi - 1;
      }
      if (xa[khi]<x)
      {
           klo = klo + 1;
           khi = khi + 1;
      }
   }

   double h = xa[khi]-xa[klo];
   double a = (xa[khi]-x)/h;
   double b = (x-xa[klo])/h;
   
   return  (a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * h*h / 6.0);
}



/**************************************
 *                                    *
 *           util_filter              *
 *                                    *
 **************************************/
void util_filter(int nray, double *g_ray, double ecut, double *v_ray)
{
   int    ncut = 15;
   double qmax = sqrt(ecut+ecut);
   double g;

   for (auto i=0; i<nray; ++i)
   {
      g = g_ray[i];
      if (g > (qmax-0.2))
      {
         v_ray[i] *= ( 1.0 - pow((1.0-exp(-pow((g/qmax),ncut))),ncut) );
      }
   }
}



#ifdef _MATPRINT_
void util_matprint(std::string matlabel, int n, double *A) {
   std::cout << "util_matprint: " << matlabel << std::endl;
   for (int i=0; i<n; ++i)
   {
      for (int j=0; j<n; ++j) std::cout << A[i+j*n] << " ";
      std::cout << std::endl;
   }
   std::cout << std::endl;
}
#endif

}

