
//#include      <iomanip>
#include        <iostream>
#include        <cstdlib>
#include        <complex>
#include        <cmath>

#include        "util_gaunt.hpp"
#include        "util_tesseral.hpp"
#include        "util_log_integrate.hpp"
#include        "util_wgaussian.hpp"

namespace pwdft {



/******************************************************
 *                                                    *
 *                 util_WGaussian                     *
 *                                                    *
 ******************************************************/
/*
     Calculates the two electron two center Gaussian integral

                                            //
    WGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                            || ------------------------------------  dr dr'
                                            //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
            /
*/
double util_WGaussian(const int la, const int ma, const double sa, 
                      const int lb, const int mb, const double sb, const double Rab[])
{
   double pi = 4.0*atan(1.0);
   double x =  util_double_factorial(2*la+1);
   double y =  util_double_factorial(2*lb+1);
   double alpha = sqrt(0.25*(sa*sa + sb*sb));
   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;

   double phi,c,fac;

   if ((std::abs(Rab[1])<1.0e-9) && (std::abs(Rab[0])<1.0e-9))
   {
      phi = 0.0;
   }
   else
   {
      phi = atan2(Rab[1],Rab[0]);
   }

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }
   

   if ((((std::abs(la-lb)+la+lb)/2)%2)==1)
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   double tmp = 0.0;
   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      double mtmp = 0.0;
      for (auto m=-l; m<=l; ++m)
      {
        mtmp += util_gaunt(false,l,m,la,ma,lb,mb)*util_Tesseral_lm(l,m,cos_theta,phi);
      }
      tmp += fac*mtmp*util_GaussBessel(la+lb,l,alpha,R);
      fac = -fac;
    }
    
   return (c*tmp);
}


/******************************************************
 *                                                    *
 *             util_dWGaussian                        *
 *                                                    *
 ******************************************************/
/*
     Calculates the two electron two center Gaussian integral and it's derivative wrt to Rab

                                            //
   dWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                            || ------------------------------------  dr dr'
                                            //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
            /
*/
void util_dWGaussian(const int la, const int ma, const double sa, 
                     const int lb, const int mb, const double sb, 
                     const double Rab[], 
                     double& W, double dW[])
{
      
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha = sqrt(0.25*(sa*sa + sb*sb));
   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;
   double phi       = atan2(Rab[1],Rab[0]);

   double c,fac;

   if (((2*la+lb)%2)==1) 
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }

   if ((((std::abs(la-lb)+la+lb)/2)%2)==1)
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   W = 0.0;
   dW[0] = 0.0;
   dW[1] = 0.0;
   dW[2] = 0.0;
   double Tx,Ty,Tz;
   for (auto l=abs(la-lb);  l<=(la+lb); l+=2)
   {
         double mtmp  = 0.0;
         double mtmpx = 0.0;
         double mtmpy = 0.0;
         double mtmpz = 0.0;
         for (auto m=-l; m<=l; ++m)
         {
            double gg1 = util_gaunt(false,l,m,la,ma,lb,mb);
            util_dTesseral_lm(l,m,cos_theta,phi,Tx,Ty,Tz);
            mtmp  += gg1*util_Tesseral_lm(l,m,cos_theta,phi);
            mtmpx += gg1*Tx;
            mtmpy += gg1*Ty;
            mtmpz += gg1*Tz;
         }
         double gg2 = util_GaussBessel(la+lb,l,alpha,R);
         double gg3 = util_dGaussBessel(la+lb,l,alpha,R);
         W     += fac*(mtmp*gg2);
         dW[0] += fac*(mtmpx*gg2 + mtmp*(Rab[0]/R)*gg3);
         dW[1] += fac*(mtmpy*gg2 + mtmp*(Rab[1]/R)*gg3);
         dW[2] += fac*(mtmpz*gg2 + mtmp*(Rab[2]/R)*gg3);
         fac = -fac;
   }
   W *= c;
   dW[0] *= c;
   dW[1] *= c;
   dW[2] *= c;
}

/******************************************************
 *                                                    *
 *             util_UGaussian                         *
 *                                                    *
 ******************************************************/
/*
     Calculates the two electron one center Gaussian integral

                                            //
    WGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r) * g(lb,mb,sb;r')  
                                            || ------------------------------------  dr dr'
                                            //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
            /


   Note - this routine is equivalent to find_self_energy_coeff in the paw code
*/
double util_UGaussian(const int la, const int ma, const double sa, 
                      const int lb, const int mb, const double sb)
{
   double U = 0.0;

   if ((la==lb) && (ma==mb))
   {
      double twopi = 8.0*atan(1.0);
      double x = (2*la+1)*util_double_factorial(2*la+1);
      double y = std::pow(sqrt(0.5*(sa*sa+sb*sb)),(2*la+1));
      U = 4.0*sqrt(twopi)/(x*y);
   }

   return U;
}


//*************************** real combo versions ***********************************

/******************************************************
 *                                                    *
 *             util_WGaussian3                        *
 *                                                    *
 ******************************************************/
/*
     Calculates the 4 terms sum

        WGaussian(la,ma,sa,lb,mb,sb,Rab)
      - WGaussian(la,ma,sa,lb,mb,sm,Rab)
      - WGaussian(la,ma,sm,lb,mb,sb,Rab)
      + WGaussian(la,ma,sm,lb,mb,sm,Rab)

    of the two electron two center Gaussian integral

                                            //
    WGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                            || ------------------------------------  dr dr'
                                            //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
            /
*/
double util_WGaussian3(const int la, const int ma, const double sa, const int lb, const int mb, const double sb, const double sm, const double Rab[])
{
 
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha_ab = sqrt(0.25*(sa*sa + sb*sb));
   double alpha_am = sqrt(0.25*(sa*sa + sm*sm));
   double alpha_mb = sqrt(0.25*(sm*sm + sb*sb));
   double alpha_mm = sqrt(0.25*(sm*sm + sm*sm));

   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;

   double phi,c,fac;

   if ((std::abs(Rab[1])<1.0e-9) && (std::abs(Rab[0])<1.0e-9))
   {
      phi = 0.0;
   }
   else
   {
      phi = atan2(Rab[1],Rab[0]);
   }

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }

   if ((((std::abs(la-lb)+la+lb)/2)%2)==1)
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   double tmp = 0.0;
   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      double mtmp = 0.0;
      for (auto m=-l; m<=l; ++m)
      {
         mtmp += util_gaunt(false,l,m,la,ma,lb,mb)*util_Tesseral_lm(l,m,cos_theta,phi);
      }
      double gg2 = util_GaussBessel(la+lb,l,alpha_ab,R)
                 - util_GaussBessel(la+lb,l,alpha_am,R)
                 - util_GaussBessel(la+lb,l,alpha_mb,R)
                 + util_GaussBessel(la+lb,l,alpha_mm,R);

      tmp += fac*mtmp*gg2;
      fac = -fac;
   }

   return (c*tmp);
}

/******************************************************
 *                                                    *
 *             util_WGaussian2                        *
 *                                                    *
 ******************************************************/
/*
     Calculates the 4 term sum

        WGaussian(la,ma,sa,lb,mb,sa,Rab)
      - WGaussian(la,ma,sa,lb,mb,sm,Rab)
      - WGaussian(la,ma,sm,lb,mb,sa,Rab)
      + WGaussian(la,ma,sm,lb,mb,sm,Rab)

    of the two electron two center Gaussian integral

                                            //
    WGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                            || ------------------------------------  dr dr'
                                            //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
            /
*/
double util_WGaussian2(const int la, const int ma, const double sa, const int lb, const int mb, const double sm, const double Rab[])
{
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha_aa = sqrt(0.25*(sa*sa + sa*sa));
   double alpha_am = sqrt(0.25*(sa*sa + sm*sm));
   double alpha_mm = sqrt(0.25*(sm*sm + sm*sm));

   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;

   double phi,c,fac;

   if ((std::abs(Rab[1])<1.0e-9) && (std::abs(Rab[0])<1.0e-9))
   {
      phi = 0.0;
   }
   else
   {
      phi = atan2(Rab[1],Rab[0]);
   }

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }

   if ((((std::abs(la-lb)+la+lb)/2)%2)==1)
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   double tmp = 0.0;
   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      double mtmp = 0.0;
      for (auto m=-l; m<=l; ++m)
         mtmp += util_gaunt(false,l,m,la,ma,lb,mb) *util_Tesseral_lm(l,m,cos_theta,phi);

      double gg2 = util_GaussBessel(la+lb,l,alpha_aa,R)
                 + util_GaussBessel(la+lb,l,alpha_mm,R)
                 - 2.0*util_GaussBessel(la+lb,l,alpha_am,R);

      tmp += fac*mtmp*gg2;
      fac = -fac;
   }

   return (c*tmp);
}

/******************************************************
 *                                                    *
 *             util_dWGaussian3                       *
 *                                                    *
 ******************************************************/
/*
     Calculates the 4 terms sum

        WGaussian(la,ma,sa,lb,mb,sb,Rab)
      - WGaussian(la,ma,sa,lb,mb,sm,Rab)
      - WGaussian(la,ma,sm,lb,mb,sb,Rab)
      + WGaussian(la,ma,sm,lb,mb,sm,Rab)

     of the two electron two center Gaussian integral and it's derivative wrt to Rab

                                            //
   dWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                            || ------------------------------------  dr dr'
                                            //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Tlm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *Tlm(rhat) dr = 1
            /
*/
void  util_dWGaussian3(const int la, const int ma, const double sa, 
                       const int lb, const int mb, const double sb, 
                       const double sm, const double Rab[], double& W, double dW[])
{
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha_ab = sqrt(0.25*(sa*sa + sb*sb));
   double alpha_am = sqrt(0.25*(sa*sa + sm*sm));
   double alpha_mb = sqrt(0.25*(sm*sm + sb*sb));
   double alpha_mm = sqrt(0.25*(sm*sm + sm*sm));

   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;
   double phi       = atan2(Rab[1],Rab[0]);

   double c,fac;

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }

   if ((((std::abs(la-lb)+la+lb)/2)%2)==1)
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   W = 0.0;
   dW[0] = 0.0;
   dW[1] = 0.0;
   dW[2] = 0.0;
   double Tx,Ty,Tz;
   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      double mtmp  = 0.0;
      double mtmpx = 0.0;
      double mtmpy = 0.0;
      double mtmpz = 0.0;
      for (auto m=-l; m<=l; ++m)
      {
         double gg1 = util_gaunt(false,l,m,la,ma,lb,mb);
         util_dTesseral_lm(l,m,cos_theta,phi,Tx,Ty,Tz);
         mtmp  += gg1*util_Tesseral_lm(l,m,cos_theta,phi);
         mtmpx += gg1*Tx;
         mtmpy += gg1*Ty;
         mtmpz += gg1*Tz;
      }

      double gg2 = util_GaussBessel(la+lb,l,alpha_ab,R)
                 - util_GaussBessel(la+lb,l,alpha_am,R)
                 - util_GaussBessel(la+lb,l,alpha_mb,R)
                 + util_GaussBessel(la+lb,l,alpha_mm,R);

      double gg3 = util_dGaussBessel(la+lb,l,alpha_ab,R)
                 - util_dGaussBessel(la+lb,l,alpha_am,R)
                 - util_dGaussBessel(la+lb,l,alpha_mb,R)
                 + util_dGaussBessel(la+lb,l,alpha_mm,R);

      W     += fac*(mtmp*gg2);
      dW[0] += fac*(mtmpx*gg2 + mtmp*(Rab[0]/R)*gg3);
      dW[1] += fac*(mtmpy*gg2 + mtmp*(Rab[1]/R)*gg3);
      dW[2] += fac*(mtmpz*gg2 + mtmp*(Rab[2]/R)*gg3);

      fac = -fac;
   }
   W *= c;
   dW[0] *= c;
   dW[1] *= c;
   dW[2] *= c;
}

//*************************** real combo versions ***********************************


//*************************** complex versions ***********************************

/******************************************************
 *                                                    *
 *             util_CWGaussian                        *
 *                                                    *
 ******************************************************/
/*
     Calculates the two electron two center Gaussian integral

                                             //
    CWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                             || ------------------------------------  dr dr'
                                             //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *conjg(Ylm(rhat)) dr = 1
            /
*/
std::complex<double> util_CWGaussian(const int la, const int ma, const double sa,
                                     const int lb, const int mb, const double sb, const double Rab[])
{
   int m = ma + mb;
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha = sqrt(0.25*(sa*sa + sb*sb));
   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;

   double phi,c,fac;

   if ((std::abs(Rab[1])<1.0e-9) && (std::abs(Rab[0])<1.0e-9))
   {
      phi = 0.0;
   }
   else
   {
      phi = atan2(Rab[1],Rab[0]);
   }

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }

   if (((mb+(std::abs(la-lb)+la+lb)/2)%2)==1) 
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   std::complex<double> tmp = std::complex<double>(0.0,0.0);
   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      if (std::abs(m)<=l)
      {
         std::complex<double> mtmp = util_gaunt(true,l,m,la,ma,lb,-mb)*util_YSpherical_lm(l,m,cos_theta,phi);
         tmp +=  fac*mtmp*util_GaussBessel(la+lb,l,alpha,R);
      }
      fac = -fac;
   }

   return (c*tmp);
}


/******************************************************
 *                                                    *
 *             util_dCWGaussian                       *
 *                                                    *
 ******************************************************/
/*
     Calculates the two electron two center Gaussian integral and it's derivative wrt to Rab

                                             //
   dCWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                             || ------------------------------------  dr dr'
                                             //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l * congj(Ylm(rhat)) dr = 1
            /
*/
void util_dCWGaussian(const int la, const int ma, const double sa,
                      const int lb, const int mb, const double sb, 
                      const double Rab[],  
                      std::complex<double>& W,  std::complex<double> dW[])
{
   int m = ma + mb;
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha = sqrt(0.25*(sa*sa + sb*sb));

   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;
   double phi       = atan2(Rab[1],Rab[0]);

   double c,fac;

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }

   if (((mb+(std::abs(la-lb)+la+lb)/2)%2)==1) 
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   W = std::complex<double> (0.0,0.0);
   dW[0] = std::complex<double>(0.0,0.0);
   dW[1] = std::complex<double>(0.0,0.0);
   dW[2] = std::complex<double>(0.0,0.0);
   std::complex<double> Tx,Ty,Tz;
   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      if (std::abs(m)<=l)
      {
         double gg1 = util_gaunt(true,l,m,la,ma,lb,-mb);
         util_dYspherical_lm(l,m,cos_theta,phi,Tx,Ty,Tz);
         std::complex<double> mtmp  = gg1*util_YSpherical_lm(l,m,cos_theta,phi);
         std::complex<double> mtmpx = gg1*Tx;
         std::complex<double> mtmpy = gg1*Ty;
         std::complex<double> mtmpz = gg1*Tz;
         double gg2 = util_GaussBessel(la+lb,l,alpha,R);
         double gg3 = util_dGaussBessel(la+lb,l,alpha,R);
         W     += fac*(mtmp*gg2);
         dW[0] += fac*(mtmpx*gg2 + mtmp*(Rab[0]/R)*gg3);
         dW[1] += fac*(mtmpy*gg2 + mtmp*(Rab[1]/R)*gg3);
         dW[2] += fac*(mtmpz*gg2 + mtmp*(Rab[2]/R)*gg3);
      }
      fac = -fac;
   }
   W *= c;
   dW[0] *= c;
   dW[1] *= c;
   dW[2] *= c;
}
//*************************** complex versions ***********************************

//*************************** complex combo versions ***********************************

/******************************************************
 *                                                    *
 *             util_CWGaussian3                       *
 *                                                    *
 ******************************************************/
/*
     Calculates the two electron two center Gaussian integral

                                             //
    CWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                             || ------------------------------------  dr dr'
                                             //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *conjg(Ylm(rhat)) dr = 1
            /
*/
std::complex<double> util_CWGaussian3(const int la, const int ma, const double sa,
                                      const int lb, const int mb, const double sb,
                                      const double sm, const double Rab[])
{
   int m = ma + mb;
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha_ab = sqrt(0.25*(sa*sa + sb*sb));
   double alpha_am = sqrt(0.25*(sa*sa + sm*sm));
   double alpha_mb = sqrt(0.25*(sm*sm + sb*sb));
   double alpha_mm = sqrt(0.25*(sm*sm + sm*sm));
   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;

   double phi,c,fac;

   if ((std::abs(Rab[1])<1.0e-9) && (std::abs(Rab[0])<1.0e-9)) 
   {
      phi = 0.0;
   } 
   else
   {
      phi = atan2(Rab[1],Rab[0]);
   }

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }
      
   if (((mb+(std::abs(la-lb)+la+lb)/2)%2)==1)
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   std::complex<double> tmp = std::complex<double>(0.0,0.0);
   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      if (std::abs(m)<=l)
      {
         std::complex<double> mtmp = util_gaunt(true,l,m,la,ma,lb,-mb)*util_YSpherical_lm(l,m,cos_theta,phi);

         double gg2 = util_GaussBessel(la+lb,l,alpha_ab,R)
                    - util_GaussBessel(la+lb,l,alpha_am,R)
                    - util_GaussBessel(la+lb,l,alpha_mb,R)
                    + util_GaussBessel(la+lb,l,alpha_mm,R);

         tmp +=  fac*mtmp*gg2;
      }
      fac = -fac;
   }

   return (c*tmp);
}


/******************************************************
 *                                                    *
 *             util_CWGaussian2                       *
 *                                                    *
 ******************************************************/
/*
     Calculates the two electron two center Gaussian integral

                                             //
    CWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                             || ------------------------------------  dr dr'
                                             //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l *conjg(Ylm(rhat)) dr = 1
            /
*/
std::complex<double> util_CWGaussian2(const int la, const int ma, const double sa,
                                      const int lb, const int mb, const double sm, const double Rab[])
{
   int m = ma + mb;
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha_aa = sqrt(0.25*(sa*sa + sa*sa));
   double alpha_am = sqrt(0.25*(sa*sa + sm*sm));
   double alpha_mm = sqrt(0.25*(sm*sm + sm*sm));
   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;

   double phi,c,fac;

   if ((std::abs(Rab[1])<1.0e-9) && (std::abs(Rab[0])<1.03-9))
   {
      phi = 0.0;
   }
   else
   {
      phi = atan2(Rab[1],Rab[0]);
   }

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }

   if (((mb+(std::abs(la-lb)+la+lb)/2)%2)==1)
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   std::complex<double> tmp = std::complex<double>(0.0,0.0);
   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      if (std::abs(m)<=l)
      {
         std::complex<double> mtmp = util_gaunt(true,l,m,la,ma,lb,-mb) * util_YSpherical_lm(l,m,cos_theta,phi);

         double gg2 = util_GaussBessel(la+lb,l,alpha_aa,R)
                    + util_GaussBessel(la+lb,l,alpha_mm,R)
                    - 2.0*util_GaussBessel(la+lb,l,alpha_am,R);

         tmp += fac*mtmp*gg2;
      }
      fac = -fac;
   }

   return (c*tmp);
}


/******************************************************
 *                                                    *
 *             util_dCWGaussian3                      *
 *                                                    *
 ******************************************************/
/*
     Calculates the two electron two center Gaussian integral and it's derivative wrt to Rab

                                             //
   dCWGaussian(la,ma,sa,lb,Ra,mb,sb,Rb)   =  || g(la,ma,sa;r-Ra) * g(lb,mb,sb;r'-Rb)  
                                             || ------------------------------------  dr dr'
                                             //                |r-r'|

     where g(l,m,s; r) = C_l * |r|**l * exp(-(r/s)**2) * Ylm(rhat) 

          and C_l = 2**(l+2) / (sqrt(pi) (2*l+1)!! s**(2*l+3) )
           
     The normalization constant C_l is defined such at
            /
            | g(l,m,s;r) * |r|**l * congj(Ylm(rhat)) dr = 1
            /
*/
void util_dCWGaussian3(const int la, const int ma, const double sa,
                       const int lb, const int mb, const double sb,
                       const double sm, const double Rab[], 
                       std::complex<double>& W,  std::complex<double> dW[])
{
   int m = ma + mb;
   double pi = 4.0*atan(1.0);
   double x = util_double_factorial(2*la+1);
   double y = util_double_factorial(2*lb+1);
   double alpha_ab = sqrt(0.25*(sa*sa + sb*sb));
   double alpha_am = sqrt(0.25*(sa*sa + sm*sm));
   double alpha_mb = sqrt(0.25*(sm*sm + sb*sb));
   double alpha_mm = sqrt(0.25*(sm*sm + sm*sm));

   double R = sqrt(Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2]);
   double cos_theta = Rab[2]/R;
   double phi       = atan2(Rab[1],Rab[0]);

   double c,fac;

   if (((2*la+lb)%2)==1)
   {
      c = -32.0*pi/(x*y);
   }
   else
   {
      c = 32.0*pi/(x*y);
   }

   if (((mb+(std::abs(la-lb)+la+lb)/2)%2)==1)
   {
      fac = -1.0;
   }
   else
   {
      fac = 1.0;
   }
      
   W = std::complex<double> (0.0,0.0);
   dW[0] = std::complex<double> (0.0,0.0);
   dW[1] = std::complex<double> (0.0,0.0);
   dW[2] = std::complex<double> (0.0,0.0);

   std::complex<double> Tx,Ty,Tz;

   for (auto l=abs(la-lb); l<=(la+lb); l+=2)
   {
      if (std::abs(m)<=l)
      {
         double gg1 = util_gaunt(true,l,m,la,ma,lb,-mb);
         util_dYspherical_lm(l,m,cos_theta,phi,Tx,Ty,Tz);
         std::complex<double> mtmp  = gg1*util_YSpherical_lm(l,m,cos_theta,phi);
         std::complex<double> mtmpx = gg1*Tx;
         std::complex<double> mtmpy = gg1*Ty;
         std::complex<double> mtmpz = gg1*Tz;

         double gg2 = util_GaussBessel(la+lb,l,alpha_ab,R)
                    - util_GaussBessel(la+lb,l,alpha_am,R)
                    - util_GaussBessel(la+lb,l,alpha_mb,R)
                    + util_GaussBessel(la+lb,l,alpha_mm,R);

         double gg3 = util_dGaussBessel(la+lb,l,alpha_ab,R)
                    - util_dGaussBessel(la+lb,l,alpha_am,R)
                    - util_dGaussBessel(la+lb,l,alpha_mb,R)
                    + util_dGaussBessel(la+lb,l,alpha_mm,R);

         W     += fac*(mtmp*gg2);
         dW[0] += fac*(mtmpx*gg2 + mtmp*(Rab[0]/R)*gg3);
         dW[1] += fac*(mtmpy*gg2 + mtmp*(Rab[1]/R)*gg3);
         dW[2] += fac*(mtmpz*gg2 + mtmp*(Rab[2]/R)*gg3);
      }
      fac = -fac;
   }
   W *= c;
   dW[0] *= c;
   dW[1] *= c;
   dW[2] *= c;
}    


//*************************** complex combo versions ***********************************


}
