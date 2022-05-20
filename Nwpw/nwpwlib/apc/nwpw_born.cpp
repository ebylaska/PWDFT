/* nwpw_born.cpp -
   Author - Eric Bylaska
*/


#include        <iostream>
#include        <cstring>
#include        <cmath>

#include        "nwpw_timing.hpp"

#include        "blas.h"

#include        "nwpw_born.hpp"

#include        "json.hpp"
using json = nlohmann::json;

namespace pwdft {


/*******************************************
 *                                         *
 *               born_radius               *
 *                                         *
 *******************************************/
static double born_radius(const int ii, const int nion, const double rion[], const double vradii[])
{
   double bornr1 = 1.0/vradii[ii];

   for (auto jj=0; jj<nion; ++jj)
   {
      double dx = rion[3*ii]-rion[3*jj];
      double dy = rion[3*ii+1]-rion[3*jj+1];
      double dz = rion[3*ii+2]-rion[3*jj+2];
      double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

      if (dist>0.1)
      {
         double L=0.0;
         if ((dist+vradii[jj])<=vradii[ii]) 
         {
            L = 1.0;
         }
         else if (((dist-vradii[jj])<=vradii[ii]) && (vradii[ii]<=(dist+vradii[jj])))
         {
            L = vradii[jj];
         }
         else if (vradii[ii]<=(dist+vradii[jj]))
         {
            L = dist - vradii[jj];
         }

         double U=0.0;
         if ((dist+vradii[jj])<=vradii[ii])
         {
            U = 1.0;
         }
         else if (vradii[ii]<(dist+vradii[jj]))
         {
            U = dist + vradii[jj];
         }

         if ((U>0.0) && (L>0.0))
         {
            double invL = 1.0/L;
            double invU = 1.0/U;
            double delta = -0.5*( (invL - invU)
                                 + 0.25*dist*(invU*invU - invL*invL)
                                 + 0.50/dist*std::log(L/U)
                                 + 0.25*vradii[jj]*vradii[jj]
                                  /dist * (invL*invL - invU*invU) );
            bornr1 = bornr1 + delta;
         }
      }
   }

   if (bornr1<1.0e-6) bornr1 = 1.0e-6;
   return (1.0/bornr1);
}

/*******************************************
 *                                         *
 *               born_energy0              *
 *                                         *
 *******************************************/

// this routine need to parallelized!!!

static double born_energy0(const int nion, const double rion[], const double bradii[],
                           const double q[], const double dielec)
{
   double Gsolv = 0.0;
   double screen = (1.0 - 1.0/dielec);

   for (auto ii=0; ii<nion; ++ii)
   for (auto jj=0; jj<nion; ++jj)
   {
      double dx = rion[3*ii]  -rion[3*jj];
      double dy = rion[3*ii+1]-rion[3*jj+1];
      double dz = rion[3*ii+2]-rion[3*jj+2];
      double dist2 = dx*dx + dy*dy + dz*dz;
      double C     = std::exp(-0.25*dist2/(bradii[ii]*bradii[jj]));
      double f     = std::sqrt(dist2 + bradii[ii]*bradii[jj]*C);
      Gsolv -=  0.5*screen*q[ii]*q[jj]/f;
   }
   return Gsolv;
}

/*******************************************
 *                                         *
 *               born_dgsolv               *
 *                                         *
 *******************************************/
static double born_dgsolv(const double screen, const double qi, const double qj, 
                        const double bi, const double bj, const double xx)
{
   double C = std::exp(-0.25*xx/(bi*bj));
   double f = std::sqrt(xx + bi*bj*C);
   double gsolv = -0.5*screen*qi*qj/f;

   return (-0.5*gsolv*(1.0-0.25*C)/(f*f));
}

/*******************************************
 *                                         *
 *               born_fion0                *
 *                                         *
 *******************************************/

// this routine need to parallelized!!!

static void born_fion0(const int nion, const double rion[], const double bradii[],
                       const double q[], const double dielec,
                       double fion[])
{
   double screen = (1.0 - 1.0/dielec);

   for (auto ii=0; ii<nion; ++ii)
   for (auto jj=0; jj<nion; ++jj)
   {
      double dx = rion[3*ii]  -rion[3*jj];
      double dy = rion[3*ii+1]-rion[3*jj+1];
      double dz = rion[3*ii+2]-rion[3*jj+2];
      double dist2 = dx*dx + dy*dy + dz*dz;

      double dGsolv = born_dgsolv(screen,q[ii],q[jj],bradii[ii],bradii[jj],dist2);

      fion[3*ii]   -= 2.0*dGsolv*dx;
      fion[3*ii+1] -= 2.0*dGsolv*dy;
      fion[3*ii+2] -= 2.0*dGsolv*dz;

      fion[3*jj]   += 2.0*dGsolv*dx;
      fion[3*jj+1] += 2.0*dGsolv*dy;
      fion[3*jj+2] += 2.0*dGsolv*dz;
   }
}




/* Constructors */

/*******************************************
 *                                         *
 *            nwpw_born::nwpw_born         *
 *                                         *
 *******************************************/
nwpw_born::nwpw_born(Ion *myionin, Control2& control)
{
   myion      = myionin;
   born_on    = control.born_on();
   born_relax = control.born_relax();

   if (born_on)
   {
      dielec = control.born_dielec();
      vradii = new (std::nothrow) double [myion->nion];
      bradii = new (std::nothrow) double [myion->nion];

      // radii for heavy elements: 1.17*1.9
      double vdwr[102] = {1.300,1.638,1.404,1.053,2.0475,2.00,
                          1.830,1.720,1.720,1.8018,1.755,1.638,
                          1.404,2.457,2.106,2.160,2.05,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.160,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.320,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223,
                          2.223,2.223,2.223,2.223,2.223,2.223};

        for (auto ii=0; ii<myion->nion; ++ii)
        {
           int iq     = (int) std::round(myion->charge[ii]);
           vradii[ii] = vdwr[iq-1]/0.529177;
           if (vradii[ii] < 1.0e-3) 
              vradii[ii] = 1.17*1.90/0.529177;

           //read in if rtdb 
           bradii[ii] = born_radius(ii,myion->nion,myion->rion1, vradii);
        }
   }
}

/*******************************************
 *                                         *
 *            nwpw_born::energy            *
 *                                         *
 *******************************************/
double nwpw_born::energy(const double q[])
{
}


/*******************************************
 *                                         *
 *            nwpw_born::fion              *
 *                                         *
 *******************************************/
void nwpw_born::fion(const double q[], double fion[])
{
}


/*******************************************
 *                                         *
 *            nwpw_born::dVdq              *
 *                                         *
 *******************************************/
void nwpw_born::dVdq(const double q[], double u[])
{
}

/*******************************************
 *                                         *
 *          nwpw_born::header_print        *
 *                                         *
 *******************************************/
std::string nwpw_born::header_print()
{
}

/*******************************************
 *                                         *
 *             nwpw_born::Qprint           *
 *                                         *
 *******************************************/
std::string nwpw_born::Qprint(const double q[])
{
}

/*******************************************
 *                                         *
 *           nwpw_born::final_print        *
 *                                         *
 *******************************************/
std::string nwpw_born::final_print(const double egas, const double esolv)
{
}


/********************************************
 *                                          *
 *          nwpw_born::writejsonstr         *
 *                                          *
 ********************************************/
void nwpw_born::writejsonstr(std::string& rtdbstring)
{
   if (born_on)
   {
      /* save restart information into json */
      auto rtdbjson = json::parse(rtdbstring);
      auto bornjson = rtdbjson["nwpw"]["Born"];

      bornjson["born_on"]    = born_on;
      bornjson["born_relax"] = born_relax;
      bornjson["bradii"]     = std::vector<double>(bradii, &bradii[myion->nion]);

      rtdbstring = rtdbjson.dump();
   }
}



}
