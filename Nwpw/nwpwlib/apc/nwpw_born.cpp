/* nwpw_born.cpp -
   Author - Eric Bylaska
*/


#include        <iostream>
#include        <cstring>
#include        <sstream>
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
                           const double q[], const double dielec, const double rcut)
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
      double gg    = std::erf(std::sqrt(dist2)/rcut);
      Gsolv -=  0.5*screen*q[ii]*q[jj]*gg/f;
   }
   return Gsolv;
}

/*******************************************
 *                                         *
 *               born_dgsolv               *
 *                                         *
 *******************************************/
static double born_dgsolv(const double screen, const double qi, const double qj, 
                          const double bi, const double bj, const double rcut, const double xx)
{
   double gg  = std::erf(std::sqrt(xx)/rcut);
   double dgg = (1.00/std::sqrt(xx*4.0*std::atan(1.0)))*std::exp(-xx/(rcut*rcut))/rcut;
   double C = std::exp(-0.25*xx/(bi*bj));
   double f = std::sqrt(xx + bi*bj*C);
   double gsolv = -0.5*screen*qi*qj*gg/f;

   return (-0.5*gsolv*(1.0-0.25*C)/(f*f) - 0.5*screen*qi*qj*dgg/f);
}

/*******************************************
 *                                         *
 *               born_fion0                *
 *                                         *
 *******************************************/

// this routine need to parallelized!!!

static void born_fion0(const int nion, const double rion[], const double bradii[],
                       const double q[], const double dielec, const double rcut,
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

      double dGsolv = born_dgsolv(screen,q[ii],q[jj],bradii[ii],bradii[jj],rcut,dist2);

      fion[3*ii]   -= 2.0*dGsolv*dx;
      fion[3*ii+1] -= 2.0*dGsolv*dy;
      fion[3*ii+2] -= 2.0*dGsolv*dz;

      fion[3*jj]   += 2.0*dGsolv*dx;
      fion[3*jj+1] += 2.0*dGsolv*dy;
      fion[3*jj+2] += 2.0*dGsolv*dz;
   }
}

/*******************************************
 *                                         *
 *               born_dVdq0                *
 *                                         *
 *******************************************/

// this routine need to parallelized!!!

static void born_dVdq0(const int nion, const double rion[], const double bradii[],
                       const double q[], const double dielec, const double rcut,
                       double u[])
{
   //double Gsolv = 0.0;
   double screen = (1.0 - 1.0/dielec);

   std::memset(u,0,nion*sizeof(double));

   for (auto ii=0; ii<nion; ++ii)
   for (auto jj=0; jj<nion; ++jj)
   {
      double dx = rion[3*ii]  -rion[3*jj];
      double dy = rion[3*ii+1]-rion[3*jj+1];
      double dz = rion[3*ii+2]-rion[3*jj+2];
      double dist2 = dx*dx + dy*dy + dz*dz;

      double gg = std::erf(std::sqrt(dist2)/rcut);
      double C  = std::exp(-0.25*dist2/(bradii[ii]*bradii[jj]));
      double f  = std::sqrt(dist2 + bradii[ii]*bradii[jj]*C);

      u[ii] += 0.5*screen*q[jj]*gg/f;
      u[jj] += 0.5*screen*q[ii]*gg/f;
      //Gsolv -=  0.5*screen*q[ii]*q[jj]/f;
   }
}


/* Constructors */

/*******************************************
 *                                         *
 *            nwpw_born::nwpw_born         *
 *                                         *
 *******************************************/
nwpw_born::nwpw_born(Ion *myionin, Parallel *myparallin, Control2& control, std::ostream& coutput)
{
   myion      = myionin;
   myparall   = myparallin;
   born_on    = control.born_on();
   born_relax = control.born_relax();

   bool oprint = ((myparall->is_master()) && (control.print_level("medium")));

   if (born_on)
   {
      dielec = control.born_dielec();
      rcut   = control.born_rcut();
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
         double  vrdii = control.born_vradii(ii);
         if (vrdii>1.0e-3)
            vradii[ii] = vrdii;
         else
         {
            int iq     = (int) std::round(myion->charge[ii]);
            vradii[ii] = vdwr[iq-1]/0.529177;
            if (vradii[ii] < 1.0e-3) 
               vradii[ii] = 1.17*1.90/0.529177;
         }
      }
      //read in bradii if in rtdb 
      for (auto ii=0; ii<myion->nion; ++ii)
      {
         double  brdii = control.born_bradii(ii);
         bradii[ii] = ((brdii>1.0e-3) ? brdii : born_radius(ii,myion->nion,myion->rion1,vradii));
      }

      /* write out Born header */
      if (oprint) 
         coutput << this->header_print();
      
   }
}

/*******************************************
 *                                         *
 *            nwpw_born::energy            *
 *                                         *
 *******************************************/
double nwpw_born::energy(const double q[])
{
   return born_energy0(myion->nion,myion->rion1,bradii,q,dielec,rcut);
}


/*******************************************
 *                                         *
 *            nwpw_born::fion              *
 *                                         *
 *******************************************/

// Needs to be parallelized!!

void nwpw_born::fion(const double q[], double fion[])
{
   born_fion0(myion->nion,myion->rion1,bradii,q,dielec,rcut,fion);
}


/*******************************************
 *                                         *
 *            nwpw_born::dVdq              *
 *                                         *
 *******************************************/
void nwpw_born::dVdq(const double q[], double u[])
{
   born_dVdq0(myion->nion,myion->rion1,bradii,q,dielec,rcut,u);
}

/*******************************************
 *                                         *
 *          nwpw_born::header_print        *
 *                                         *
 *******************************************/
std::string nwpw_born::header_print()
{
   std::stringstream stream;

   stream << std::endl;
   stream << " extended Born solvation model:" << std::endl;
   stream << "    G.D. Hawkins, C.R. Cramer, D.G. Truhlar (1995) Pairwise solute descreening of solue" << std::endl;
   stream << "     charge from a dielectric medium, Chem. Phys. Lett., vol 246, pages 122-129." << std::endl;
   stream << std::endl;
   stream << "    dielectric constant=" << std::fixed << std::setw(11) << std::setprecision(6) << dielec << std::endl;
   stream << "    dielectric rcut    =" << std::fixed << std::setw(11) << std::setprecision(6) << rcut   << std::endl;
   if (born_relax)
      stream << "    self-consistent solvation" << std::endl;
   else
      stream << "    perturbative solvation" << std::endl;
   stream << "    generalized Born radii:" << std::endl;
   for (auto ii=0; ii<myion->nion; ++ii)
      stream << "      - Born radius: " << std::left << std::setw(4) << myion->symbol(ii)
             << "- a(" << std::right << std::setw(4) << (ii+1)
             << ") = " << std::fixed << std::setw(6) << std::setprecision(3) << bradii[ii]*0.529177
             << " Angstroms (1.17*vdw radius="
             << std::fixed << std::setw(6) << std::setprecision(3) << vradii[ii]*0.529177
             << " Angstroms)" << std::endl;

   return stream.str();
}

/*******************************************
 *                                         *
 *             nwpw_born::Qprint           *
 *                                         *
 *******************************************/
std::string nwpw_born::Qprint(const double q[])
{
   double Gsolv = born_energy0(myion->nion,myion->rion1,bradii,q,dielec,rcut);

   std::stringstream stream;

   stream << std::endl << std::endl;
   stream << " generalized Born Solvation" << std::endl;
   stream << " --------------------------" << std::endl;
   stream << "   - Radii defined by G.D. Hawkins, C.R. Cramer, D.G. Truhlar (1995) Pairwise" << std::endl;
   stream << "     solute descreening of solute charge from a dielectric medium, Chem. Phys. Lett.," << std::endl;
   stream << "     vol. 246, pages 122-129." << std::endl;
   stream << std::endl;
   stream << "   - dielectric constant -eps- =" << std::fixed << std::setw(11) << std::setprecision(6) << dielec << std::endl;
   stream << "   - dielectric rcut           =" << std::fixed << std::setw(11) << std::setprecision(6) << rcut   << std::endl;
   stream << std::endl;
   for (auto ii=0; ii<myion->nion; ++ii)
      stream << "   - Born radius: " << std::left << std::setw(4) << myion->symbol(ii)
             << "- a(" << std::right << std::setw(4) << (ii+1)
             << ") = " << std::fixed << std::setw(8) << std::setprecision(3) << bradii[ii]*0.529177
             << " Angstroms (1.17*vdw radius="
             << std::fixed << std::setw(8) << std::setprecision(3) << vradii[ii]*0.529177
             << " Angstroms) - atom charge = "
             << std::fixed << std::setw(8) << std::setprecision(3) << q[ii] << std::endl;
   stream << std::endl;
   stream << "   - Gsolvation(Born) = " << std::fixed << std::setw(14) << std::setprecision(6) << Gsolv 
          << " (" << std::fixed << std::setw(8) << std::setprecision(3) << Gsolv*27.2116*23.06
          << " kcal/mol)" << std::endl;

   return stream.str();
}

/*******************************************
 *                                         *
 *        nwpw_born::shortprint            *
 *                                         *
 *******************************************/
std::string nwpw_born::shortprint(const double egas, const double esol)
{
   std::stringstream stream;

   stream << "     extended Born solvation results    " << std::endl;
   stream << "     -------------------------------    " << std::endl;
   if (egas>0.0)
   {
      stream << "     gas phase energy                 = " << std::scientific << std::setw(19) << std::setprecision(10) << egas << std::endl;
      stream << "     sol phase energy                 = " << std::scientific << std::setw(19) << std::setprecision(10) << esol << std::endl;
      stream << "     (electrostatic) solvation energy = " << std::scientific << std::setw(19) << std::setprecision(10) << egas-esol
             << "(" << std::fixed << std::setw(8) << std::setprecision(3) << (egas-esol)*27.2116*23.06 << " kcal/mol)" << std::endl;
   }
   else
   {
      stream << "     skipped: no gas phase energy " << std::endl;
      stream << "     sol phase energy                 = " << std::scientific << std::setw(19) << std::setprecision(10) << esol << std::endl;
   }
   stream << std::endl << std::endl;

   return stream.str();
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

      rtdbjson["nwpw"]["born"]["on"]     = born_on;
      rtdbjson["nwpw"]["born"]["relax"]  = born_relax;
      rtdbjson["nwpw"]["born"]["bradii"] = std::vector<double>(bradii, &bradii[myion->nion]);
      rtdbjson["nwpw"]["born"]["vradii"] = std::vector<double>(vradii, &vradii[myion->nion]);

      rtdbstring = rtdbjson.dump();
   }
}


}
