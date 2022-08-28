
#include        <iomanip>
#include        <cstring>
#include        <cmath>
#include        <iostream>
#include        <fstream>
#include        <cstdio>
#include        <string>
#include        <unistd.h>
#include        <sys/stat.h>
#include        <string>
#include        <vector>
#include        <set>
#include        <map>
#include        <algorithm>

#include "parsestring.hpp"

#include	"qmmm.hpp"


#define MASTER          0
#define ANGTOBOHR       1.88972687777


/**************************************************
 *                                                *
 *                symboltomass                    *
 *                                                *
 **************************************************/
static double symboltomass(std::string symbol)
{
   double mass = 1.008;
   if (mystring_contains(symbol,"H")) mass = 1.008;
   if (mystring_contains(symbol,"B")) mass = 11.00931;
   if (mystring_contains(symbol,"C")) mass = 12.0;
   if (mystring_contains(symbol,"N")) mass = 14.00307;
   if (mystring_contains(symbol,"O")) mass = 15.99491;
   if (mystring_contains(symbol,"F")) mass = 18.9984;
   if (mystring_contains(symbol,"P")) mass = 30.97376;
   if (mystring_contains(symbol,"S")) mass = 31.97207;
   if (mystring_contains(symbol,"K")) mass = 38.96371;

   if (mystring_contains(symbol,"He")) mass = 4.0026;
   if (mystring_contains(symbol,"Li")) mass = 7.016;
   if (mystring_contains(symbol,"Be")) mass = 9.01218;
   if (mystring_contains(symbol,"Ne")) mass = 19.99244;

   if (mystring_contains(symbol,"Na")) mass = 22.9898;
   if (mystring_contains(symbol,"Mg")) mass = 23.98504;
   if (mystring_contains(symbol,"Al")) mass = 26.98154;
   if (mystring_contains(symbol,"Si")) mass = 27.97693;
   if (mystring_contains(symbol,"Cl")) mass = 34.96885;
   return mass;
}


/**************************************************
 *                                                *
 *             symboltoepsilon                    *
 *                                                *
 **************************************************/
static double symboltoepsilon(std::string symbol)
{
   double epsilon = 0.044;
   if (mystring_contains(symbol,"H")) epsilon = 0.044;
   if (mystring_contains(symbol,"B")) epsilon = 0.180;
   if (mystring_contains(symbol,"C")) epsilon = 0.105;
   if (mystring_contains(symbol,"N")) epsilon = 0.069;
   if (mystring_contains(symbol,"O")) epsilon = 0.060;
   if (mystring_contains(symbol,"F")) epsilon = 0.050;
   if (mystring_contains(symbol,"P")) epsilon = 0.305;
   if (mystring_contains(symbol,"S")) epsilon = 0.274;
   if (mystring_contains(symbol,"K")) epsilon = 0.035;

   if (mystring_contains(symbol,"He")) epsilon = 0.056;
   if (mystring_contains(symbol,"Li")) epsilon = 0.025;
   if (mystring_contains(symbol,"Be")) epsilon = 0.085;
   if (mystring_contains(symbol,"Ne")) epsilon = 0.042;

   if (mystring_contains(symbol,"Na")) epsilon = 0.030;
   if (mystring_contains(symbol,"Mg")) epsilon = 0.111;
   if (mystring_contains(symbol,"Al")) epsilon = 0.505;
   if (mystring_contains(symbol,"Si")) epsilon = 0.402;
   if (mystring_contains(symbol,"Cl")) epsilon = 0.227;
   return epsilon;
}


/**************************************************
 *                                                *
 *              symboltosigma                     *
 *                                                *
 **************************************************/
static double symboltosigma(std::string symbol)
{
   double sigma = 2.886;
   if (mystring_contains(symbol,"H")) sigma = 2.886;
   if (mystring_contains(symbol,"B")) sigma = 4.083;
   if (mystring_contains(symbol,"C")) sigma = 3.851;
   if (mystring_contains(symbol,"N")) sigma = 3.66;
   if (mystring_contains(symbol,"O")) sigma = 3.50;
   if (mystring_contains(symbol,"F")) sigma = 3.364;
   if (mystring_contains(symbol,"P")) sigma = 4.147;
   if (mystring_contains(symbol,"S")) sigma = 4.035;
   if (mystring_contains(symbol,"K")) sigma = 3.812;

   if (mystring_contains(symbol,"He")) sigma = 2.362;
   if (mystring_contains(symbol,"Li")) sigma = 2.451;
   if (mystring_contains(symbol,"Be")) sigma = 2.745;
   if (mystring_contains(symbol,"Ne")) sigma = 3.243;

   if (mystring_contains(symbol,"Na")) sigma = 2.983;
   if (mystring_contains(symbol,"Mg")) sigma = 3.021;
   if (mystring_contains(symbol,"Al")) sigma = 4.499;
   if (mystring_contains(symbol,"Si")) sigma = 4.295;
   if (mystring_contains(symbol,"Cl")) sigma = 3.947;
   return sigma;
}






/**************************************************
 *                                                *
 *                spring_bond_frag                *
 *                                                *
 **************************************************/
static double spring_bond_frag(const int nbond, const int indx[], const double Kr0[], const double rion[])
{
   double E = 0.0;
   for (auto b=0; b<nbond; ++b)
   {
      int ii = indx[0+2*b];
      int jj = indx[1+2*b];

      double x = rion[3*ii]   - rion[3*jj];
      double y = rion[3*ii+1] - rion[3*jj+1];
      double z = rion[3*ii+2] - rion[3*jj+2];
      double r = std::sqrt(x*x + y*y + z*z);
      double dr = r - Kr0[1+2*b];
      E += Kr0[2*b]*dr*dr;

   }
   return E;
}

/**************************************************
 *                                                *
 *                spring_bond_force_frag          *
 *                                                *
 **************************************************/
static void spring_bond_force_frag(const int nbond, const int indx[], const double Kr0[], const double rion[], double fion[])
{
   for (auto b=0; b<nbond; ++b)
   {

      int ii = indx[0+2*b];
      int jj = indx[1+2*b];

      double x = rion[3*ii]   - rion[3*jj];
      double y = rion[3*ii+1] - rion[3*jj+1];
      double z = rion[3*ii+2] - rion[3*jj+2];
      double r = std::sqrt(x*x + y*y + z*z);
      double dr = r - Kr0[1+2*b];
      double dE = 2.0*Kr0[2*b]*dr/r;

      fion[3*ii]   -= x*dE;
      fion[3*ii+1] -= y*dE;
      fion[3*ii+2] -= z*dE;

      fion[3*jj]   += x*dE;
      fion[3*jj+1] += y*dE;
      fion[3*jj+2] += z*dE;
   }
}

/**************************************************
 *                                                *
 *               spring_angle_frag                *
 *                                                *
 **************************************************/
static double spring_angle_frag(const int nangle, const int indx[], const double Kr0[], const double rion[])
{
   double E = 0.0;
   for (auto a=0; a<nangle; ++a)
   {
      int ii = indx[0+3*a];
      int jj = indx[1+3*a];
      int kk = indx[2+3*a];

      double x1 = rion[3*ii]   - rion[3*jj];
      double y1 = rion[3*ii+1] - rion[3*jj+1];
      double z1 = rion[3*ii+2] - rion[3*jj+2];
      double r1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);

      double x2 = rion[3*kk]   - rion[3*jj];
      double y2 = rion[3*kk+1] - rion[3*jj+1];
      double z2 = rion[3*kk+2] - rion[3*jj+2];
      double r2 = std::sqrt(x2*x2 + y2*y2 + z2*z2);

      double denom = r1*r2;
      if (denom>1.0e-11)
      {
         double ctheta = (x1*x2+y1*y2+z1*z2)/(denom);
         if (ctheta > 1.0)  ctheta = 1.0;
         if (ctheta < -1.0) ctheta = -1.0;
         double q = std::acos(ctheta)-Kr0[1+2*a];
         E += Kr0[2*a]*q*q;
      }
   }
   return E;
}

/**************************************************
 *                                                *
 *              spring_angle_force_frag           *
 *                                                *
 **************************************************/
static void spring_angle_force_frag(const int nangle, const int indx[], const double Kr0[], const double rion[], double fion[])
{
   for (auto a=0; a<nangle; ++a)
   {
      int ii = indx[0+3*a];
      int jj = indx[1+3*a];
      int kk = indx[2+3*a];

      double x1 = rion[3*ii]   - rion[3*jj];
      double y1 = rion[3*ii+1] - rion[3*jj+1];
      double z1 = rion[3*ii+2] - rion[3*jj+2];
      double r1sq = x1*x1 + y1*y1 + z1*z1;
      double r1 = std::sqrt(r1sq);

      double x2 = rion[3*kk]   - rion[3*jj];
      double y2 = rion[3*kk+1] - rion[3*jj+1];
      double z2 = rion[3*kk+2] - rion[3*jj+2];
      double r2sq = x2*x2 + y2*y2 + z2*z2;
      double r2 = std::sqrt(r2sq);

      double denom = r1*r2;
      if (denom>1.0e-11)
      {
         double ctheta = (x1*x2+y1*y2+z1*z2)/(denom);
         if (ctheta > 1.0)  ctheta = 1.0;
         if (ctheta < -1.0) ctheta = -1.0;
         double stheta = std::sqrt(1.0 - ctheta*ctheta);
         if (stheta < 0.001) stheta = 0.001;
         stheta = 1.0/stheta;

         double q  = std::acos(ctheta)-Kr0[1+2*a];
         double tk = Kr0[2*a]*q;

         double  aa     = 2.0*tk*stheta;
         double  a11    =  aa*ctheta/r1sq;
         double  a12    = -aa/(denom);
         double  a22    =  aa*ctheta/r2sq;

         double  vx1 = a11*x1 + a12*x2;
         double  vx2 = a22*x2 + a12*x1;

         double  vy1 = a11*y1 + a12*y2;
         double  vy2 = a22*y2 + a12*y1;

         double  vz1 = a11*z1 + a12*z2;
         double  vz2 = a22*z2 + a12*z1;

         fion[3*ii]   -= vx1;
         fion[3*ii+1] -= vy1;
         fion[3*ii+2] -= vz1;

         fion[3*jj]   += vx1 + vx2;
         fion[3*jj+1] += vy1 + vy2;
         fion[3*jj+2] += vz1 + vz2;

         fion[3*kk]   -= vx2;
         fion[3*kk+1] -= vy2;
         fion[3*kk+2] -= vz2;
      }
   }
}


/**************************************************
 *                                                *
 *                   LJ_energy                    *
 *                                                *
 **************************************************/
static double LJ_energy(const double epsilon12, const double sigma12,
                        const double r1[], const double r2[])
{

   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double r = std::sqrt(x*x + y*y + z*z);
   double u = (sigma12/r);
   double u6  = u*u*u*u*u*u;
   double u12 = u6*u6;

   return (4.0*epsilon12*(u12-u6));
}

/**************************************************
 *                                                *
 *                   LJ_force                     *
 *                                                *
 **************************************************/
static void LJ_force(const double epsilon12, const double sigma12,
                     const double r1[], double f1[],
                     const double r2[], double f2[])
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double r = std::sqrt(x*x + y*y + z*z);
   double u = (sigma12/r);
   double u6  = u*u*u*u*u*u;
   double u12 = u6*u6;

   double dVLJ = -(4.00*epsilon12/r)*(12.0*u12-6.0*u6);

   f1[0] += (x/r)*dVLJ;
   f1[1] += (y/r)*dVLJ;
   f1[2] += (z/r)*dVLJ;

   f2[0] -= (x/r)*dVLJ;
   f2[1] -= (y/r)*dVLJ;
   f2[2] -= (z/r)*dVLJ;
}


/**************************************************
 *                                                *
 *               Q_Switching                      *
 *                                                *
 **************************************************/
static double Q_Switching(const double Rin, const double Rout, const double r)
{
   double s;

   if (r<=Rin)  return 0.0;
   if (r>=Rout) return 1.0;

   double c1 = (Rout-Rin);
   double c3 = c1*c1*c1;
   double c2 = 3.0*Rout-Rin;
   return(((r-Rin)*(r-Rin)) * (c2-2.0*r)/c3);
}

/**************************************************
 *                                                *
 *               Q_dSwitching                     *
 *                                                *
 **************************************************/
static void Q_dSwitching(const double Rin, const double Rout, const double r,
                         double *S, double *dS)
{
   if (r<=Rin)  { *S = 0.0; *dS = 0.0; return; }
   if (r>=Rout) { *S = 1.0; *dS = 0.0; return; }

   double c1 = (Rout-Rin);
   double c3 = c1*c1*c1;
   double c2 = 3.0*Rout-Rin;

   *S  = ((r-Rin)*(r-Rin)) * (c2-2.0*r)/c3;
   *dS = 2.0*(r-Rin)*((c2-2.0*r)/c3) - 2.0*((r-Rin)*(r-Rin)) /c3;
   return;
}

/**************************************************
 *                                                *
 *            Q_Electrostatic_potential           *
 *                                                *
 **************************************************/
static double Q_Electrostatic_potential(const double r1[], const double q1,
                                 const double r2[])
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double r = std::sqrt(x*x + y*y + z*z);

   return (q1/r);
}

/**************************************************
 *                                                *
 *            Q_Electrostatic_self                *
 *                                                *
 **************************************************/
static double Q_Electrostatic_self(const double r1[], const double q1,
                            const double r2[], const double q2)
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double r = std::sqrt(x*x + y*y + z*z);

   return (q1*q2/r);
}

/**************************************************
 *                                                *
 *          Q_Electrostatic_Force_self            *
 *                                                *
 **************************************************/

static void Q_Electrostatic_Force_self(const double r1[], const double q1, double f1[],
                                       const double r2[], const double q2, double f2[])
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double rr = (x*x + y*y + z*z);
   double r  = std::sqrt(rr);

   double der = -q1*q2/rr;
   double fx = -(x/r)*der;
   double fy = -(y/r)*der;
   double fz = -(z/r)*der;

   f2[0] += fx;
   f2[1] += fy;
   f2[2] += fz;

   f1[0] -= fx;
   f1[1] -= fy;
   f1[2] -= fz;
}


/**************************************************
 *                                                *
 *          Q_Electrostatic_SmForce               *
 *                                                *
 **************************************************/

// Calculates the Electrostatic force between two centers

static void Q_Electrostatic_SmForce(const double Sm, 
                                    const double r1[], const double q1, double f1[],
                                    const double r2[], const double q2, double f2[])
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double rr = (x*x + y*y + z*z);
   double r  = std::sqrt(rr);

   double der = -Sm*q1*q2/rr;
   double fx = -(x/r)*der;
   double fy = -(y/r)*der;
   double fz = -(z/r)*der;

   f2[0] += fx;
   f2[1] += fy;
   f2[2] += fz;

   f1[0] -= fx;
   f1[1] -= fy;
   f1[2] -= fz;
}





/**************************************************
 *                                                *
 *         Q_Electrostatic_mForce_self            *
 *                                                *
 **************************************************/

static void Q_Electrostatic_mForce_self(const double r1[], const double q1, double f1[],
                                        const double r2[], const double q2, double f2[])
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double rr = (x*x + y*y + z*z);
   double r  = std::sqrt(rr);

   double der = -q1*q2/rr;
   double fx = -(x/r)*der;
   double fy = -(y/r)*der;
   double fz = -(z/r)*der;

   f2[0] -= fx;
   f2[1] -= fy;
   f2[2] -= fz;

   f1[0] += fx;
   f1[1] += fy;
   f1[2] += fz;
}


/****************************************************
 *                                                  *
 *                    Q_cm                          *
 *                                                  *
 ****************************************************/
static void Q_cm(const int n, const int ks, const double amass[], const double rion[],
                 double rcm[])
{
   rcm[0] = 0.0; rcm[1] = 0.0; rcm[2] = 0.0;
   double m = 0.0;
   int kk = ks;
   for (auto k=0; k<n; ++k)
   {
      rcm[0] += amass[kk]*rion[3*kk];
      rcm[1] += amass[kk]*rion[3*kk+1];
      rcm[2] += amass[kk]*rion[3*kk+1];
      m += amass[kk];
      ++kk;
   }
   rcm[0] /= m;
   rcm[1] /= m;
   rcm[2] /= m;
}

/****************************************************
 *                                                  *
 *                    Q_cmm                         *
 *                                                  *
 ****************************************************/
static double Q_cmm(const int n, const int ks, const double amass[], const double rion[],
                    double rcm[])
{
   rcm[0] = 0.0; rcm[1] = 0.0; rcm[2] = 0.0;
   double m = 0.0;
   int kk = ks;
   for (auto k=0; k<n; ++k)
   {
      rcm[0] += amass[kk]*rion[3*kk];
      rcm[1] += amass[kk]*rion[3*kk+1];
      rcm[2] += amass[kk]*rion[3*kk+1];
      m += amass[kk];
      ++kk;
   }
   rcm[0] /= m;
   rcm[1] /= m;
   rcm[2] /= m;

   return m;
}





/****************************************************
 *                                                  *
 *                  Q_E_frag_frag                   *
 *                                                  *
 ****************************************************/

// This function takes out frag-frag Coulomb energies if rcm<Rin,
// takes out part of frag-frag Coulomb energies if Rcin<rcm<Rcout,
// or does not subtract frag-frag Coulomb energy if rcm>=Rcout

static double Q_E_frag_frag(const int n1, const int ks1, const double rw1_cm[],
                            const int n2, const int ks2, const double rw2_cm[],
                            const double Rin, const double Rout,
                            const double qion[], const double rion[])
{
   double E = 0.0;
   double x = rw1_cm[0] - rw2_cm[0];
   double y = rw1_cm[1] - rw2_cm[1];
   double z = rw1_cm[2] - rw2_cm[2];
   double rcm = std::sqrt(x*x + y*y + z*z);

   if (rcm < Rout)
   {
      double Sm  = Q_Switching(Rin,Rout,rcm)-1.0;
      for (auto a=0; a<n1; ++a)
      for (auto b=0; b<n2; ++b)
         E += Sm*Q_Electrostatic_self(&rion[3*(ks1+a)],qion[(ks1+a)],
                                      &rion[3*(ks2+b)],qion[(ks2+b)]);
   }

   return E;
}



/****************************************************
 *                                                  *
 *                  Q_fion_frag_frag                *
 *                                                  *
 ****************************************************/
static void Q_fion_frag_frag(const int n1, const int ks1, const double rw1_cm[], const double m1,
                             const int n2, const int ks2, const double rw2_cm[], const double m2,
                             const double Rin, const double Rout,
                             const double qion[], const double mass[], const double rion[], double fion[])
{
   double x = rw1_cm[0] - rw2_cm[0];
   double y = rw1_cm[1] - rw2_cm[1];
   double z = rw1_cm[2] - rw2_cm[2];
   double rcm = std::sqrt(x*x + y*y + z*z);

   if (rcm < Rout)
   {
      double E=0.0;
      double Sm,dSm;
      Q_dSwitching(Rin,Rout,rcm,&Sm,&dSm);
      Sm = Sm - 1.0;

      // calculate E, and -Sm*grad(E)
      for (auto a=0; a<n1; ++a)
      for (auto b=0; b<n2; ++b)
      {
         E += Q_Electrostatic_self(&rion[3*(ks1+a)],qion[(ks1+a)],
                                   &rion[3*(ks2+b)],qion[(ks2+b)]);
         Q_Electrostatic_SmForce(Sm,&rion[3*(ks1+a)],qion[(ks1+a)],&fion[3*(ks1+a)],
                                    &rion[3*(ks2+b)],qion[(ks2+b)],&fion[3*(ks2+b)]);
      }

      //calculate -E*grad(Sm) 
      for (auto a=0; a<n1; ++a)
      {
         fion[3*(ks1+a)]   -= E*dSm*(x/rcm)*mass[(ks1+a)]/m1;
         fion[3*(ks1+a)+1] -= E*dSm*(y/rcm)*mass[(ks1+a)]/m1;
         fion[3*(ks1+a)+2] -= E*dSm*(z/rcm)*mass[(ks1+a)]/m1;
      }
      for (auto b=0; b<n2; ++b)
      {
         fion[3*(ks2+b)]   += E*dSm*(x/rcm)*mass[(ks2+b)]/m2;
         fion[3*(ks2+b)+1] += E*dSm*(y/rcm)*mass[(ks2+b)]/m2;
         fion[3*(ks2+b)+2] += E*dSm*(z/rcm)*mass[(ks2+b)]/m2;
      }
   }
}







/*******************************************
 *                                         *
 *        QMMM_Operator::QMMM_Operator     *
 *                                         *
 *******************************************/
QMMM_Operator::QMMM_Operator(std::string nwinput)
{



   // set up geometry sizes and indexes: nion, nkatm, katm[nion] 
   std::string geomblock;
   if (mystring_contains(mystring_lowercase(nwinput),"geometry"))
   {
      geomblock  = mystring_rtrim(mystring_ltrim(mystring_split(mystring_split(nwinput,"nocenter")[1],"end")[0]));
   }
   std::vector<std::string> geomlines = mystring_split(geomblock,"\n");
   nion  = geomlines.size();
   nkatm = 0;
   katm   = new (std::nothrow) int [nion];
   symbol = new (std::nothrow) std::string [nion];
   rion   = new (std::nothrow) double [3*nion];
   qion   = new (std::nothrow) double [nion];
   mass   = new (std::nothrow) double [nion];
   epsilon   = new (std::nothrow) double [nion];
   sigma     = new (std::nothrow) double [nion];
   {  std::set<std::string> anameset;
      for (auto & line: geomlines)
      {
        std::vector<std::string> ss = mystring_split0(line);
        anameset.emplace(ss[0]);
      }
      nkatm = anameset.size();
   }
   aname = new (std::nothrow) std::string [nkatm];

   {  int ia=0;
      int ii=0;
      int n = sizeof(aname)/sizeof(aname[0]);
      std::set<std::string> anameset;
      for (auto & line: geomlines)
      {
         std::vector<std::string> ss = mystring_split0(line);
         if (anameset.count(ss[0])==0)
         {
            anameset.emplace(ss[0]);
            aname[ia] = ss[0];
            ++ia;
         }
         auto itr = find(aname, aname + n, ss[0]);
         katm[ii] = std::distance(aname,itr);
         ++ii;
      }
   } 

   std::memset(qion,0,   nion*sizeof(double));




   //set up qmmm fragments and sizes
   std::vector<std::string> fragments;
   if  (mystring_contains(mystring_lowercase(nwinput),"fragment"))
   {
      fragments = mystring_split(nwinput,"fragment");
      fragments.erase(fragments.begin());
      nkfrag = fragments.size();
      for (auto & frag: fragments)
      {
         int nfrag0=0;
         frag = mystring_rtrim(mystring_ltrim(mystring_split(frag,"end")[0]));

         if (mystring_contains(mystring_lowercase(frag),"index_start"))
            nfrag0 = mystring_split0(mystring_trim(mystring_split(mystring_split(frag,"index_start")[1],"\n")[0])).size();

         if (mystring_contains(mystring_lowercase(frag),"bond_spring"))
            nbond += (mystring_split(frag,"bond_spring").size() - 1);

         if (mystring_contains(mystring_lowercase(frag),"angle_spring"))
            nangle += (mystring_split(frag,"angle_spring").size() - 1);

         nfrag += nfrag0;
      }
   }

   // set QM and MM LJ parameters
   if  (mystring_contains(mystring_lowercase(nwinput),"lj_qm_parameters"))
   {
      std::vector<std::string> ss;
      std::vector<std::string> lj_fragments = mystring_split(nwinput,"lj_qm_parameters");
      lj_fragments.erase(lj_fragments.begin());
      for (auto & lj_frag: lj_fragments)
      {
         ss = mystring_split0(lj_frag);
         lj_qm_data[ss[0]].push_back(std::stod(ss[1]));
         lj_qm_data[ss[0]].push_back(std::stod(ss[2]));
      }
   }
   if  (mystring_contains(mystring_lowercase(nwinput),"lj_mm_parameters"))
   {
      std::vector<std::string> ss;
      std::vector<std::string> lj_fragments = mystring_split(nwinput,"lj_mm_parameters");
      lj_fragments.erase(lj_fragments.begin());
      for (auto & lj_frag: lj_fragments)
      {
         ss = mystring_split0(lj_frag);
         lj_mm_data[ss[0]].push_back(std::stod(ss[1]));
         lj_mm_data[ss[0]].push_back(std::stod(ss[2]));
      }
   }

   // get the qm atoms

   nion_qm = 0;
   for (auto ii=0; ii<nion; ++ii)
   {
      if (!mystring_contains(mystring_lowercase(geomlines[ii]),"#"))
      {
         std::vector<std::string> ss = mystring_split0(geomlines[ii]);
         symbol[nion_qm]  = ss[0];
         if (lj_qm_data.count(symbol[nion_qm]) > 0)
         {
            epsilon[nion_qm] = lj_qm_data[symbol[nion_qm]][0]/23.06/27.2116;
            sigma[nion_qm]   = lj_qm_data[symbol[nion_qm]][1]*ANGTOBOHR;
         }
         else
         {
            epsilon[nion_qm] = symboltoepsilon(symbol[nion_qm])/23.06/27.2116;
            sigma[nion_qm]   = symboltosigma(symbol[nion_qm])*ANGTOBOHR;
         }
         mass[nion_qm]    = symboltomass(symbol[nion_qm])*1822.89;
         rion[3*nion_qm]   = std::stod(ss[1])*ANGTOBOHR;
         rion[3*nion_qm+1] = std::stod(ss[2])*ANGTOBOHR;
         rion[3*nion_qm+2] = std::stod(ss[3])*ANGTOBOHR;
         ++nion_qm;
      }
   }

   // get the mm atoms
   nion_mm = 0;
   for (auto ii=0; ii<nion; ++ii)
   {
      if (mystring_contains(mystring_lowercase(geomlines[ii]),"#"))
      {
         std::vector<std::string> ss = mystring_split0(geomlines[ii]);
         std::string str = ss[0];
         str.erase(std::remove(str.begin(),str.end(),'#'),str.end());
         symbol[(nion_qm+nion_mm)]  = str;
         if (lj_mm_data.count(str) > 0)
         {
            epsilon[(nion_qm+nion_mm)] = lj_mm_data[symbol[(nion_qm+nion_mm)]][0]/23.06/27.2116;
            sigma[(nion_qm+nion_mm)]   = lj_mm_data[symbol[(nion_qm+nion_mm)]][1]*ANGTOBOHR;
         }
         else
         {
            epsilon[(nion_qm+nion_mm)] = 0.0;
            sigma[(nion_qm+nion_mm)]   = 1.0;
         }
         mass[(nion_qm+nion_mm)]      = symboltomass(symbol[(nion_qm+nion_mm)])*1822.89;
         rion[3*(nion_qm+nion_mm)]   = std::stod(ss[1])*ANGTOBOHR;
         rion[3*(nion_qm+nion_mm)+1] = std::stod(ss[2])*ANGTOBOHR;
         rion[3*(nion_qm+nion_mm)+2] = std::stod(ss[3])*ANGTOBOHR;
         qion[(nion_qm+nion_mm)]      = std::stod(ss[5]);
         ++nion_mm;
      }
   }



   // switching parameters
   switch_Rin  = new (std::nothrow) double [nkfrag]; //  = { (2.0160/0.529177) };
   switch_Rout = new (std::nothrow) double [nkfrag]; //  = { (3.1287/0.529177) };

   // Allocate QMMM indexing arrays 
   size_frag        = new (std::nothrow) int [nkfrag];
   nfrag_frag       = new (std::nothrow) int [nkfrag];
   nbond_frag       = new (std::nothrow) int [nkfrag];
   nangle_frag      = new (std::nothrow) int [nkfrag];
   bond_start_frag  = new (std::nothrow) int [nkfrag];
   angle_start_frag = new (std::nothrow) int [nkfrag];

   indxfrag_start = new (std::nothrow) int [nfrag];
   kfrag          = new (std::nothrow) int [nfrag];

   self_interaction = new (std::nothrow) bool [nkfrag];

   // bond and angle spring parameters
   bond_indx  = new (std::nothrow) int [2*nbond],
   angle_indx = new (std::nothrow) int [3*nangle];
   Krb        = new (std::nothrow) double [2*nbond];
   Kra        = new (std::nothrow) double [2*nangle];

   double pi = 4.0*std::atan(1.0);

   int bond_start_sum  = 0;
   int angle_start_sum = 0;
   int ntf = 0;
   int ia  = 0;
   for (auto & frag: fragments)
   {
      frag = mystring_rtrim(mystring_ltrim(mystring_split(frag,"end")[0]));

      int tsize = std::stoi(mystring_split0(mystring_split(frag,"size")[1])[0]);

      size_frag[ia] = tsize;

      if (mystring_contains(mystring_lowercase(frag),"index_start"))
      {
         std::vector<std::string> ss = mystring_split0(mystring_trim(mystring_split(mystring_split(frag,"index_start")[1],"\n")[0]));
         nfrag_frag[ia] = ss.size();
         for (auto i=0; i<ss.size(); ++i)
         {
            indxfrag_start[ntf] = std::stoi(ss[i])-1;
            kfrag[ntf] = ia;
            ++ntf;
         }
      }

      if (mystring_contains(mystring_lowercase(frag),"bond_spring"))
      {
         std::string line;
         std::vector<std::string> ss;
         nbond_frag[ia]      = (mystring_split(frag,"bond_spring").size() - 1);
         bond_start_frag[ia] = bond_start_sum;
         bond_start_sum     += nbond_frag[ia];
         for (auto b=0; b<nbond_frag[ia]; ++b)
         {
            int b1 = b + bond_start_frag[ia];
            line = mystring_split(frag,"bond_spring")[b+1];
            ss   = mystring_split0(line);
            int i = std::stoi(ss[0]);
            int j = std::stoi(ss[1]);
            double k = std::stod(ss[2]);
            double r = std::stod(ss[3]);

            bond_indx[2*b1  ] = i-1;
            bond_indx[2*b1+1] = j-1;
            Krb[2*b1  ] = k/27.2116/23.06/ANGTOBOHR/ANGTOBOHR;
            Krb[2*b1+1] = r*ANGTOBOHR;
         }
      }

      if (mystring_contains(mystring_lowercase(frag),"angle_spring"))
      {
         std::string line;
         std::vector<std::string> ss;
         nangle_frag[ia]      = (mystring_split(frag,"angle_spring").size() - 1);
         angle_start_frag[ia] = angle_start_sum;
         angle_start_sum     += nangle_frag[ia];
         for (auto a=0; a<nangle_frag[ia]; ++a)
         {
            int a1 = a + angle_start_frag[ia];
            line = mystring_split(frag,"angle_spring")[a+1];
            ss   = mystring_split0(line);
            int i = std::stoi(ss[0]);
            int j = std::stoi(ss[1]);
            int k = std::stoi(ss[2]);
            double ks    = std::stod(ss[3]);
            double theta = std::stod(ss[4]);
            std::cout << "i=" << i << " j=" << j << " k=" << k
                      << " ks =" << ks << " theta=" << theta << std::endl;
            angle_indx[3*a1  ] = i-1;
            angle_indx[3*a1+1] = j-1;
            angle_indx[3*a1+2] = k-1;
            Kra[2*a1  ] = ks/27.2116/23.06;
            Kra[2*a1+1] = theta*pi/180.0;
         }
      }
      ++ia;
   }
}


/*******************************************
 *                                         *
 *      QMMM_Operator::spring_Energy       *
 *                                         *
 *******************************************/
double QMMM_Operator::spring_Energy(const double rion[])
{
   double   E = 0.0;
   for (auto w1=0; w1<nfrag; ++w1)
   {
      int ks1 = indxfrag_start[w1];
      int ia  = kfrag[w1];
      int nbs = nbond_frag[ia];
      int nas = nangle_frag[ia];
      if (nbs>0)
      {
         int bs1 = bond_start_frag[ia];
         E += spring_bond_frag(nbs,&bond_indx[2*bs1],&Krb[2*bs1],&rion[3*ks1]);
      }
      if (nas>0)
      {
         int ba1 = angle_start_frag[ia];
         E += spring_angle_frag(nas,&angle_indx[3*ba1],&Kra[2*ba1],&rion[3*ks1]);
      }
   }
   return E;
}

/**************************************************
 *                                                *
 *           QMMM_Operator::spring_Force          *
 *                                                *
 **************************************************/
void QMMM_Operator::spring_Force( const double rion[], double fion[])
{
   for (auto w1=0; w1<nfrag; ++w1)
   {
      int ks1 = indxfrag_start[w1];
      int ia  = kfrag[w1];
      int nbs = nbond_frag[ia];
      int nas = nangle_frag[ia];
      if (nbs>0)
      {
         int bs1 = bond_start_frag[ia];
         spring_bond_force_frag(nbs,&bond_indx[2*bs1],&Krb[2*bs1],&rion[3*ks1],&fion[3*ks1]);
      }
      if (nas>0)
      {
         int ba1 = angle_start_frag[ia];
         spring_angle_force_frag(nas,&angle_indx[3*ba1],&Kra[2*ba1],&rion[3*ks1],&fion[3*ks1]);
      }
   }
}

/****************************************************
 *                                                  *
 *        QMMM_Operator::Coulomb_Energy             *
 *                                                  *
 ****************************************************/
double QMMM_Operator::Coulomb_Energy(const double qion[], const double rion[])
{
   double E=0.0;

   // total q-q interactions
   for (auto ii=0; ii<(nion-1); ++ii)
      for (auto jj=ii+1; jj<nion; ++jj)
         E += Q_Electrostatic_self(&rion[3*ii],qion[ii], 
                                   &rion[3*jj],qion[jj]);

   return E;
}

/****************************************************
 *                                                  *
 *        QMMM_Operator::Coulomb_Force              *
 *                                                  *
 ****************************************************/
void QMMM_Operator::Coulomb_Force(const double qion[], const double rion[], double fion[])
{
   // total q-q interactions
   for (auto ii=0;       ii<(nion-1); ++ii)
      for (auto jj=ii+1; jj<nion;     ++jj)
         Q_Electrostatic_Force_self(&rion[3*ii],qion[ii],&fion[3*ii],
                                    &rion[3*jj],qion[jj],&fion[3*jj]);
}



/****************************************************
 *                                                  *
 *        MMMM_electrostatic_energy                 *
 *                                                  *
 ****************************************************/

// Function partially removes Coulomb frag-frag energies, and completely removes Coulomb frag self energies.
// This function is intended to be called in conjunction with Coulomb_Energy.

double QMMM_Operator::MMMM_electrostatic_energy(const double qion[], const double rion[])
{
   double E=0.0;
   double rw1_cm[3],rw2_cm[3];

   // Partialy removes MM/MM Coulomb frag-frag energies 
   for (auto w1=0; w1<(nfrag-1); ++ w1)
   {
      int ks1 = indxfrag_start[w1];
      int n1  = size_frag[w1];
      double Rin1  = switch_Rin[kfrag[w1]];
      double Rout1 = switch_Rin[kfrag[w1]];
      Q_cm(n1,ks1,mass,rion,rw1_cm);

      for (auto w2=w1+1; w2<nfrag; ++w2)
      {
         int ks2 = indxfrag_start[w2];
         int n2  = size_frag[w2];
         double Rin2  = switch_Rin[kfrag[w2]];
         double Rout2 = switch_Rin[kfrag[w2]];
         Q_cm(n2,ks2,mass,rion,rw2_cm);

         double Rin  = 0.5*(Rin1 +Rin2);
         double Rout = 0.5*(Rout1+Rout2);
         E += Q_E_frag_frag(n1,ks1,rw1_cm,
                            n2,ks2,rw2_cm,
                            Rin,Rout,
                            qion,rion);
      }
   }


   // take out MM-MM Coulomb self energy ****
   for (auto w1=0; w1<nfrag; ++w1)
   {
      if (!self_interaction[kfrag[w1]])
      {
         int ks1 = indxfrag_start[w1];
         for (auto a=0;   a<(size_frag[w1]-1); ++a)
         for (auto b=a+1; b<size_frag[w1];     ++b)
         {
            int kk1 = ks1+a;
            int kk2 = ks1+b;
            E -= Q_Electrostatic_self(&rion[3*kk1],qion[kk1],
                                      &rion[3*kk2],qion[kk2]);
         }
      }
   }

   return E;
}

/****************************************************
 *                                                  *
 *         MMMM_Coulomb_electrostatic_force         *
 *                                                  *
 ****************************************************/

void QMMM_Operator::MMMM_electrostatic_force(const double qion[], const double rion[], double fion[])
{
   double m1,m2;
   double rw1_cm[3],rw2_cm[3];

   // Partialy removes MM/MM Coulomb frag-frag energies 
   for (auto w1=0; w1<(nfrag-1); ++ w1)
   {
      int ks1 = indxfrag_start[w1];
      int n1  = size_frag[w1];
      double Rin1  = switch_Rin[kfrag[w1]];
      double Rout1 = switch_Rin[kfrag[w1]];
      m1 = Q_cmm(n1,ks1,mass,rion,rw1_cm);

      for (auto w2=w1+1; w2<nfrag; ++w2)
      {
         int ks2 = indxfrag_start[w2];
         int n2  = size_frag[w2];
         double Rin2  = switch_Rin[kfrag[w2]];
         double Rout2 = switch_Rin[kfrag[w2]];
         m2 = Q_cmm(n2,ks2,mass,rion,rw2_cm);

         double Rin  = 0.5*(Rin1 +Rin2);
         double Rout = 0.5*(Rout1+Rout2);
         Q_fion_frag_frag(n1,ks1,rw1_cm,m1,
                          n2,ks2,rw2_cm,m2,
                          Rin,Rout,
                          qion,mass,rion,fion);
      }
   }

   // take out MM-MM Coulomb self energy ****
   for (auto w1=0; w1<nfrag; ++w1)
   {
      if (!self_interaction[kfrag[w1]])
      {
         int ks1 = indxfrag_start[w1];
         for (auto a=0;   a<(size_frag[w1]-1); ++a)
         for (auto b=a+1; b<size_frag[w1];     ++b)
         {
            int kk1 = ks1+a;
            int kk2 = ks1+b;
            Q_Electrostatic_mForce_self(&rion[3*kk1],qion[kk1],&fion[3*kk1],
                                        &rion[3*kk2],qion[kk2],&fion[3*kk2]);
         }
      }
   }

}








/****************************************************
 *                                                  *
 *        QMMM_electrostatic_potential              *
 *                                                  *
 ****************************************************/
void QMMM_Operator::QMMM_electrostatic_potential(const double qion[], const double rion[], double uion[])
{
   for (auto qm=0; qm<nion_qm; ++qm)
   {
      uion[qm] = 0.0;
      for (auto mm=nion_qm; mm<nion; ++mm)
         uion[qm] += Q_Electrostatic_potential(&rion[3*mm],qion[mm], &rion[3*qm]);
   }
}


/****************************************************
 *                                                  *
 *        QMMM_electrostatic_energy                 *
 *                                                  *
 ****************************************************/
double QMMM_Operator::QMMM_electrostatic_energy(const double qion[], const double rion[])
{
   double E=0.0;

   // total qmmm interactions
   for (auto ii=0; ii<(nion_qm); ++ii)
      for (auto jj=nion_qm; jj<nion; ++jj)
         E += Q_Electrostatic_self(&rion[3*ii],qion[ii], 
                                   &rion[3*jj],qion[jj]);

   return E;
}

/****************************************************
 *                                                  *
 *        QMMM_electrostatic_force                  *
 *                                                  *
 ****************************************************/
void QMMM_Operator::QMMM_electrostatic_force(const double qion[], const double rion[], double fion[])
{
   // total qmmm interactions
   for (auto ii=0; ii<(nion_qm); ++ii)
      for (auto jj=nion_qm; jj<nion; ++jj)
         Q_Electrostatic_Force_self(&rion[3*ii],qion[ii],&fion[3*ii],
                                    &rion[3*jj],qion[jj],&fion[3*jj]);

}


/****************************************************
 *                                                  *
 *        QMQM_electrostatic_energy                 *
 *                                                  *
 ****************************************************/
double QMMM_Operator::QMQM_electrostatic_energy(const double qion[], const double rion[])
{
   double E=0.0;

   // total qmqm interactions
   for (auto ii=0; ii<(nion_qm-1); ++ii)
      for (auto jj=ii+1; jj<nion_qm; ++jj)
         E += Q_Electrostatic_self(&rion[3*ii],qion[ii],
                                   &rion[3*jj],qion[jj]);

   return E;
}


/****************************************************
 *                                                  *
 *        QMQM_electrostatic_force                  *
 *                                                  *
 ****************************************************/
void QMMM_Operator::QMQM_electrostatic_force(const double qion[], const double rion[], double fion[])
{
   // total qmqm interactions
   for (auto ii=0; ii<(nion_qm-1); ++ii)
      for (auto jj=ii+1; jj<nion_qm; ++jj)
         Q_Electrostatic_Force_self(&rion[3*ii],qion[ii],&fion[3*ii],
                                    &rion[3*jj],qion[jj],&fion[3*jj]);

}




/****************************************************
 *                                                  *
 *               QMMM_LJ_energy                     *
 *                                                  *
 ****************************************************/
double QMMM_Operator::QMMM_LJ_energy(const double rion[])
{
   double E = 0.0;
   // QMQM LJ == 0
   // QMMM LJ
   // qm and mm interactions
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm)
      {
         double epsilon12 = std::sqrt(epsilon[qm]*epsilon[mm]);
         double sigma12 = 0.5*(sigma[qm] + sigma[mm]);
         double elj = LJ_energy(epsilon12,sigma12,&rion[3*qm],&rion[3*mm]);
         E += elj;
      }

   return E;
}

/****************************************************
 *                                                  *
 *               QMMM_LJ_force                      *
 *                                                  *
 ****************************************************/
void QMMM_Operator::QMMM_LJ_force(const double rion[], double fion[])
{
   // qm and mm interactions
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm) {
         double epsilon12 = std::sqrt(epsilon[qm]*epsilon[mm]);
         double sigma12 = 0.5*(sigma[qm] + sigma[mm]);
         LJ_force(epsilon12,sigma12,&rion[3*qm],&fion[3*qm],&rion[3*mm],&fion[3*mm]);
      }
}


/****************************************************
 *                                                  *
 *               MMMM_LJ_energy                     *
 *                                                  *
 ****************************************************/
double QMMM_Operator::MMMM_LJ_energy(const double rion[])
{
   double E = 0.0;

   // mm interactions
   for (auto mm1=nion_qm; mm1<(nion-1); ++mm1)
      for (auto mm2=mm1+1; mm2<nion; ++mm2)
      {
         double epsilon12 = std::sqrt(epsilon[mm1]*epsilon[mm2]);
         double sigma12 = 0.5*(sigma[mm1] + sigma[mm2]);
         double elj = LJ_energy(epsilon12,sigma12,&rion[3*mm1],&rion[3*mm2]);
         E += elj;
      }
   return E;
}


/****************************************************
 *                                                  *
 *               MMMM_LJ_force                      *
 *                                                  *
 ****************************************************/
void QMMM_Operator::MMMM_LJ_force(const double rion[], double fion[])
{
   // mm interactions
   for (auto mm1=nion_qm; mm1<(nion-1); ++mm1)
      for (auto mm2=mm1+1; mm2<nion; ++mm2)
      {
         double epsilon12 = std::sqrt(epsilon[mm1]*epsilon[mm2]);
         double sigma12 = 0.5*(sigma[mm1] + sigma[mm2]);
         LJ_force(epsilon12,sigma12,&rion[3*mm1],&fion[3*mm1],&rion[3*mm2],&fion[3*mm2]);
      }
}


/****************************************************
 *                                                  *
 *                    KE_ion                        *
 *                                                  *
 ****************************************************/
double QMMM_Operator::KE_ion(const double vion[])
{
      double KE = 0.0;
      for (auto ii=0; ii<nion; ++ii)
      {
         double vx = vion[3*ii];
         double vy = vion[3*ii+1];
         double vz = vion[3*ii+2];
         KE += mass[ii]*(vx*vx + vy*vy + vz*vz);
      }
      return 0.5*KE;
}










