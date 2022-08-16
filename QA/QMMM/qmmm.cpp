
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















/*******************************************
 *                                         *
 *        QMMM_Operator::QMMM_Operator     *
 *                                         *
 *******************************************/
QMMM_Operator::QMMM_Operator(std::string nwinput)
{
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


