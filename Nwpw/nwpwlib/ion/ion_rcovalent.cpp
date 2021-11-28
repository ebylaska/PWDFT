/* ion_rcovalent.cpp
   Author - Eric Bylaska
*/

#include 	<iostream>
#include 	<vector>
#include 	<map>
#include 	<list>
#include 	<string>
#include 	<iomanip>
#include 	<sstream>

#include <algorithm> 
#include <cctype>
#include <locale>

#include	<cmath>

namespace pwdft {
using namespace pwdft;


static std::map<std::string, std::vector<double> > rcovalent = { 
{"H",  {0.32,0.00,0.00,0.00}},
{"He", {0.46,0.00,0.00,0.00}},
{"Li", {1.33,1.24,0.00,0.00}},
{"Be", {1.02,0.90,0.85,0.00}},
{"B",  {0.85,0.78,0.73,0.00}},
{"C",  {0.75,0.67,0.60,0.68}},
{"N",  {0.71,0.60,0.54,0.00}},
{"O",  {0.63,0.57,0.53,0.00}},
{"F",  {0.64,0.59,0.53,0.00}},
{"Ne", {0.67,0.96,0.00,0.00}},
{"Na", {1.55,1.60,0.00,0.00}},
{"Mg", {1.39,1.32,1.27,0.00}},
{"Al", {1.26,1.13,1.11,0.00}},
{"Si", {1.16,1.07,1.02,0.00}},
{"P",  {1.11,1.02,0.94,0.00}},
{"S",  {1.03,0.94,0.95,0.00}},
{"Cl", {0.99,0.95,0.93,0.00}},
{"Ar", {0.96,1.07,0.96,0.00}},
{"K",  {1.96,1.93,0.00,0.00}},
{"Ca", {1.71,1.47,1.33,0.00}},
{"Sc", {1.48,1.16,1.14,0.00}},
{"Ti", {1.36,1.17,1.08,0.00}},
{"V",  {1.34,1.12,1.06,0.00}},
{"Cr", {1.22,1.11,1.03,0.00}},
{"Mn", {1.19,1.05,1.03,0.00}},
{"Fe", {1.16,1.09,1.02,0.00}},
{"Co", {1.11,1.03,0.96,0.00}},
{"Ni", {1.10,1.01,1.01,0.00}},
{"Cu", {1.12,1.15,1.20,0.00}},
{"Zn", {1.18,1.20,0.00,0.00}},
{"Ga", {1.24,1.16,1.21,0.00}},
{"Ge", {1.21,1.11,1.14,0.00}},
{"As", {1.21,1.14,1.06,0.00}},
{"Se", {1.16,1.07,1.07,0.00}},
{"Br", {1.14,1.09,1.10,0.00}},
{"Kr", {1.17,1.21,1.08,0.00}},
{"Rb", {2.10,2.02,0.00,0.00}},
{"Sr", {1.85,1.57,1.39,0.00}},
{"Y",  {1.63,1.30,1.24,0.00}},
{"Zr", {1.54,1.27,1.21,0.00}},
{"Nb", {1.47,1.25,1.16,0.00}},
{"Mo", {1.38,1.21,1.13,0.00}},
{"Tc", {1.28,1.20,1.10,0.00}},
{"Ru", {1.25,1.14,1.03,0.00}},
{"Rh", {1.25,1.10,1.06,0.00}},
{"Pd", {1.20,1.17,1.12,0.00}},
{"Ag", {1.28,1.39,1.37,0.00}},
{"Cd", {1.36,1.44,0.00,0.00}},
{"In", {1.42,1.36,1.46,0.00}},
{"Sn", {1.40,1.30,1.32,0.00}},
{"Sb", {1.40,1.33,1.27,0.00}},
{"Te", {1.36,1.28,1.21,0.00}},
{"I",  {1.33,1.29,1.25,0.00}},
{"Xe", {1.31,1.35,1.22,0.00}},
{"Cs", {2.32,1.96,0.00,0.00}},
{"Ba", {1.96,1.61,1.49,0.00}},
{"La", {1.80,1.39,1.39,0.00}},
{"Ce", {1.63,1.37,1.31,0.00}},
{"Pr", {1.76,1.38,1.28,0.00}},
{"Nd", {1.74,1.37,0.00,0.00}},
{"Pm", {1.73,1.35,0.00,0.00}},
{"Sm", {1.72,1.34,0.00,0.00}},
{"Eu", {1.68,1.34,0.00,0.00}},
{"Gd", {1.69,1.35,1.32,0.00}},
{"Tb", {1.68,1.35,0.00,0.00}},
{"Dy", {1.67,1.33,0.00,0.00}},
{"Ho", {1.66,1.33,0.00,0.00}},
{"Er", {1.65,1.33,0.00,0.00}},
{"Tm", {1.64,1.31,0.00,0.00}},
{"Yb", {1.70,1.29,0.00,0.00}},
{"Lu", {1.62,1.31,1.31,0.00}},
{"Hf", {1.52,1.28,1.22,0.00}},
{"Ta", {1.46,1.26,1.19,0.00}},
{"W",  {1.37,1.20,1.15,0.00}},
{"Re", {1.31,1.19,1.10,0.00}},
{"Os", {1.29,1.16,1.09,0.00}},
{"Ir", {1.22,1.15,1.07,0.00}},
{"Pt", {1.23,1.12,1.10,0.00}},
{"Au", {1.24,1.21,1.23,0.00}},
{"Hg", {1.33,1.42,0.00,0.00}},
{"Tl", {1.44,1.42,1.50,0.00}},
{"Pb", {1.44,1.35,1.37,0.00}},
{"Bi", {1.51,1.41,1.35,0.00}},
{"Po", {1.45,1.35,1.29,0.00}},
{"At", {1.47,1.38,1.38,0.00}},
{"Rn", {1.42,1.45,1.33,0.00}},
{"Fr", {2.23,2.18,0.00,0.00}},
{"Ra", {2.01,1.73,1.59,0.00}},
{"Ac", {1.86,1.53,1.40,0.00}},
{"Th", {1.75,1.43,1.36,0.00}},
{"Pa", {1.69,1.38,1.29,0.00}},
{"U",  {1.70,1.34,1.18,0.00}},
{"Np", {1.71,1.36,1.16,0.00}},
{"Pu", {1.72,1.35,0.00,0.00}},
{"Am", {1.66,1.35,0.00,0.00}},
{"Cm", {1.66,1.36,0.00,0.00}},
{"Bk", {1.68,1.39,0.00,0.00}},
{"Cf", {1.68,1.40,0.00,0.00}},
{"Es", {1.65,1.40,0.00,0.00}},
{"Fm", {1.67,0.00,0.00,0.00}},
{"Md", {1.73,1.39,0.00,0.00}},
{"No", {1.76,0.00,0.00,0.00}},
{"Lr", {1.61,1.41,0.00,0.00}},
{"Rf", {1.57,1.40,1.31,0.00}},
{"Db", {1.49,1.36,1.26,0.00}},
{"Sg", {1.43,1.28,1.21,0.00}},
{"Bh", {1.41,1.28,1.19,0.00}},
{"Hs", {1.34,1.25,1.18,0.00}},
{"Mt", {1.29,1.25,1.13,0.00}},
{"Ds", {1.28,1.16,1.12,0.00}},
{"Rg", {1.21,1.16,1.18,0.00}},
{"Cn", {1.22,1.37,1.30,0.00}},
{"Fl",  {1.43,0.00,0.00,0.00}},
{"Lv",  {1.75,0.00,0.00,0.00}},
};


// *******************************************
// *                                         *
// *           ion_bond_order                *
// *                                         *
// *******************************************
double ion_bond_order(std::vector<double> rc1, std::vector<double> rc2, double r12)
{ 
   double dd     = 0.0001;
   double cov[4] = {abs(r12-(rc1[0]+rc2[0]))/(rc1[0]+rc2[0]+dd),
                    abs(r12-(rc1[1]+rc2[1]))/(rc1[1]+rc2[1]+dd),
                    abs(r12-(rc1[2]+rc2[2]))/(rc1[2]+rc2[2]+dd),
                    abs(r12-(rc1[3]+rc2[3]))/(rc1[3]+rc2[3]+dd)};
   int imin = 0;
   double dmin = cov[0];
   if (cov[1]<dmin) {
      dmin = cov[1];
      imin = 1;
   }
   if (cov[2]<dmin) {
      dmin = cov[2];
      imin = 2;
   }
   if (cov[3]<dmin) {
      dmin = cov[3];
      imin = 3;
   }
   double b = 0;
   if (cov[imin]<.10) {
      b = 1+imin;
      if (imin==3) {b = 1.5;}
   }
   if ((b==0) && (r12-(rc1[0]+rc2[0])) < 0.0) b = 1;

   return b;
}



// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// *******************************************
// *                                         *
// *      ion_print_bond_angle_torsions      *
// *                                         *
// *******************************************
std::string ion_print_bond_angle_torsions(int nion, int *katm, char *symbol, double *rion) {

   double autoang = 0.529177;
   std::string msg = "";
   msg += "\n";
   msg += " Number of Atoms = " + std::to_string(nion) + "\n\n";

   msg += "\n";
   msg += " Units are Angstrom for bonds and degrees for angles\n";
   msg += "\n";
   msg += "      Type          I     J     K     L     M      Value\n";
   msg += "      ----------- ----- ----- ----- ----- ----- ----------\n";

   int count = 1;
   int ia,ja,ka,la;
   double xi,xj,xk,xl;
   double yi,yj,yk,yl;
   double zi,zj,zk,zl;
   double dx,dy,dz,r;
   std::vector<double> rci,rcj,rck,rcl;

   /* determining bonds */
   for (auto jj=0; jj<nion; ++jj)
      for (auto ii=jj+1; ii<nion; ++ii)
      {
         ia = katm[ii]; ja = katm[jj];

         xj = autoang*rion[3*jj]; yj = autoang*rion[3*jj+1]; zj = autoang*rion[3*jj+2];
         xi = autoang*rion[3*ii]; yi = autoang*rion[3*ii+1]; zi = autoang*rion[3*ii+2];
         dx = xi-xj;
         dy = yi-yj;
         dz = zi-zj;
         r = std::sqrt(dx*dx + dy*dy + dz*dz);

         std::string sym1 = &symbol[3*ia];
         std::string sym2 = &symbol[3*ja];
         rtrim(sym1);
         rtrim(sym2);
         rci   = rcovalent[sym1];
         rcj   = rcovalent[sym2];

         if (ion_bond_order(rci,rcj,r)>0) {
            sym1 += std::to_string(ii+1);
            sym2 += std::to_string(jj+1);
            std::stringstream stream;
            stream << std::setw(5) << count   << " Stretch "
                   << std::setw(9) << sym2
                   << std::setw(6) << sym1
                   << std::fixed   << std::setprecision(5) << std::setw(29) << r << "\n";;

            msg += stream.str();
            ++count;
         }
      }

   /* determining angles */
   double pi = std::atan(1.0)*4.0;
   double cstheta,theta;
   double dxij,dyij,dzij,rij;
   double dxkj,dykj,dzkj,rkj;
   for (auto jj=0; jj<nion; ++jj)
      for (auto ii=0; ii<nion; ++ii)
         for (auto kk=ii+1; kk<nion; ++kk)
            if ((kk!=ii) && (kk!=jj) && (ii!=jj)) {
               ia = katm[ii]; ja = katm[jj]; ka = katm[kk];

               xi = autoang*rion[3*ii]; yi = autoang*rion[3*ii+1]; zi = autoang*rion[3*ii+2];
               xj = autoang*rion[3*jj]; yj = autoang*rion[3*jj+1]; zj = autoang*rion[3*jj+2];
               xk = autoang*rion[3*kk]; yk = autoang*rion[3*kk+1]; zk = autoang*rion[3*kk+2];
               dxij = xi-xj;
               dyij = yi-yj;
               dzij = zi-zj;
               rij = std::sqrt(dxij*dxij + dyij*dyij + dzij*dzij);

               dxkj = xk-xj;
               dykj = yk-yj;
               dzkj = zk-zj;
               rkj = std::sqrt(dxkj*dxkj + dykj*dykj + dzkj*dzkj);

               std::string sym1 = &symbol[3*ia];
               std::string sym2 = &symbol[3*ja];
               std::string sym3 = &symbol[3*ka];
               rtrim(sym1);
               rtrim(sym2);
               rtrim(sym3);

               rci   = rcovalent[sym1];
               rcj   = rcovalent[sym2];
               rck   = rcovalent[sym3];
           
               if ((ion_bond_order(rci,rcj,rij)>0) && (ion_bond_order(rck,rcj,rkj)>0)) {
                  sym1 += std::to_string(ii+1);
                  sym2 += std::to_string(jj+1);
                  sym3 += std::to_string(kk+1);
                  cstheta = (dxij*dxkj + dyij*dykj + dzij*dzkj)/(rij*rkj);
                  theta = std::acos(cstheta)*180.0/pi;

                  std::stringstream stream;
                  stream << std::setw(5) << count   << " Bend    "
                         << std::setw(9) << sym1
                         << std::setw(6) << sym2
                         << std::setw(6) << sym3
                         << std::fixed   << std::setprecision(5) << std::setw(23) << theta << "\n";
                  msg += stream.str();
                  ++count;
               }
            }


   /* determining torsions */
   double f;
   double dxjk,dyjk,dzjk,rjk;
   double dxkl,dykl,dzkl,rkl;
   for (auto ii=0; ii<nion; ++ii)
      for (auto jj=0; jj<nion; ++jj)
         for (auto kk=0; kk<nion; ++kk)
            for (auto ll=ii+1; ll<nion; ++ll)
               if ((ii!=jj) && (ii!=kk) && (ii!=ll) && (jj!=kk) && (jj!=ll) && (kk!=ll)) {
                  ia = katm[ii]; ja = katm[jj]; ka = katm[kk]; la = katm[ll];

                  xi = autoang*rion[3*ii]; yi = autoang*rion[3*ii+1]; zi = autoang*rion[3*ii+2];
                  xj = autoang*rion[3*jj]; yj = autoang*rion[3*jj+1]; zj = autoang*rion[3*jj+2];
                  xk = autoang*rion[3*kk]; yk = autoang*rion[3*kk+1]; zk = autoang*rion[3*kk+2];
                  xl = autoang*rion[3*ll]; yl = autoang*rion[3*ll+1]; zl = autoang*rion[3*ll+2];

                  dxij = xj-xi;
                  dyij = yj-yi;
                  dzij = zj-zi;
                  rij = std::sqrt(dxij*dxij + dyij*dyij + dzij*dzij);

                  dxjk = xk-xj;
                  dyjk = yk-yj;
                  dzjk = zk-zj;
                  rjk = std::sqrt(dxjk*dxjk + dyjk*dyjk + dzjk*dzjk);

                  dxkl = xl-xk;
                  dykl = yl-yk;
                  dzkl = zl-zk;
                  rkl = std::sqrt(dxkl*dxkl + dykl*dykl + dzkl*dzkl);

                  std::string sym1 = &symbol[3*ia];
                  std::string sym2 = &symbol[3*ja];
                  std::string sym3 = &symbol[3*ka];
                  std::string sym4 = &symbol[3*la];
                  rtrim(sym1);
                  rtrim(sym2);
                  rtrim(sym3);
                  rtrim(sym4);

                  rci   = rcovalent[sym1];
                  rcj   = rcovalent[sym2];
                  rck   = rcovalent[sym3];
                  rcl   = rcovalent[sym4];

                  if ((ion_bond_order(rci,rcj,rij)>0) && (ion_bond_order(rcj,rck,rjk)>0) && (ion_bond_order(rck,rcl,rkl)>0)) {
                     dy = dxij*(dyjk*dzkl-dykl*dzjk)+dyij*(dzjk*dxkl-dzkl*dxjk)+dzij*(dxjk*dykl-dxkl*dyjk);
                     dy = dy*std::sqrt(dxjk*dxjk+dyjk*dyjk+dzjk*dzjk);
                     dx = (dyjk*dzkl-dykl*dzjk)*(dyij*dzjk-dyjk*dzij)+(dzjk*dxkl-dzkl*dxjk)*(dzij*dxjk-dzjk*dxij)+(dxjk*dykl-dxkl*dyjk)*(dxij*dyjk-dxjk*dyij);
                     f = std::atan2(dy,dx)*180.0/pi;

                     sym1 += std::to_string(ii+1);
                     sym2 += std::to_string(jj+1);
                     sym3 += std::to_string(kk+1);
                     sym4 += std::to_string(ll+1);

                     std::stringstream stream;
                     stream << std::setw(5) << count   << " Dihedral"
                            << std::setw(9) << sym1
                            << std::setw(6) << sym2
                            << std::setw(6) << sym3
                            << std::setw(6) << sym4
                            << std::fixed   << std::setprecision(5) << std::setw(17) << f << "\n";
                     msg += stream.str();
                     ++count;
                  }
               }

   msg += "\n\n";
   msg += "  XYZ format geometry\n";
   msg += "  -------------------\n";
   msg +=  std::to_string(nion) + "\n\n";

   for (auto ii=0; ii<nion; ++ii)
   {
      int ia = katm[ii];
      std::string sym1 = &symbol[3*ia];
      msg += sym1;
      msg += "    ";
      if (sym1.length()==1) msg += " ";

      std::stringstream stream;
      //stream <<  symbol[3*ia] << symbol[3*ia+1] << symbol[3*ia+2] << "     ";
      stream << std::fixed << std::setprecision(6) << std::setw(12) << autoang*rion[3*ii]   << " "
                                                   << std::setw(12) << autoang*rion[3*ii+1] << " "
                                                   << std::setw(12) << autoang*rion[3*ii+2] << "\n";
      msg += stream.str();
   }


   return msg;
}

}

