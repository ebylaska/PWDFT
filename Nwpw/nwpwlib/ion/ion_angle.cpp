
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "iofmt.hpp"
#include "ion_angle.hpp"
#include "Control2.hpp"

namespace pwdft {

static void get_ub(double *unita, double *ub)
{
   ub[0] = unita[4]*unita[8] - unita[5]*unita[7];
   ub[1] = unita[5]*unita[6] - unita[3]*unita[8];
   ub[2] = unita[3]*unita[7] - unita[4]*unita[6];
   ub[3] = unita[7]*unita[2] - unita[8]*unita[1];
   ub[4] = unita[8]*unita[0] - unita[6]*unita[2];
   ub[5] = unita[6]*unita[1] - unita[7]*unita[0];
   ub[6] = unita[1]*unita[5] - unita[2]*unita[4];
   ub[7] = unita[2]*unita[3] - unita[0]*unita[5];
   ub[8] = unita[0]*unita[4] - unita[1]*unita[3];
   double volume = unita[0]*ub[0] + unita[1]*ub[1] + unita[2]*ub[2];
   for (auto i=0; i<9; ++i)
      ub[i] /= volume;
}


  /* Constructors */
ion_angle::ion_angle(double *rion1, Control2 &control) 
{
   rion = rion1;
   nha = control.nha_angle();
   angle_exists = (nha>0);
   periodic = (control.version==3);

   if (angle_exists) {
      ua[0] = control.unita(0,0); ua[1] = control.unita(1,0); ua[2] = control.unita(2,0);
      ua[3] = control.unita(0,1); ua[4] = control.unita(1,1); ua[5] = control.unita(2,1);
      ua[6] = control.unita(0,2); ua[7] = control.unita(1,2); ua[8] = control.unita(2,2);
      get_ub(ua,ub);

      i0 = new (std::nothrow) int[nha]();
      j0 = new (std::nothrow) int[nha]();
      k0 = new (std::nothrow) int[nha]();
      Kspring0 = new (std::nothrow) double [nha]();
      Theta0   = new (std::nothrow) double [nha]();
      for (auto i=0; i<nha; ++i)
      {
         i0[i] = control.i0_angle(i);
         j0[i] = control.j0_angle(i);
         k0[i] = control.k0_angle(i);
         Kspring0[i] = control.Kspring0_angle(i);
         Theta0[i]   = control.Theta0_angle(i);
      }
   }
}


  /* functions */

  /********************************
   *                              *
   *         min_diff_xyz         *
   *                              *
   ********************************/
   void ion_angle::min_diff_xyz(double *x, double *y, double *z) 
   {
      if (periodic)
      {
         double c0 = (*x)*ub[0] + (*y)*ub[1] + (*z)*ub[2];
         double c1 = (*x)*ub[3] + (*y)*ub[4] + (*z)*ub[5];
         double c2 = (*x)*ub[6] + (*y)*ub[7] + (*z)*ub[8];
         *x = ua[0]*c0 + ua[3]*c1 + ua[6]*c2;
         *y = ua[1]*c0 + ua[4]*c1 + ua[7]*c2;
         *z = ua[2]*c0 + ua[5]*c1 + ua[8]*c2;
      }
   }


  /********************************************
   *                                          *
   *        ion_angle::spring_energy          *
   *                                          *
   ********************************************/
   /**
    * Calculates the energy for the specified spring bonding.
    *
    * @param i The index of the spring bonding.
    * @return The calculated energy of the spring bonding.
    */
   double ion_angle::spring_energy(const int i)
   {
      double ea = 0.0;
      auto ii0 = i0[i]-1;
      auto jj0 = j0[i]-1;
      auto kk0 = k0[i]-1;

      auto x1 = rion[3*ii0]   - rion[3*jj0];
      auto y1 = rion[3*ii0+1] - rion[3*jj0+1];
      auto z1 = rion[3*ii0+2] - rion[3*jj0+2];
      min_diff_xyz(&x1,&y1,&z1);
      double R1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);

      auto x2 = rion[3*kk0]   - rion[3*jj0];
      auto y2 = rion[3*kk0+1] - rion[3*jj0+1];
      auto z2 = rion[3*kk0+2] - rion[3*jj0+2];
      min_diff_xyz(&x2,&y2,&z2);
      double R2 = std::sqrt(x2*x2 + y2*y2 + z2*z2);

      auto denom = R1*R2;
      if (denom>1.0e-11) 
      {
          auto ctheta = (x1*x2+y1*y2+z1*z2)/(denom);
          if (ctheta >  1.0) ctheta = 1.0;
          if (ctheta < -1.0) ctheta = -1.0;
          auto Theta = std::acos(ctheta);
          ea =  Kspring0[i]*std::pow((Theta-Theta0[i]),2);
      }
      return ea; 
   }


  /**************************************
   *                                    *
   *        ion_angle::energy           *
   *                                    *
   **************************************/
  /**
   * Calculates the energy of the system based on the input parameters.
   *
   * @return The calculated energy for the bond spring constraint of the system;
   */
   double ion_angle::energy()
   {
      double ea = 0.0;
      for (auto i=0; i<nha; ++i)
         ea += spring_energy(i);
      return ea;
   }

  /**************************************
   *                                    *
   *      ion_angle::energyfion         *
   *                                    *
   **************************************/
  /**
   * Calculates the energy and forces on a set of ions based for the bond spring constraint of the system.
   *
   * @param fion updates the values in for the pointer to a double array representing the forces on ions.
   * @return The calculated energy of the system and adds the bond bondforce to fion.. 
   */
  double ion_angle::energyfion(double *fion) 
  {
     double ea=0.0;
     for (auto i=0; i<nha; ++i)
     {
        auto ii0 = i0[i]-1;
        auto jj0 = j0[i]-1;
        auto kk0 = k0[i]-1;

        auto x1 = rion[3*ii0]   - rion[3*jj0];
        auto y1 = rion[3*ii0+1] - rion[3*jj0+1];
        auto z1 = rion[3*ii0+2] - rion[3*jj0+2];
        min_diff_xyz(&x1,&y1,&z1);
        double R1sqr = (x1*x1 + y1*y1 + z1*z1);
        double R1 = std::sqrt(R1sqr);

        auto x2 = rion[3*kk0]   - rion[3*jj0];
        auto y2 = rion[3*kk0+1] - rion[3*jj0+1];
        auto z2 = rion[3*kk0+2] - rion[3*jj0+2];
        min_diff_xyz(&x2,&y2,&z2);
        double R2sqr = (x2*x2 + y2*y2 + z2*z2);
        double R2 = std::sqrt(R2sqr);
       
        auto denom = R1*R2;
        if (denom > 1.0e-11)
        {
           auto ctheta = (x1*x2+y1*y2+z1*z2)/(denom);
           if (ctheta >  1.0) ctheta =  1.0;
           if (ctheta < -1.0) ctheta = -1.0;
           auto stheta = std::sqrt(1.0 - ctheta*ctheta);
           if (stheta < 0.001) stheta = 0.001;
           stheta = 1.0/stheta;

           auto q      = std::acos(ctheta) - Theta0[i];
           auto tk     = Kspring0[i]*q;
           auto aa     = 2.0*tk*stheta;
           auto a11    =  aa*ctheta/R1sqr;
           auto a12    = -aa/(denom);
           auto a22    =  aa*ctheta/R2sqr;

           auto vx1 = a11*x1 + a12*x2;
           auto vx2 = a22*x2 + a12*x1;

           auto vy1 = a11*y1 + a12*y2;
           auto vy2 = a22*y2 + a12*y1;

           auto vz1 = a11*z1 + a12*z2;
           auto vz2 = a22*z2 + a12*z1;

           ea +=  Kspring0[i]*q*q;
       
           fion[3*ii0]   -= vx1;
           fion[3*ii0+1] -= vy1;
           fion[3*ii0+2] -= vz1;

           fion[3*jj0]   += vx1 + vx2;
           fion[3*jj0+1] += vy1 + vy2;
           fion[3*jj0+2] += vz1 + vz2;

           fion[3*kk0]   -= vx2;
           fion[3*kk0+1] -= vy2;
           fion[3*kk0+2] -= vz2;
        }
     }
     return ea;
  }


  /**
   * Overloaded output stream operator for printing the ion_angle object to an output stream.
   **
   ** @param os The output stream to write to.
   ** @param bc The ion_angle object to be printed.
   ** @return The modified output stream.
   **/
   /*******************************************
    *                                         *
    *         ion_angle::print_all            *
    *                                         *
    *******************************************/
   std::string ion_angle::print_all(const int opt)
   {
     if (nha==0) 
        return "";
     else
     {
        std::stringstream stream;

        for (auto i=0; i<nha; ++i)
        {
           auto ii0 = i0[i]-1;
           auto jj0 = j0[i]-1;
           auto kk0 = k0[i]-1;
           auto x1 = rion[3*ii0]   - rion[3*jj0];
           auto y1 = rion[3*ii0+1] - rion[3*jj0+1];
           auto z1 = rion[3*ii0+2] - rion[3*jj0+2];
           min_diff_xyz(&x1,&y1,&z1);
           double R1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);
           auto x2 = rion[3*kk0]   - rion[3*jj0];
           auto y2 = rion[3*kk0+1] - rion[3*jj0+1];
           auto z2 = rion[3*kk0+2] - rion[3*jj0+2];
           min_diff_xyz(&x2,&y2,&z2);
           double R2 = std::sqrt(x2*x2 + y2*y2 + z2*z2);

           double Theta = 0.0;
           auto denom = R1*R2;
           if (denom>1.0e-11) 
           {
               auto ctheta = (x1*x2+y1*y2+z1*z2)/(denom);
               if (ctheta >  1.0) ctheta = 1.0;
               if (ctheta < -1.0) ctheta = -1.0;
               Theta = std::acos(ctheta);
           }
           double espring=this->spring_energy(i);
 
           if (opt==0)
           {
              stream << "      bond spring #" << Ifmt(5) << i+1 << std::endl;
              stream << "      spring parameters:" << std::endl;
              stream << "         i0         =" << Ifmt(12) << this->i0[i] << std::endl;
              stream << "         j0         =" << Ifmt(12) << this->j0[i] << std::endl;
              stream << "         k0         =" << Ifmt(12) << this->k0[i] << std::endl;
              stream << "         Kspring0   =" << Ffmt(12,6) << this->Kspring0[i] << std::endl;
              stream << "         Theta0     =" << Ffmt(12,6) << this->Theta0[i] << std::endl;
              stream << "      Theta         =" << Ffmt(12,6) << Theta   << std::endl;
              stream << "      spring energy =" << Ffmt(12,6) << espring << std::endl;
           }
           if (opt==1)
           {
              stream << " bond spring #" << Ifmt(5) << i+1 << std::endl;
              stream << " spring parameters:" << std::endl;
              stream << "    i0         =" << Ifmt(12) << this->i0[i] << std::endl;
              stream << "    j0         =" << Ifmt(12) << this->j0[i] << std::endl;
              stream << "    k0         =" << Ifmt(12) << this->k0[i] << std::endl;
              stream << "    Kspring0   =" << Ffmt(12,6) << this->Kspring0[i] << std::endl;
              stream << "    Theta0     =" << Ffmt(12,6) << this->Theta0[i]   << std::endl;
              stream << " Theta         =" << Ffmt(12,6) << Theta   << std::endl;
              stream << " spring energy =" << Ffmt(12,6) << espring << std::endl;
           }

           stream << std::endl;
        }
        return stream.str();
     }
  }

} // namespace pwdft

