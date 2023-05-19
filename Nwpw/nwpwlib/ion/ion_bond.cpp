
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "iofmt.hpp"
#include "ion_bond.hpp"
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
ion_bond::ion_bond(double *rion1, Control2 &control) 
{
   rion = rion1;
   nhb = control.nhb_bond();
   bond_exists = (nhb>0);
   periodic = (control.version==3);

   if (bond_exists) {
      ua[0] = control.unita(0,0); ua[1] = control.unita(1,0); ua[2] = control.unita(2,0);
      ua[3] = control.unita(0,1); ua[4] = control.unita(1,1); ua[5] = control.unita(2,1);
      ua[6] = control.unita(0,2); ua[7] = control.unita(1,2); ua[8] = control.unita(2,2);
      get_ub(ua,ub);

      i0 = new (std::nothrow) int[nhb]();
      j0 = new (std::nothrow) int[nhb]();
      K0 = new (std::nothrow) double [nhb]();
      R0 = new (std::nothrow) double [nhb]();
      for (auto i=0; i<nhb; ++i)
      {
         i0[i] = control.i0_bond(i);
         j0[i] = control.j0_bond(i);
         K0[i] = control.K0_bond(i);
         R0[i] = control.R0_bond(i);
      }
   }
}


  /* functions */

  /********************************
   *                              *
   *         min_diff_xyz         *
   *                              *
   ********************************/
   void ion_bond::min_diff_xyz(double *x, double *y, double *z) 
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
   *        ion_bond::spring_energy           *
   *                                          *
   ********************************************/
   /**
    * Calculates the energy for the specified spring bonding.
    *
    * @param i The index of the spring bonding.
    * @return The calculated energy of the spring bonding.
    */
   double ion_bond::spring_energy(const int i)
   {
       auto ii0 = i0[i]-1;
       auto jj0 = j0[i]-1;
       auto x = rion[3*ii0]   - rion[3*jj0];
       auto y = rion[3*ii0+1] - rion[3*jj0+1];
       auto z = rion[3*ii0+2] - rion[3*jj0+2];
       min_diff_xyz(&x,&y,&z);
       double R = std::sqrt(x*x + y*y + z*z);
    
       return 0.5*K0[i]*std::pow((R-R0[i]),2);
   }


  /**************************************
   *                                    *
   *        ion_bond::energy            *
   *                                    *
   **************************************/
  /**
   * Calculates the energy of the system based on the input parameters.
   *
   * @return The calculated energy for the bond spring constraint of the system;
   */
   double ion_bond::energy()
   {
      double eb = 0.0;
      for (auto i=0; i<nhb; ++i)
         eb += spring_energy(i);
      return eb;
   }

  /**************************************
   *                                    *
   *      ion_bond::energyfion          *
   *                                    *
   **************************************/
  /**
   * Calculates the energy and forces on a set of ions based for the bond spring constraint of the system.
   *
   * @param fion updates the values in for the pointer to a double array representing the forces on ions.
   * @return The calculated energy of the system and adds the bond bondforce to fion.. 
   */
  double ion_bond::energyfion(double *fion) 
  {
     double eb=0.0;
     for (auto i=0; i<nhb; ++i)
     {
        auto ii0 = i0[i]-1;
        auto jj0 = j0[i]-1;
        auto x = rion[3*ii0]   - rion[3*jj0];
        auto y = rion[3*ii0+1] - rion[3*jj0+1];
        auto z = rion[3*ii0+2] - rion[3*jj0+2];
        min_diff_xyz(&x,&y,&z);
        double R = std::sqrt(x*x + y*y + z*z);
       
        eb +=  0.5*K0[i]*std::pow((R-R0[i]),2);
        double deb = K0[i]*(R-R0[i]);
       
        fion[3*ii0]   += deb*(x/R);
        fion[3*ii0+1] += deb*(y/R);
        fion[3*ii0+2] += deb*(z/R);
        fion[3*jj0]   -= deb*(x/R);
        fion[3*jj0+1] -= deb*(y/R);
        fion[3*jj0+2] -= deb*(z/R);
     }
     return eb;
  }


  /**
   * Overloaded output stream operator for printing the ion_bond object to an output stream.
   **
   ** @param os The output stream to write to.
   ** @param bc The ion_bond object to be printed.
   ** @return The modified output stream.
   **/
   /*******************************************
    *                                         *
    *         ion_bond::print_all             *
    *                                         *
    *******************************************/
   std::string ion_bond::print_all(const int opt)
   {
     if (nhb==0) 
        return "";
     else
     {
        std::stringstream stream;

        for (auto i=0; i<nhb; ++i)
        {
           auto ii0 = i0[i]-1;
           auto jj0 = j0[i]-1;
           auto x = rion[3*ii0]   - rion[3*jj0];
           auto y = rion[3*ii0+1] - rion[3*jj0+1];
           auto z = rion[3*ii0+2] - rion[3*jj0+2];
           min_diff_xyz(&x,&y,&z);
           double R = std::sqrt(x*x + y*y + z*z);
           double espring=this->spring_energy(i);
 
           if (opt==0)
           {
              stream << "      bond spring #" << Ifmt(5) << i+1 << std::endl;
              stream << "      spring parameters:" << std::endl;
              stream << "         i0         =" << Ifmt(12) << this->i0[i] << std::endl;
              stream << "         j0         =" << Ifmt(12) << this->j0[i] << std::endl;
              stream << "         K0         =" << Ffmt(12,6) << this->K0[i] << std::endl;
              stream << "         R0         =" << Ffmt(12,6) << this->R0[i] << std::endl;
              stream << "      R             =" << Ffmt(12,6) << R   << std::endl;
              stream << "      spring energy =" << Ffmt(12,6) << espring << std::endl;
           }
           if (opt==1)
           {
              stream << " bond spring #" << Ifmt(5) << i+1 << std::endl;
              stream << " spring parameters:" << std::endl;
              stream << "    i0         =" << Ifmt(12) << this->i0[i] << std::endl;
              stream << "    j0         =" << Ifmt(12) << this->j0[i] << std::endl;
              stream << "    K0         =" << Ffmt(12,6) << this->K0[i] << std::endl;
              stream << "    R0         =" << Ffmt(12,6) << this->R0[i] << std::endl;
              stream << " R             =" << Ffmt(12,6) << R   << std::endl;
              stream << " spring energy =" << Ffmt(12,6) << espring << std::endl;
           }

           stream << std::endl;
        }
        return stream.str();
     }
  }

} // namespace pwdft

