
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "iofmt.hpp"
#include "ion_cbond.hpp"
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
ion_cbond::ion_cbond(double *rion1, Control2 &control) 
{
   rion = rion1;
   nhcb = control.nhcb_cbond();
   cbond_exists = (nhcb>0);
   periodic = (control.version==3);

   if (cbond_exists) {
      ua[0] = control.unita(0,0); ua[1] = control.unita(1,0); ua[2] = control.unita(2,0);
      ua[3] = control.unita(0,1); ua[4] = control.unita(1,1); ua[5] = control.unita(2,1);
      ua[6] = control.unita(0,2); ua[7] = control.unita(1,2); ua[8] = control.unita(2,2);
      get_ub(ua,ub);

      i0 = new (std::nothrow) int[nhcb]();
      j0 = new (std::nothrow) int[nhcb]();
      k0 = new (std::nothrow) int[nhcb]();
      l0 = new (std::nothrow) int[nhcb]();
      Kspring0 = new (std::nothrow) double [nhcb]();
      Rij0     = new (std::nothrow) double [nhcb]();
      Rkl0     = new (std::nothrow) double [nhcb]();
      for (auto i=0; i<nhcb; ++i)
      {
         i0[i] = control.i0_cbond(i);
         j0[i] = control.j0_cbond(i);
         k0[i] = control.k0_cbond(i);
         l0[i] = control.l0_cbond(i);
         Kspring0[i] = control.Kspring0_cbond(i);
         Rij0[i]     = control.Rij0_cbond(i);
         Rkl0[i]     = control.Rkl0_cbond(i);
      }
   }
}


  /* functions */

  /********************************
   *                              *
   *         min_diff_xyz         *
   *                              *
   ********************************/
   void ion_cbond::min_diff_xyz(double *x, double *y, double *z) 
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
   *        ion_cbond::spring_energy          *
   *                                          *
   ********************************************/
   /**
    * Calculates the energy for the specified spring bonding.
    *
    * @param i The index of the spring bonding.
    * @return The calculated energy of the spring bonding.
    */
   double ion_cbond::spring_energy(const int i)
   {
       auto ii0 = i0[i]-1;
       auto jj0 = j0[i]-1;
       auto kk0 = k0[i]-1;
       auto ll0 = l0[i]-1;

       auto x1 = rion[3*ii0]   - rion[3*jj0];
       auto y1 = rion[3*ii0+1] - rion[3*jj0+1];
       auto z1 = rion[3*ii0+2] - rion[3*jj0+2];
       min_diff_xyz(&x1,&y1,&z1);
       double Rij = std::sqrt(x1*x1 + y1*y1 + z1*z1);

       auto x2 = rion[3*kk0]   - rion[3*ll0];
       auto y2 = rion[3*kk0+1] - rion[3*ll0+1];
       auto z2 = rion[3*kk0+2] - rion[3*ll0+2];
       min_diff_xyz(&x2,&y2,&z2);
       double Rkl = std::sqrt(x2*x2 + y2*y2 + z2*z2);

       return Kspring0[i]*(Rij-Rij0[i])*(Rkl-Rkl0[i]);
   }


  /**************************************
   *                                    *
   *        ion_cbond::energy           *
   *                                    *
   **************************************/
  /**
   * Calculates the energy of the system based on the input parameters.
   *
   * @return The calculated energy for the bond spring constraint of the system;
   */
   double ion_cbond::energy()
   {
      double ecb = 0.0;
      for (auto i=0; i<nhcb; ++i)
         ecb += spring_energy(i);
      return ecb;
   }

  /**************************************
   *                                    *
   *      ion_cbond::energyfion         *
   *                                    *
   **************************************/
  /**
   * Calculates the energy and forces on a set of ions based for the bond spring constraint of the system.
   *
   * @param fion updates the values in for the pointer to a double array representing the forces on ions.
   * @return The calculated energy of the system and adds the bond bondforce to fion.. 
   */
  double ion_cbond::energyfion(double *fion) 
  {
     double ecb=0.0;
     for (auto i=0; i<nhcb; ++i)
     {
        auto ii0 = i0[i]-1;
        auto jj0 = j0[i]-1;
        auto kk0 = k0[i]-1;
        auto ll0 = l0[i]-1;

        auto x1 = rion[3*ii0]   - rion[3*jj0];
        auto y1 = rion[3*ii0+1] - rion[3*jj0+1];
        auto z1 = rion[3*ii0+2] - rion[3*jj0+2];
        min_diff_xyz(&x1,&y1,&z1);
        double Rij = std::sqrt(x1*x1 + y1*y1 + z1*z1);

        auto x2 = rion[3*kk0]   - rion[3*ll0];
        auto y2 = rion[3*kk0+1] - rion[3*ll0+1];
        auto z2 = rion[3*kk0+2] - rion[3*ll0+2];
        min_diff_xyz(&x2,&y2,&z2);
        double Rkl = std::sqrt(x2*x2 + y2*y2 + z2*z2);
       
        ecb +=  Kspring0[i]*(Rij-Rij0[i])*(Rkl-Rkl0[i]);

        double decb1 = Kspring0[i]*(Rij-Rij0[i]);
        double decb2 = Kspring0[i]*(Rkl-Rkl0[i]);
       
        fion[3*ii0]   -= decb1*(x1/Rij);
        fion[3*ii0+1] -= decb1*(y1/Rij);
        fion[3*ii0+2] -= decb1*(z1/Rij);

        fion[3*jj0]   += decb1*(x1/Rij);
        fion[3*jj0+1] += decb1*(y1/Rij);
        fion[3*jj0+2] += decb1*(z1/Rij);

        fion[3*kk0]   -= decb2*(x2/Rkl);
        fion[3*kk0+1] -= decb2*(y2/Rkl);
        fion[3*kk0+2] -= decb2*(z2/Rkl);

        fion[3*ll0]   += decb2*(x2/Rkl);
        fion[3*ll0+1] += decb2*(y2/Rkl);
        fion[3*ll0+2] += decb2*(z2/Rkl);
     }
     return ecb;
  }

  /**
   * Overloaded output stream operator for printing the ion_cbond object to an output stream.
   **
   ** @param os The output stream to write to.
   ** @param bc The ion_cbond object to be printed.
   ** @return The modified output stream.
   **/
   /*******************************************
    *                                         *
    *         ion_cbond::print_all             *
    *                                         *
    *******************************************/
   std::string ion_cbond::print_all(const int opt)
   {
     if (nhcb==0) 
        return "";
     else
     {
        std::stringstream stream;

        for (auto i=0; i<nhcb; ++i)
        {
           auto ii0 = i0[i]-1;
           auto jj0 = j0[i]-1;
           auto kk0 = k0[i]-1;
           auto ll0 = l0[i]-1;

           auto x1 = rion[3*ii0]   - rion[3*jj0];
           auto y1 = rion[3*ii0+1] - rion[3*jj0+1];
           auto z1 = rion[3*ii0+2] - rion[3*jj0+2];
           min_diff_xyz(&x1,&y1,&z1);
           double Rij = std::sqrt(x1*x1 + y1*y1 + z1*z1);

           auto x2 = rion[3*kk0]   - rion[3*ll0];
           auto y2 = rion[3*kk0+1] - rion[3*ll0+1];
           auto z2 = rion[3*kk0+2] - rion[3*ll0+2];
           min_diff_xyz(&x2,&y2,&z2);
           double Rkl = std::sqrt(x2*x2 + y2*y2 + z2*z2);

           double espring=this->spring_energy(i);
 
           if (opt==0)
           {
              stream << "      cbond spring #" << Ifmt(5) << i+1 << std::endl;
              stream << "      spring parameters:" << std::endl;
              stream << "         i0         =" << Ifmt(12) << this->i0[i] << std::endl;
              stream << "         j0         =" << Ifmt(12) << this->j0[i] << std::endl;
              stream << "         k0         =" << Ifmt(12) << this->k0[i] << std::endl;
              stream << "         l0         =" << Ifmt(12) << this->l0[i] << std::endl;
              stream << "         Ksping0    =" << Ffmt(12,6) << this->Kspring0[i] << std::endl;
              stream << "         Rij0       =" << Ffmt(12,6) << this->Rij0[i] << std::endl;
              stream << "         Rkl0       =" << Ffmt(12,6) << this->Rkl0[i] << std::endl;
              stream << "      Rij           =" << Ffmt(12,6) << Rij   << std::endl;
              stream << "      Rkl           =" << Ffmt(12,6) << Rkl   << std::endl;
              stream << "      spring energy =" << Ffmt(12,6) << espring << std::endl;
           }
           if (opt==1)
           {
              stream << " cbond spring #" << Ifmt(5) << i+1 << std::endl;
              stream << " spring parameters:" << std::endl;
              stream << "    i0         =" << Ifmt(12) << this->i0[i] << std::endl;
              stream << "    j0         =" << Ifmt(12) << this->j0[i] << std::endl;
              stream << "    k0         =" << Ifmt(12) << this->k0[i] << std::endl;
              stream << "    l0         =" << Ifmt(12) << this->l0[i] << std::endl;
              stream << "    Ksping0    =" << Ffmt(12,6) << this->Kspring0[i] << std::endl;
              stream << "    Rij0       =" << Ffmt(12,6) << this->Rij0[i] << std::endl;
              stream << "    Rkl0       =" << Ffmt(12,6) << this->Rkl0[i] << std::endl;
              stream << " Rij           =" << Ffmt(12,6) << Rij << std::endl;
              stream << " Rkl           =" << Ffmt(12,6) << Rkl << std::endl;
              stream << " spring energy =" << Ffmt(12,6) << espring << std::endl;
           }

           stream << std::endl;
        }
        return stream.str();
     }
  }

} // namespace pwdft

