
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "iofmt.hpp"
#include "ion_bondings.hpp"
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
ion_bondings::ion_bondings(double *rion1, Control2 &control) 
{
   rion = rion1;
   nhc = control.nhc_bondings();
   bondings_exists = (nhc>0);
   periodic = (control.version==3);

   if (bondings_exists) {
      ua[0] = control.unita(0,0); ua[1] = control.unita(1,0); ua[2] = control.unita(2,0);
      ua[3] = control.unita(0,1); ua[4] = control.unita(1,1); ua[5] = control.unita(2,1);
      ua[6] = control.unita(0,2); ua[7] = control.unita(1,2); ua[8] = control.unita(2,2);
      get_ub(ua,ub);

      n0     = new (std::nothrow) int[nhc]();
      K0     = new (std::nothrow) double [nhc]();
      gamma0 = new (std::nothrow) double [nhc]();
      coef   = new (std::nothrow) double* [nhc]();
      indx   = new (std::nothrow) int* [nhc]();
      for (auto i=0; i<nhc; ++i) {
         K0[i]     = control.K0_bondings(i);
         gamma0[i] = control.gamma0_bondings(i);

         std::vector<double> tmp = control.coef_bondings(i);
         std::vector<int>   itmp = control.indx_bondings(i);
         n0[i]   = tmp.size();

         coef[i] = new (std::nothrow) double [n0[i]]();
         indx[i] = new (std::nothrow) int [2*n0[i]]();

         for (auto j=0; j<n0[i]; ++j)
         {
            coef[i][j]     = tmp[j];
            indx[i][2*j]   = itmp[2*j];
            indx[i][2*j+1] = itmp[2*j+1];
         }
      }
   }
}


  /* functions */

  /********************************
   *                              *
   *         min_diff_xyz         *
   *                              *
   ********************************/
void ion_bondings::min_diff_xyz(double *x, double *y, double *z) 
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

  /**************************************
   *                                    *
   *  ion_bondings::ij_spring_distance  *
   *                                    *
   **************************************/
   double ion_bondings::ij_spring_distance(const int i, const int j)
   {
      auto coef_tmp = coef[i][j];
      auto ii0 = indx[i][2*j]   - 1;
      auto ii1 = indx[i][2*j+1] - 1;
      auto x = rion[3*ii0]   - rion[3*ii1];
      auto y = rion[3*ii0+1] - rion[3*ii1+1];
      auto z = rion[3*ii0+2] - rion[3*ii1+2];
      min_diff_xyz(&x,&y,&z);
      double dist = std::sqrt(x*x + y*y + z*z);
      return dist;
   }

  /**************************************
   *                                    *
   *      ion_bondings::spring_gamma    *
   *                                    *
   **************************************/
   /**
    * Calculates the gamma value for the specified spring bonding.
    *
    * @param i The index of the spring bonding.
    * @return The calculated gamma value.
    */
   double ion_bondings::spring_gamma(const int i)
   {
      double gamma = 0.0;
      for (auto ii=0; ii<n0[i]; ++ii)
      {
          auto coef_tmp = coef[i][ii];
          auto ii0 = indx[i][2*ii]   - 1;
          auto ii1 = indx[i][2*ii+1] - 1;
          auto x = rion[3*ii0]   - rion[3*ii1];
          auto y = rion[3*ii0+1] - rion[3*ii1+1];
          auto z = rion[3*ii0+2] - rion[3*ii1+2];
          min_diff_xyz(&x,&y,&z);
          double dist = std::sqrt(x*x + y*y + z*z);
          gamma += coef_tmp*dist;
      }
      return gamma;
   }
  /********************************************
   *                                          *
   *      ion_bondings::spring_energy         *
   *                                          *
   ********************************************/
   /**
    * Calculates the energy for the specified spring bonding.
    *
    * @param i The index of the spring bonding.
    * @return The calculated energy of the spring bonding.
    */
   double ion_bondings::spring_energy(const int i)
   {
       double gamma = this->spring_gamma(i);
       return K0[i]*std::pow((gamma-gamma0[i]),2);
   }


  /**************************************
   *                                    *
   *      ion_bondings::energy          *
   *                                    *
   **************************************/
  /**
   * Calculates the energy of the system based on the input parameters.
   *
   * @return The calculated energy for the bondings spring constraint of the system;
   */
   double ion_bondings::energy()
   {
      double eb = 0.0;
      for (auto i=0; i<nhc; ++i)
      {
         double gamma  = this->spring_gamma(i);
         eb += K0[i]*std::pow((gamma-gamma0[i]),2);
      }
      return eb;
   }

  /**************************************
   *                                    *
   *      ion_bondings::energyfion      *
   *                                    *
   **************************************/
  /**
   * Calculates the energy and forces on a set of ions based for the bondings spring constraint of the system.
   *
   * @param fion updates the values in for the pointer to a double array representing the forces on ions.
   * @return The calculated energy of the system and adds the bond bondings force to fion.. 
   */
  double ion_bondings::energyfion(double *fion) 
  {
     double eb=0.0;
     for (auto i=0; i<nhc; ++i)
     {
        double gamma = 0.0;
        for (auto ii=0; ii<n0[i]; ++ii)
        {
            auto coef_tmp = coef[i][ii];
            auto ii0 = indx[i][2*ii]   - 1;
            auto ii1 = indx[i][2*ii+1] - 1;
            auto x = rion[3*ii0]   - rion[3*ii1];
            auto y = rion[3*ii0+1] - rion[3*ii1+1];
            auto z = rion[3*ii0+2] - rion[3*ii1+2];
            min_diff_xyz(&x,&y,&z);
            double dist = std::sqrt(x*x + y*y + z*z);
            gamma += coef_tmp*dist;
        }
        eb  +=  K0[i]*std::pow((gamma-gamma0[i]),2);
        double deb = 2*K0[i]*(gamma-gamma0[i]);
        for (auto ii=0; ii<n0[i]; ++ii)
        {
            auto coef_tmp = coef[i][ii];
            auto ii0 = indx[i][2*ii]   - 1;
            auto ii1 = indx[i][2*ii+1] - 1;
            auto x = rion[3*ii0]   - rion[3*ii1];
            auto y = rion[3*ii0+1] - rion[3*ii1+1];
            auto z = rion[3*ii0+2] - rion[3*ii1+2];
            min_diff_xyz(&x,&y,&z);
            double dist = std::sqrt(x*x + y*y + z*z);

            fion[3*ii0]   += deb*coef_tmp*(x/dist);
            fion[3*ii0+1] += deb*coef_tmp*(y/dist);
            fion[3*ii0+2] += deb*coef_tmp*(z/dist);
            fion[3*ii1]   -= deb*coef_tmp*(x/dist);
            fion[3*ii1+1] -= deb*coef_tmp*(y/dist);
            fion[3*ii1+2] -= deb*coef_tmp*(z/dist);
        }
     }
     return eb;
  }


  /**
   * Overloaded output stream operator for printing the nwpw_bondings object to an output stream.
   **
   ** @param os The output stream to write to.
   ** @param bc The nwpw_bondings object to be printed.
   ** @return The modified output stream.
   **/
   /*******************************************
    *                                         *
    *       ion_bondings::print_all          *
    *                                         *
    *******************************************/
   std::string ion_bondings::print_all(const int opt)
   {
       std::stringstream stream;

     for (auto i=0; i<this->nhc; ++i)
     {
        double espring=this->spring_energy(i);
        double gamma=this->spring_gamma(i);
        if (opt==0)
        {
           stream << "      bondings spring #" << Ifmt(5) << i+1 << std::endl;
           stream << "      spring parameters:" << std::endl;
           stream << "         coefficient index1 index2       distance" << std::endl;
           for (auto j=0; j<this->n0[i]; ++j)
              stream << Ffmt(20,6) << this->coef[i][j]  
                     << Ifmt(7) << this->indx[i][2*j]
                     << Ifmt(7) << this->indx[i][2*j+1]
                     << Ffmt(15,6) << this->ij_spring_distance(i,j) << std::endl;
           stream << "         K0         =" << Ffmt(12,6) << this->K0[i]     << std::endl;
           stream << "         gamma0     =" << Ffmt(12,6) << this->gamma0[i] << std::endl;
           stream << "      gamma         =" << Ffmt(12,6) << gamma   << std::endl;
           stream << "      spring energy =" << Ffmt(12,6) << espring << std::endl;
        }
        if (opt==1)
        {
           stream << " bondings spring #" << Ifmt(5) << i+1 << std::endl;
           stream << " spring parameters:" << std::endl;
           stream << "    coefficient index1 index2       distance" << std::endl;
           for (auto j=0; j<this->n0[i]; ++j)
              stream << Ffmt(15,6) << this->coef[i][j]  
                     << Ifmt(7) << this->indx[i][2*j]
                     << Ifmt(7) << this->indx[i][2*j+1]
                     << Ffmt(15,6) << this->ij_spring_distance(i,j) << std::endl;
           stream << "    K0         =" << Ffmt(12,6) << this->K0[i]     << std::endl;
           stream << "    gamma0     =" << Ffmt(12,6) << this->gamma0[i] << std::endl;
           stream << " gamma         =" << Ffmt(12,6) << gamma   << std::endl;
           stream << " spring energy =" << Ffmt(12,6) << espring << std::endl;
        }
        stream << std::endl;
     }
     return stream.str();
  }

} // namespace pwdft

