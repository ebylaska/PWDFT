#ifndef _nwpw_lmbfgs_HPP_
#define _nwpw_lmbfgs_HPP_

#pragma once

/* nwpw_lmbfgs.hpp - this class is used to keep perform lmfgs steps
   Author - Eric Bylaska
*/

#include "blas.h"
#include <cmath>

namespace pwdft {

class nwpw_lmbfgs {

   int one = 1;
   double rone = 1.0;
   double mrone = -1.0;
 
   int m, max_m, nsize;
   int *indx;
   double *rho;
   double *sylist;

public:
   /* constructor */
   nwpw_lmbfgs(int nsize0, int max_m0, double *x0, double *g0) 
   {
      m = 0;
      max_m = max_m0;
      nsize = nsize0;
     
      indx = new int[max_m];
      rho  = new double[max_m];
      sylist = new double[2*(max_m)*nsize];
     
      for (auto k=0; k<max_m; ++k)
         indx[k] = k;
     
      // DCOPY_PWDFT(nsize,x0,one,&sylist[2*m*nsize],    one);
      // DCOPY_PWDFT(nsize,g0,one,&sylist[(2*m+1)*nsize],one);
      std::memcpy(sylist + (2*m)  *nsize, x0, nsize*sizeof(double));
      std::memcpy(sylist + (2*m+1)*nsize, g0, nsize*sizeof(double));
   }

   /* destructor */
   ~nwpw_lmbfgs() 
   {
      delete[] sylist;
      delete[] rho;
      delete[] indx;
   }

   void lmbfgs(double *x, double *g, double *q) 
   {
      // DCOPY_PWDFT(nsize,g,one,q,one);
      // DCOPY_PWDFT(nsize,x,one,&sylist[(2*(indx[m+1]))  *nsize],one);
      // DCOPY_PWDFT(nsize,g,one,&sylist[(2*(indx[m+1])+1)*nsize],one);
      std::memcpy(q, g, nsize*sizeof(double));
      std::memcpy(sylist+(2*(indx[m+1]))  *nsize, x, nsize*sizeof(double));
      std::memcpy(sylist+(2*(indx[m+1])+1)*nsize, g, nsize*sizeof(double));
     
      DAXPY_PWDFT(nsize, mrone, x, one, &sylist[(2*indx[m])  * nsize], one);
      DAXPY_PWDFT(nsize, mrone, g, one, &sylist[(2*indx[m]+1)* nsize], one);
     
      double threshold_factor = 1.0e-12; // Adjust this based on how "strict" you want to be
      double safety_floor = 1.0e-16;
      double sum = DDOT_PWDFT(nsize, &sylist[(2*indx[m]+1)*nsize], one, &sylist[(2*indx[m])*nsize], one);
      double yabs = 0.0;
      double sabs = 0.0;
      for (int i=0; i<nsize; ++i)
      {
         yabs += std::abs(sylist[(2*indx[m])*nsize+i]);
         sabs += std::abs(sylist[(2*indx[m]+1)*nsize+i]);
      }

      //if (std::fabs(sum) > 1.0e-11) 
      if (std::fabs(sum) > (threshold_factor * yabs * sabs + safety_floor))
      {
         rho[indx[m]] = 1.0 / sum;
        
         double alpha[m + 1];
         for (auto k=m-1; k>=0; --k) 
         {
            alpha[k] = rho[indx[k]] * DDOT_PWDFT(nsize, &sylist[(2*indx[k])*nsize], one, q, one); //dx terms
            double tscal = -alpha[k];
            DAXPY_PWDFT(nsize, tscal, &sylist[(2*indx[k]+1)*nsize], one, q, one);
         }
        
         /* add preconditioner here */
         if (m > 200) 
         {
            double sumsy = DDOT_PWDFT(nsize, &sylist[(2*indx[m-1])  *nsize], one, &sylist[(2*indx[m-1]+1)*nsize], one);
            double sumyy = DDOT_PWDFT(nsize, &sylist[(2*indx[m-1]+1)*nsize], one, &sylist[(2*indx[m-1]+1)*nsize], one);
            if (sumyy > 1.0e-10) 
            {
               double gmma = sumsy / sumyy;
               if (gmma > 2.0)
                  gmma = 2.0;
               DSCAL_PWDFT(nsize, gmma, q, one);
            }
         }
        
         for (auto k = 0; k < m; ++k) 
         {
            double beta = rho[indx[k]] * DDOT_PWDFT(nsize, &sylist[(2*indx[k]+1)*nsize], one, q, one);
            double tscal = -(beta - alpha[k]);
            DAXPY_PWDFT(nsize, tscal, &sylist[(2*indx[k])*nsize], one, q, one); //dx terms
         }
        
         if (m < (max_m-2)) 
         {
            ++m;
         } 
         else
         {
            int itmp = indx[0];
            for (auto k=0; k<(max_m-1); ++k)
               indx[k] = indx[k+1];
            indx[max_m-1] = itmp;
         }
      }
   }

   void lmbfgs_cartesian(double *x, double *g, double *q, const double *unita) 
   {
      std::memcpy(q, g, nsize*sizeof(double));
      std::memcpy(sylist+(2*(indx[m+1]))  *nsize, x, nsize*sizeof(double));
      std::memcpy(sylist+(2*(indx[m+1])+1)*nsize, g, nsize*sizeof(double));

      double* historical_x = &sylist[(2*indx[m])*nsize];

      // ========================================================================
      // 1. COMPUTE THE MATRIC INVERSE OF UNITA ON-THE-FLY (Cramer's Rule)
      // ========================================================================
      // unita layout: [0]=a1.x, [1]=a1.y, [2]=a1.z, [3]=a2.x, [4]=a2.y ...

      double ub[9];
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


      // ========================================================================
      // 2. APPLY MINIMUM IMAGE FILTER VIA PROJECTION SPACE
      // ========================================================================
      //DAXPY_PWDFT(nsize, mrone, x, one, &sylist[(2*indx[m])*nsize], one); 
      const int n_atoms = nsize / 3;
      std::cout << "NSIZE=" << nsize << std::endl;
      std::cout << "NATOMs=" << n_atoms << std::endl;
     
      for (int ii=0; ii<n_atoms; ++ii)
      {
         // Extract raw Cartesian displacement components for this specific atom
         double dx = historical_x[3*ii]   - x[3*ii];
         double dy = historical_x[3*ii+1] - x[3*ii+1];
         double dz = historical_x[3*ii+2] - x[3*ii+2];

         // Project the Cartesian displacement step vector into Fractional space
         double s1 = ub[0]*dx + ub[1]*dy + ub[2]*dz;
         double s2 = ub[3]*dx + ub[4]*dy + ub[5]*dz;
         double s3 = ub[6]*dx + ub[7]*dy + ub[8]*dz;

         // Apply Minimum Image Convention rounding shift in fractional coordinates
         //c1 = x*ub(1,1) + y*ub(2,1) + z*ub(3,1)
         //c2 = x*ub(1,2) + y*ub(2,2) + z*ub(3,2)
         //c3 = x*ub(1,3) + y*ub(2,3) + z*ub(3,3)
         //c1 = c1 - DNINT(c1)
         //c2 = c2 - DNINT(c2)
         //c3 = c3 - DNINT(c3)
         //x = ua(1,1)*c1 + ua(1,2)*c2 + ua(1,3)*c3
         //y = ua(2,1)*c1 + ua(2,2)*c2 + ua(2,3)*c3
         //z = ua(3,1)*c1 + ua(3,2)*c2 + ua(3,3)*c3

         s1 -= std::round(s1);
         s2 -= std::round(s2);
         s3 -= std::round(s3);

         // Transform back into clean, phase-aligned Cartesian displacements
         historical_x[3*ii]   = unita[0]*s1 + unita[3]*s2 + unita[6]*s3;
         historical_x[3*ii+1] = unita[1]*s1 + unita[4]*s2 + unita[7]*s3;
         historical_x[3*ii+2] = unita[2]*s1 + unita[5]*s2 + unita[8]*s3;
      }

      DAXPY_PWDFT(nsize, mrone, g, one, &sylist[(2*indx[m]+1)*nsize], one);

      double threshold_factor = 1.0e-12; // Adjust this based on how "strict" you want to be
      double safety_floor = 1.0e-16;
      double sum = DDOT_PWDFT(nsize, &sylist[(2*indx[m]+1)*nsize], one, &sylist[(2*indx[m])*nsize], one);
      double yabs = 0.0;
      double sabs = 0.0;
      for (int i=0; i<nsize; ++i)
      {
         yabs += std::abs(sylist[(2*indx[m])*nsize+i]);
         sabs += std::abs(sylist[(2*indx[m]+1)*nsize+i]);
      }

      //if (std::fabs(sum) > 1.0e-11)
      if (std::fabs(sum) > (threshold_factor * yabs * sabs + safety_floor))
      {
         rho[indx[m]] = 1.0/sum;
         std::cout << "m=" << m << " indx_m=" << indx[m] << std::scientific << std::setprecision(14) << " rho="  << rho[indx[m]] << std::endl;;
         
         
         double alpha[m+1];
         for (auto k=m-1; k>=0; --k) 
         {         
            //int k = indx[i];
            //double jj =  DDOT_PWDFT(nsize, &sylist[(2*indx[k])*nsize], one, q, one);
            //std::cout << "k=" << k << " indx_k=" << indx[k] << std::scientific << std::setprecision(14) << " jj=" << jj << " rho="  << rho[indx[k]] << std::endl;;

            alpha[k] = rho[indx[k]]*DDOT_PWDFT(nsize, &sylist[(2*indx[k])*nsize], one, q, one); //dx terms
            double tscal = -alpha[k];
            std::cout << "k=" << k << " indx_k=" << indx[k] << " alpha=" << alpha[k] << " tscal=" << tscal << std::endl;
            if (std::abs(tscal) < 1000.0)
            {
               DAXPY_PWDFT(nsize, tscal, &sylist[(2*indx[k]+1)*nsize], one, q, one);
               //std::cout << "include" << std::endl;
            }
         }

         /* add preconditioner here */
         if (m > 200) 
         {
            double sumsy = DDOT_PWDFT(nsize, &sylist[(2*indx[m-1])  *nsize], one, &sylist[(2*indx[m-1]+1)*nsize], one);
            double sumyy = DDOT_PWDFT(nsize, &sylist[(2*indx[m-1]+1)*nsize], one, &sylist[(2*indx[m-1]+1)*nsize], one);
            if (sumyy > 1.0e-10)
            {
               double gmma = sumsy / sumyy;
               if (gmma > 2.0)
                  gmma = 2.0;
               DSCAL_PWDFT(nsize, gmma, q, one);
            }
         }

         for (auto k = 0; k < m; ++k)
         {
            //int k = indx[i];
            double beta = rho[indx[k]] * DDOT_PWDFT(nsize, &sylist[(2*indx[k]+1)*nsize], one, q, one);
            double tscal = -(beta - alpha[k]);
            std::cout << "k=" << k << " indx_k=" << indx[k] << " beta=" << beta << " tscal=" << tscal << std::endl;
            if (std::abs(tscal) < 1000.0)
            {
               DAXPY_PWDFT(nsize, tscal, &sylist[(2*indx[k])*nsize], one, q, one); //dx terms
               //std::cout << "include2" << std::endl;
            }
         }

         if (m < (max_m-2))
         {
            ++m;
         }
         else
         {
            int itmp = indx[0];
            for (auto k=0; k<(max_m-1); ++k)
               indx[k] = indx[k+1];
            indx[max_m-1] = itmp;
         }
      }



   }
};

} // namespace pwdft

#endif
