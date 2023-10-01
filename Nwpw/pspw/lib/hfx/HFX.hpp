#ifndef _HFX_HPP_
#define _HFX_HPP_

#pragma once

#include "Control2.hpp"
#include "Coulomb2.hpp"
#include "Pneb.hpp"

namespace pwdft {

class HFX_Operator {

   Pneb *mypneb;

   bool has_coulomb2 = false;
   bool new_coulomb2 = false;
   Coulomb2_Operator *mycoulomb2;

   int ispin;
   int norbs[2];
   int *orbital_list[2];
   double *vg;
   double *ehfx_orb[2];

public:
   bool hfx_on = false; 
   bool hfx_virtual_on = true;
   bool relaxed = true;
   bool orb_contribution = false;
   bool butterfly = false;
   bool replicated = false;
 
   int solver_type = 0;
   int screening_type = 0;
   double rcut = 8.0;
   double pp = 8.0;
   double hfx_parameter = 1.0;
   double ehfx = 0.0;
   double phfx = 0.0;
 
 
   /* Constructors */
   HFX_Operator(Pneb *, bool, Coulomb2_Operator *, Control2 &);
 
   /* destructor */
   ~HFX_Operator() { 
      if (hfx_on)
      {
         if (solver_type==0)
            delete[] vg; 

         if (new_coulomb2)
            delete mycoulomb2;

         for (auto ms=0; ms<ispin; ++ms)
         {
            delete[] ehfx_orb[ms];
            delete[] orbital_list[ms];
         }
      }
   }
 
   void v_exchange(const double *, double *);
   void e_exchange(const double *, double &, double &);
   // void   vcoulomb_dielec(const double *, double *);
   // void   vcoulomb_dielec2(const double *, double *, double *);

   friend std::ostream &operator<<(std::ostream &os, const HFX_Operator &hfx) {
      os << "";
      if (hfx.hfx_on)
      {
         if (hfx.relaxed)
            os << "    - HFX relaxed" << std::endl;
         else
            os << "    - HFX unrelaxed" << std::endl;
     
         if (hfx.ispin==1)
         {
            int stages = hfx.norbs[0] / 10;
            int extra  = hfx.norbs[0] % 10;
            for (auto s=0; s<stages; ++s)
            {
               os << "    - HFX restricted orbitals :" << Ifmt(5)
                  << " " << hfx.orbital_list[0][10*s+0]+1 << " " 
                  << " " << hfx.orbital_list[0][10*s+1]+1 << " "
                  << " " << hfx.orbital_list[0][10*s+2]+1 << " "
                  << " " << hfx.orbital_list[0][10*s+3]+1 << " "
                  << " " << hfx.orbital_list[0][10*s+4]+1 << " "
                  << " " << hfx.orbital_list[0][10*s+5]+1 << " "
                  << " " << hfx.orbital_list[0][10*s+6]+1 << " "
                  << " " << hfx.orbital_list[0][10*s+7]+1 << " "
                  << " " << hfx.orbital_list[0][10*s+8]+1 << " "
                  << " " << hfx.orbital_list[0][10*s+9]+1 << std::endl;
            }
            if (extra>0)
            {
               os << "    - HFX restricted orbitals :";
               for (auto s=0; s<extra; ++s)
                  os << " " << Ifmt(5) << hfx.orbital_list[0][10*stages+s]+1;
               os << std::endl;
            }
         }
         else
         {
            os << "    - HFX alpha orbitals :" << std::endl;
            os << "    - HFX beta orbitals  :" << std::endl;
         }
     
         if (hfx.solver_type==0)
         {
            os << "    - HFX screened coulomb solver (-periodic)" << std::endl;
            os << "    - HFX screening radius  (-screening_radius)  = " << Efmt(8,3) << hfx.rcut  << std::endl;
            os << "    - HFX screening power   (-screening_power)   = " << Efmt(8,3) << hfx.pp << std::endl;
            os << "    - HFX screening type    (-screening_type)    = " << hfx.screening_type << std::endl;
         }
         else
            os << "    - HFX free-space coulomb solver (-aperiodic)" << std::endl;
     
         if (hfx.hfx_parameter!=1.0) 
            os << "    - HFX scaling parameter (-scaling_parameter) = " << Efmt(8,3) << hfx.hfx_parameter << std::endl;
      }
      return os;
   }
};


} // namespace pwdft

#endif

