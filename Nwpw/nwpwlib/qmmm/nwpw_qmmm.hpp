#ifndef _NWPW_QMMM_HPP_
#define _NWPW_QMMM_HPP_

#pragma once

// ********************************************************************
// *                                                                  *
// *             nwpw_qmmm: top level object for qmmm                 *
// *                                                                  *
// *   The algorithms used in this module are based on the work of    *
// *   the nwchem-PSPW code developed by Eric J. Bylaska.             *
// *                                                                  *
// *   Citation:                                                      *
// *   Cauët E,Bogatko S,Weare JH,Fulton JL,Schenter GK, Bylaska EJ.  *
// *   Structure and dynamics of the hydration shells of the Zn(2+)   *
// *   ion from ab initio molecular dynamics and combined ab initio   *
// *   and classical molecular dynamics simulations. J Chem Phys.     * 
// *   2010 May 21;132(19):194502. doi: 10.1063/1.3421542.            *
// *                                                                  *
// *                                                                  *
// ********************************************************************
#include "blas.h"
#include <iostream>

#include "Control2.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Strfac.hpp"
//#include "gdevice.hpp"

namespace pwdft {

class nwpw_qmmm {

   Pneb *mypneb;
   Ion *myion;
   Strfac *mystrfac;
   double dv;
   int n2ft3d, ispin, *ne;
 
   //qmmm data
   int    *indx_frag_start, *size_frag, *kfrag;
   double *switch_Rin, *switch_Rout, *rcell;
   bool   *self_interaction, *incell_frag;
   int    nshl3d,nfrag,nkfrag;
   double qmmm_lmbda;
   bool   qmmm_found, shake_found, auxiliary_only, periodic, lmbda_flag;
 
   //shake qmmm data
   int *nshake, *na, *nb, *nab_shake_start, *indx_shake, *nindx_shake_start;
   int *ndsq_shake_start;
   double *dsq_shake;

   //bond spring qmmm data
      integer nbond(2)
      integer indx_bond(2),nindx_bond_start(2)
      integer Kr0_bond(2), nKr0_bond_start(2)
      logical bond_found
      common / pspw_qmmm3 /  nbond,
     >                       indx_bond,nindx_bond_start,
     >                       Kr0_bond,nKr0_bond_start,
     >                       bond_found

   //angle spring qmmm data
      integer nangle(2)
      integer indx_angle(2),nindx_angle_start(2)
      integer Kr0_angle(2), nKr0_angle_start(2)
      logical angle_found
      common / pspw_qmmm4 /  nangle,
     >                       indx_angle,nindx_angle_start,
     >                       Kr0_angle,nKr0_angle_start,
     >                       angle_found

   //cross bond spring qmmm data
      integer ncbond(2)
      integer indx_cbond(2),nindx_cbond_start(2)
      integer Kr0_cbond(2), nKr0_cbond_start(2)
      logical cbond_found
      common / pspw_qmmm5 /  ncbond,
     >                       indx_cbond,nindx_cbond_start,
     >                       Kr0_cbond,nKr0_cbond_start,
     >                       cbond_found


   //bond morse common qmmm data
      integer nmbond(2)
      integer indx_mbond(2),nindx_mbond_start(2)
      integer Kr0_mbond(2), nKr0_mbond_start(2)
      logical mbond_found
      common / pspw_qmmm6 /  nmbond,
     >                       indx_mbond,nindx_mbond_start,
     >                       Kr0_mbond,nKr0_mbond_start,
     >                       mbond_found

   //dihedral spring qmmm data
      integer ndihedral_spring(2)
      integer indx_dihedral_spring(2),nindx_dihedral_spring_start(2)
      integer Kr0_dihedral_spring(2), nKr0_dihedral_spring_start(2)
      logical dihedral_spring_found
      common / pspw_qmmm7 /  ndihedral_spring,
     >               indx_dihedral_spring,nindx_dihedral_spring_start,
     >               Kr0_dihedral_spring,nKr0_dihedral_spring_start,
     >               dihedral_spring_found

   //dihedral qmmm data
      integer ndihedral(2)
      integer indx_dihedral(2),nindx_dihedral_start(2)
      integer Kr0_dihedral(2), nKr0_dihedral_start(2)
      logical dihedral_found
      common / pspw_qmmm8 /  ndihedral,
     >               indx_dihedral,nindx_dihedral_start,
     >               Kr0_dihedral,nKr0_dihedral_start,
     >               dihedral_found

   //coordination number spring qmmm data
      integer ncoord(2)
      integer nindx0_coord_start(2)
      integer indx1_coord(2),nindx1_coord_start(2),nsize1_coord(2)
      integer indx2_coord(2),nindx2_coord_start(2),nsize2_coord(2)
      integer Kr0_coord(2), nKr0_coord_start(2)
      logical coord_found
      common / pspw_qmmm9 / ncoord,
     >                      nindx0_coord_start,
     >                      indx1_coord,nindx1_coord_start,nsize1_coord,
     >                      indx2_coord,nindx2_coord_start,nsize2_coord,
     >                      Kr0_coord,nKr0_coord_start,
     >                      coord_found


   //bond spring qmmm data
      integer nlink
      integer indx_link(2),param_link(2)
      logical link_found
      common / pspw_qmmm10 / nlink,
     >                       indx_link,param_link,
     >                       link_found

   //bondings spring qmmm data
      integer nbondings(2)
      integer nindx0_bondings_start(2)
      integer indx1_bondings(2),nindx1_bondings_start(2)
      integer nsize1_bondings(2)
      integer Kr0_bondings(2), nKr0_bondings_start(2)
      logical bondings_found
      common / pspw_qmmm11 / nbondings,
     >                       nindx0_bondings_start,
     >                       indx1_bondings,nindx1_bondings_start,
     >                       nsize1_bondings,
     >                       Kr0_bondings, nKr0_bondings_start,
     >                       bondings_found


public:
  double autoDebye = 2.5416;

  bool efield_on;
  int efield_type = 2; // default efield_type = rgrid
  double efield_vector[3] = {0.0, 0.0, 0.0};
  double efield_center[3] = {0.0, 0.0, 0.0};

  double *v_field;

  /* constructor */
  nwpw_qmmm(Ion *, Pneb *, Strfac *, Control2 &, std::ostream &);

  /* destructor */
  ~nwpw_qmmm() {
    if (efield_on)
      mypneb->r_dealloc(v_field);
  }

  // void gen_dipole(const double *);
  // void gen_Resta_dipole(const double *, double *);

  // double Qtot_APC(const int);
  // std::string print_APC(const double *);

  /**********************************
   *                                *
   *        efield_ion_energy       *
   *                                *
   **********************************/
  // calculates the energy between the QM ions and efield. Note the ions are
  // positive.
  double efield_ion_energy() {
    double eetmp = 0.0;
    if (efield_on) {
       for (auto ii=0; ii<myion->nion; ++ii) {
          double qii = myion->zv_psp[myion->katm[ii]];
          eetmp += qii*efield_vector[0]*(myion->rion1[3*ii]   - efield_center[0]) 
                 + qii*efield_vector[1]*(myion->rion1[3*ii+1] - efield_center[1])
                 + qii*efield_vector[2]*(myion->rion1[3*ii+2] - efield_center[2]);
       }
    }
    return eetmp;
  }

  /**********************************
   *                                *
   *         efield_ion_fion        *
   *                                *
   **********************************/
  // calculates the forces between the QM ions and efield. Note the ions are
  // positive.
  void efield_ion_fion(double *fion) {
     if (efield_on && ((efield_type==0)||(efield_type==2)) ) {
        for (auto ii=0; ii<myion->nion; ++ii) {
           double qii = myion->zv_psp[myion->katm[ii]];
           fion[3*ii]   -= qii*efield_vector[0];
           fion[3*ii+1] -= qii*efield_vector[1];
           fion[3*ii+2] -= qii*efield_vector[2];
        }
     }
  }

  /**********************************
   *                                *
   *         shortprint_efield      *
   *                                *
   **********************************/
  std::string shortprint_efield();
};

} // namespace pwdft

#endif
