#ifndef _MOLECULE_HPP_
#define _MOLECULE_HPP_

#pragma once

#include <iomanip>
#include <iostream>
#include <vector>

#include "Control2.hpp"
#include "Electron.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "psi.hpp"
#include "psi_H.hpp"

namespace pwdft {

#define ionstream(A, B, C)                                                     \
  (A) << std::scientific << std::setw(19) << std::setprecision(10) << (B)      \
      << std::setw(0) << " (" << std::setw(15) << std::setprecision(5) << (C)  \
      << std::setw(0) << " /ion)"                                              \
      << "\n"
#define elcstream(A, B, C)                                                     \
  (A) << std::scientific << std::setw(19) << std::setprecision(10) << (B)      \
      << std::setw(0) << " (" << std::setw(15) << std::setprecision(5) << (C)  \
      << std::setw(0) << " /electron)"                                         \
      << "\n"

#define eig1stream(A, B)                                                       \
  std::scientific << std::setw(18) << std::setprecision(7) << (A)              \
                  << std::setw(0) << " (" << std::fixed << std::setw(8)        \
                  << std::setprecision(3) << (B) << std::setw(0) << "eV)\n"
#define eig2stream(A, B, C, D)                                                 \
  std::scientific << std::setw(18) << std::setprecision(7) << (A)              \
                  << std::setw(0) << " (" << std::fixed << std::setw(8)        \
                  << std::setprecision(3) << (B) << std::setw(0) << "eV) "     \
                  << std::scientific << std::setw(18) << std::setprecision(7)  \
                  << (C) << std::setw(0) << " (" << std::fixed << std::setw(8) \
                  << std::setprecision(3) << (D) << std::setw(0) << " eV)\n"

/*
E[0] =
E[1] =
E[2] =
E[3] =
E[4] =
E[5] =
E[6] =
*/

class Molecule {

  double omega, scal2, scal1, dv;

  int ispin, ne[2], neall, n2ft3d, shift1, shift2;
  int nextra[2];
  int ne_excited[2] = {0,0};
  int nfft[3];
  int version = 3;

public:
   Pneb *mygrid;
   Ion *myion;
   Strfac *mystrfac;
   Ewald *myewald;
   Electron_Operators *myelectron;
   Pseudopotential *mypsp;
 
   double *psi1,*rho1,*rho1_all,*dng1;
   double *psi2,*rho2,*rho2_all,*dng2;
   double *lmbda,*hml,*eig;
   double *occ1 = nullptr;
   double *occ2 = nullptr;

   double *psi1_excited,*psi2_excited;
   double *hml_excited,*eig_excited;

   int multiplicity;
   double total_charge;

   // psi smearing block
   bool fractional=false;
   int smearoccupation, smeartype;
   double smearfermi[2], smearcorrection, smearkT;

   double E[80],en[2],ep,sp,tole;
 
   bool newpsi;
 
   /* Constructors */
   Molecule(char *,bool,Pneb *,Ion *,Strfac *,Ewald *,Electron_Operators *,Pseudopotential *, Control2 &, std::ostream &);
 
   /* destructor */
   ~Molecule() {
      delete[] psi1;
      delete[] rho1;
      delete[] rho1_all;
      delete[] dng1;
     
      delete[] psi2;
      delete[] rho2;
      delete[] rho2_all;
      delete[] dng2;
     
      delete[] lmbda;
      delete[] hml;
      delete[] eig;

      delete[] occ1; occ1 = nullptr;
      delete[] occ2; occ2 = nullptr;
 
      if ((ne_excited[0] + ne_excited[1])>0)
      {
         delete[] psi1_excited;
         delete[] psi2_excited;
         delete[] hml_excited;
         delete[] eig_excited;
      }
   }

   void epsi_initialize(char *,bool, const int *, std::ostream &);
   void epsi_finalize(char *, std::ostream &);
   void epsi_minimize(double *, std::ostream &);
   void epsi_get_gradient(double *, double *, double *);
   double epsi_KS_update_virtual(const int, const int, const int, const double, const double, double *, double *, double *, std::ostream &);

   void epsi_linesearch_update(double, double, double *, double *, double *, double *);
   void epsi_sort_virtual();

   void psi_minimize(double *, std::ostream &);
   void psi_get_gradient(double *, double *, double *);
   double psi_KS_update_orb(const int, const int, const int, const double, const double, double *, double *, double *, std::ostream &);


   double psi_KS_update(const int, const double, const double, double *, 
                        const int, const int *, double *,
                        double *, std::ostream &);


   void psi_linesearch_update(double, double, double *, double *, double *, double *);
   void psi_sort();
 
   /* write psi molecule */
   void writepsi(char *output_filename, std::ostream &coutput) {
      psi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne,
                psi1,output_filename,coutput);
   }

   void writepsi_excited(char *output_filename, std::ostream &coutput) {
      psi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne,
                psi1_excited,output_filename,coutput);
   }

 
   /* molecule energy */
   double energy() {
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
         
      E[0] = (myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy());
      return E[0];
   }
 
   double psi2_energy() {
      myelectron->run(psi2, rho2, dng2, rho2_all, occ2);
      E[0] = (myelectron->energy(psi2, rho2, dng2, rho2_all) + myewald->energy());
      return E[0];
   }
 
   /* molecule energy and eigenvalues */
   double energy_eigenvalues() {
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
      E[0] = (myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy());
      
      /* generate eigenvalues */
      myelectron->gen_hml(psi1, hml);
      mygrid->m_diagonalize(hml, eig);
      
      return E[0];
   }
 
   /* molecule energy and eigenvalues and other energies and en */
   double gen_all_energies() 
   {
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
      myelectron->gen_energies_en(psi1, rho1, dng1, rho1_all, E, en, occ1);
      
      /*  ion-ion energy */
      if (myelectron->is_periodic())
        E[4] = myewald->energy();
      if (myelectron->is_aperiodic())
        E[4] = myion->ion_ion_energy();
      
      E[0] += E[4];

      /* get contraints energies */
      if (myion->has_ion_bond_constraints()) {
         E[70] = myion->energy_ion_bond_constraints();
         E[0] += E[70];
      }
      if (myion->has_ion_bondings_constraints()) {
         E[71] = myion->energy_ion_bondings_constraints();
         E[0] +=  E[71];
      }

      
      /* generate eigenvalues */
      myelectron->gen_hml(psi1, hml);
      mygrid->m_diagonalize(hml, eig);

      if (fractional)
      {
         mygrid->m_0define_occupation(0.0, false,
                                     multiplicity,
                                     myion->total_zv(),total_charge,
                                     eig, hml, occ1,
                                     smeartype,smearkT,smearfermi,&smearcorrection);
        E[28] = smearcorrection;
        E[0]  +=  E[28];
      }

      
      /* generate dipole */
      mypsp->mydipole->gen_dipole(rho1);
      
      return E[0];
   }
 
   /* various molecule energies */
   double eorbit() { return myelectron->eorbit(psi1); }
   double psi2_eorbit() { return myelectron->eorbit(psi2); }
   double ehartree() { return myelectron->ehartree(dng1); }
   double exc() { return myelectron->exc(rho1_all); }
   double pxc() { return myelectron->pxc(rho1); }
   double eke() { return myelectron->eke(psi1,occ1); }
   double vl_ave() { return myelectron->vl_ave(dng1); }
   double vlr_ave() { return myelectron->vlr_ave(rho1); }
   double vnl_ave() { return myelectron->vnl_ave(psi1,occ1); }
   double eion() {
     double ee = 0.0;
     if (myelectron->is_periodic())
       ee = myewald->energy();
     if (myelectron->is_aperiodic())
       ee = myion->ion_ion_energy();
     return ee;
   }
 
   /* molecule - generate current hamiltonian */
   void gen_hml() { myelectron->gen_hml(psi1, hml); }

   void gen_vall() { myelectron->gen_vall(); }
   void get_vall(double *vall_out) { myelectron->get_vall(vall_out); }
   void set_vall(const double *vall_in) { myelectron->set_vall(vall_in); }
 
   /* molecule - diagonalize the current hamiltonian */
   void diagonalize() { mygrid->m_diagonalize(hml, eig); }
   void rotate1to2() { mygrid->fmf_Multiply(-1,psi1,hml,1.0,psi2,0.0); }
   
 
   /* molecule - call phafacs and gen_vl_potential and semicore */
   void phafacs_vl_potential_semicore() {
      mystrfac->phafac();
      myewald->phafac();
      myelectron->gen_vl_potential();
      myelectron->semicore_density_update();
   }
 
   /* apply psi2 = psi1 - dte*Hpsi1 + lmbda*psi1*/
   void sd_update(double dte) {
 
      /* apply psi2 = psi1 + dte*Hpsi1 */
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
      
      // myelectron->add_dteHpsi((-dte),psi1,psi2);
      myelectron->add_dteHpsi((dte), psi1, psi2);
      
      /* lagrange multiplier - Expensive */
      mygrid->ggm_lambda(dte, psi1, psi2, lmbda);
      
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }

   /* apply psi2 = psi1 - dte*Hpsi1 + lmbda*psi1*/
   void sd_update2(double dte) {

      /* apply psi2 = psi1 + dte*Hpsi1 */
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
      
      // myelectron->add_dteHpsi((-dte),psi1,psi2);
      myelectron->add_dteHpsi((dte), psi1, psi2);
      
      /* carry out direct ortho - Expensive */
      mygrid->g_ortho(-1,psi2);
      
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }
 
   /* apply psi2 = psi1 - dte*Hpsi1 + lmbda*psi1*/
   void sd_update_sic(double dte) {
 
      /* apply psi2 = psi1 + dte*Hpsi1 */
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
      
      // myelectron->add_dteHpsi((-dte),psi1,psi2);
      myelectron->add_dteHpsi((dte), psi1, psi2);
      
      /* lagrange multiplier - Expensive */
      mygrid->ggm_lambda_sic(dte, psi1, psi2, lmbda);
      
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }
 
   double psi_1get_Tgradient(double *G1) {
      double total_energy;
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
      total_energy = myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy();
      myelectron->gen_hml(psi1, hml);
      myelectron->get_Tgradient(psi1, hml, G1);
     
      return total_energy;
   }
 
   double psi_1get_TSgradient(double *G1) {
      double total_energy;
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
      total_energy =
          myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy();
      myelectron->gen_hmlt(psi1, hml);
      myelectron->get_Tgradient(psi1, hml, G1);
      
      return total_energy;
   }
 
   void swap_psi1_psi2() {
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }
 
   /* calculates the difference squared  error between rho1 and rho2 */
   double rho_error() {
      double x;
      double sumxx = 0.0;
      // mygrid->r_zero_ends(rho2);
      // mygrid->r_zero_ends(rho1);
      // if (ispin==2)
      //{
      //    mygrid->r_zero_ends(&rho1[n2ft3d]);
      //    mygrid->r_zero_ends(&rho2[n2ft3d]);
      // }
      for (int i = 0; i < n2ft3d; ++i) {
        x = (rho2[i] - rho1[i]);
        x += (rho2[i + (ispin - 1) * n2ft3d] - rho1[i + (ispin - 1) * n2ft3d]);
        sumxx += x * x;
      }
      return mygrid->d3db::parall->SumAll(1, sumxx) * dv;
   }
 
   std::vector<double> eig_vector() {
      return std::vector<double>(eig, &eig[neall]);
   }
 
   void psi_1local_force(double *grad_ion) {
      myelectron->vl_force(dng1, grad_ion);
      if (myelectron->is_aperiodic())
         myelectron->vlr_force(rho1, grad_ion);
   }
 
   void psi_1nonlocal_force(double *grad_ion) {
      myelectron->vnl_force(psi1, grad_ion);
   }
   void semicore_force(double *grad_ion) {
      myelectron->semicore_force(grad_ion);
   }
   void ewald_fion(double *grad_ion) { myewald->force(grad_ion); }
   void ion_fion(double *grad_ion) {
      if (myelectron->is_periodic())
         myewald->force(grad_ion);
      if (myelectron->is_aperiodic())
         myion->ion_ion_force(grad_ion);

      myion->add_contraint_force(grad_ion);
   }
 
   void psi_1apc_force(double *grad_ion) 
   {
      myelectron->apc_force(dng1, grad_ion);
   }

   void psi_1dielec_force(double *grad_ion) 
   {
      myelectron->dielectric_force(grad_ion);
   }

   std::string print_virtual() 
   {
      std::stringstream stream;
 
      std::ios init(NULL);
      init.copyfmt(stream);
      std::string eoln = "\n";

      stream << eoln;
      stream << eoln;
      stream << " virtual orbital energies:" << eoln;
      int nn = ne_excited[0] - ne_excited[1];
      double ev = 27.2116;

      // Print the first set of excited states in reverse order without symmetry considerations
      for (int i = ne_excited[0]-1; i>=ne_excited[1]; --i) 
         stream << eig1stream(eig_excited[i], eig_excited[i]*ev);

      // Print the second set of excited states in reverse order without symmetry considerations
      for (int i = ne_excited[1]-1; i>=0; --i) 
      {
         stream << eig2stream(
                     eig_excited[i + nn], eig_excited[i + nn]*ev,
                     eig_excited[i + (ispin - 1) * ne_excited[0]],
                     eig_excited[i + (ispin - 1) * ne_excited[0]]*ev);
      }
      //for (int i=0; i<nn; ++i)
      //   stream << eig1stream(eig_excited[i], eig_excited[i] * ev);
      //for (int i=0; i<ne_excited[1]; ++i)
      //   stream << eig2stream(eig_excited[i+nn], 
      //                        eig_excited[i+nn]*ev,
      //                        eig_excited[i+(ispin-1)*ne_excited[0]],
      //                        eig_excited[i+(ispin-1)*ne_excited[0]]*ev);

      stream << eoln;

      return stream.str();
   }
 
   friend std::ostream &operator<<(std::ostream &os, const Molecule &mymolecule) {
      /* using old style c++ formatting */
      std::ios init(NULL);
      init.copyfmt(os);
      std::string eoln = "\n";
      os << "     =============  energy results (Molecule object)  =============" << eoln;
      os << eoln << eoln;
     
      os << std::fixed << " number of electrons: spin up= " << std::setw(11)
         << std::setprecision(5) << mymolecule.en[0]
         << "  down= " << std::setw(11) << std::setprecision(5)
         << mymolecule.en[mymolecule.ispin - 1] << " (real space)";
      os << eoln << eoln;
      if (mymolecule.mypsp->myapc->v_apc_on)
         os << mymolecule.mypsp->myapc->shortprint_APC();
      os << eoln;
      os << ionstream(" total     energy    : ", mymolecule.E[0],mymolecule.E[0]/mymolecule.myion->nion);
      os << elcstream(" total orbital energy: ", mymolecule.E[1],mymolecule.E[1]/mymolecule.neall);
      os << elcstream(" hartree energy      : ", mymolecule.E[2],mymolecule.E[2]/mymolecule.neall);
      os << elcstream(" exc-corr energy     : ", mymolecule.E[3],mymolecule.E[3]/mymolecule.neall);

      if (mymolecule.myelectron->is_dielectric_on())
         os << " dielectric energy   : "
            << Efmt(19,10) << mymolecule.E[61] << " ("
            << Efmt(15,5)  << mymolecule.E[61]/(mymolecule.neall) << " /electron)" << std::endl;

      if (mymolecule.myelectron->is_v_apc_on())
         os << ionstream(" APC energy          : ", mymolecule.E[51],mymolecule.E[51]/mymolecule.myion->nion);

      os << ionstream(" ion-ion energy      : ", mymolecule.E[4], mymolecule.E[4]/mymolecule.myion->nion);

      if (mymolecule.fractional)
          os << elcstream(" smearing energy     : ", mymolecule.E[28],mymolecule.E[28]/mymolecule.neall);

      os << eoln;
      os << elcstream(" kinetic (planewave) : ", mymolecule.E[5],mymolecule.E[5]/mymolecule.neall);
      os << elcstream(" V_local (planewave) : ", mymolecule.E[6],mymolecule.E[6]/mymolecule.neall);
      os << elcstream(" V_nl    (planewave) : ", mymolecule.E[7],mymolecule.E[7]/mymolecule.neall);
      os << elcstream(" V_Coul  (planewave) : ", mymolecule.E[8],mymolecule.E[8]/mymolecule.neall);
      os << elcstream(" V_xc    (planewave) : ", mymolecule.E[9],mymolecule.E[9]/mymolecule.neall);

      if (mymolecule.myelectron->is_dielectric_on())
         os << " K.S. V_dielec energy: "
            << Efmt(19,10) << mymolecule.E[62] << " ("
            << Efmt(15,5)  << mymolecule.E[62]/(mymolecule.neall) << " /electron)" << std::endl;

      if (mymolecule.myelectron->is_v_apc_on())
         os << ionstream(" K.S. V_APC energy   : ",mymolecule.E[52],mymolecule.E[52]/mymolecule.myion->nion);
     
      os << " Viral Coefficient   : " << std::setw(19) << std::setprecision(10)
         << (mymolecule.E[9]+mymolecule.E[8]+mymolecule.E[7]+mymolecule.E[6])/mymolecule.E[5];
      os << std::endl;

      if (mymolecule.fractional)
      {
         os << std::endl;
         double ev = 27.2116;
         if (mymolecule.ispin==1)
            os << " Fermi energy        : "
               << Efmt(19,10) << mymolecule.smearfermi[0] << " (" << Ffmt(8,3) << mymolecule.smearfermi[0]*ev << " eV)" << std::endl;
         else
            os << "  Fermi energy = "
               << Efmt(19,10) << mymolecule.smearfermi[0] << " (" << Ffmt(8,3) << mymolecule.smearfermi[0]*ev << " eV)"
               << Efmt(19,10) << mymolecule.smearfermi[1] << " (" << Ffmt(8,3) << mymolecule.smearfermi[1]*ev << " eV)" << std::endl;
      }

      if (mymolecule.myion->has_ion_constraints())
      {
         os << std::endl;
         if (mymolecule.myion->has_ion_bond_constraints())
            os << " spring bond         : " << Efmt(19,10) << mymolecule.E[70] << " ("
                                             << Efmt(15,5)  << mymolecule.E[70]/mymolecule.myion->nion << " /ion)" << std::endl;
         if (mymolecule.myion->has_ion_bondings_constraints())
            os << " spring bondings     : " << Efmt(19,10) << mymolecule.E[71] << " ("
                                            << Efmt(15,5)  << mymolecule.E[71]/mymolecule.myion->nion << " /ion)" << std::endl;
      }
      os << eoln;
      os << eoln;
      os << " orbital energies:" << eoln;
      int nn = mymolecule.ne[0] - mymolecule.ne[1];
      double ev = 27.2116;
      for (int i=0; i<nn; ++i)
      {
         //os << eig1stream(mymolecule.eig[i], mymolecule.eig[i] * ev);
         if (mymolecule.fractional)
            os << Efmt(18,7) << mymolecule.eig[i] << " (" << Ffmt(8,3) << mymolecule.eig[i] * ev << "eV) occ="
               << Ffmt(5,3) << mymolecule.occ1[i] << std::endl;
         else
            os << Efmt(18,7) << mymolecule.eig[i] << " (" << Ffmt(8,3) << mymolecule.eig[i] * ev << "eV)" << std::endl;
      }
      for (int i=0; i<mymolecule.ne[1]; ++i)
      {
         //os << eig2stream(mymolecule.eig[i+nn], 
         //                 mymolecule.eig[i+nn]*ev,
         //                 mymolecule.eig[i+(mymolecule.ispin-1)*mymolecule.ne[0]],
         //                 mymolecule.eig[i+(mymolecule.ispin-1)*mymolecule.ne[0]]*ev);
         if (mymolecule.fractional)
            os << Efmt(18,7) << mymolecule.eig[i+nn]    << " ("
               << Ffmt(8,3)  << mymolecule.eig[i+nn]*ev << "eV)  occ="
               << Ffmt(5,3)  << mymolecule.occ1[i+nn]   << " "
               << Efmt(18,7) << mymolecule.eig[i+(mymolecule.ispin-1)*mymolecule.ne[0]]    << " ("
               << Ffmt(8,3)  << mymolecule.eig[i+(mymolecule.ispin-1)*mymolecule.ne[0]]*ev << "eV) occ="
               << Ffmt(5,3)  << mymolecule.occ1[i+(mymolecule.ispin-1)*mymolecule.ne[0]] << std::endl;
         else
            os << Efmt(18,7) << mymolecule.eig[i+nn] << " ("
               << Ffmt(8,3)  << mymolecule.eig[i + nn] * ev << "eV) "
               << Efmt(18,7) << mymolecule.eig[i+(mymolecule.ispin-1)*mymolecule.ne[0]] << " ("
               << Ffmt(8,3)  << mymolecule.eig[i+(mymolecule.ispin-1)*mymolecule.ne[0]]*ev << "eV)" << std::endl;
      }
      os << eoln;
     
      // write dipoles
      os << mymolecule.mypsp->mydipole->shortprint_dipole();
     
      os.copyfmt(init);
     
      return os;
   }

};

} // namespace pwdft

#endif
