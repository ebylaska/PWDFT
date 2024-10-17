#ifndef _SOLID_HPP_
#define _SOLID_HPP_

#pragma once

#include <iomanip>
#include <iostream>
#include <vector>

#include "Control2.hpp"
#include "cElectron.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"
#include "CStrfac.hpp"
#include "cpsi.hpp"

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

class Solid {

  double omega, scal2, scal1, dv;

  int ispin, ne[2], neall, nbrillq, nbrillouin, n2ft3d, nfft3d, shift1, shift2;
  int ne_excited[2] = {0,0};
  int nfft[3];
  int version = 3;

public:
   Cneb *mygrid;
   Ion *myion;
   CStrfac *mystrfac;
   Ewald *myewald;
   cElectron_Operators *myelectron;
   CPseudopotential *mypsp;
 
   double *psi1,*rho1,*rho1_all,*dng1;
   double *psi2,*rho2,*rho2_all,*dng2;
   double *lmbda,*hml,*eig;

   double *psi1_excited,*psi2_excited;
   double *hml_excited,*eig_excited;
 
   double E[80],en[2],ep,sp,tole;
 
   bool newpsi;
 
   /* Constructors */
   Solid(char *,bool,Cneb *,Ion *,CStrfac *,Ewald *,cElectron_Operators *,CPseudopotential *,Control2 &, std::ostream &);
 
   /* destructor */
   ~Solid() {
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
   }


   void epsi_initialize(char *,bool, const int *, std::ostream &);
   void epsi_finalize(char *, std::ostream &);
   void epsi_minimize(double *, std::ostream &);
   void epsi_get_gradient(const int, double *, double *, double *);
   double epsi_KS_update_virtual(const int, const int, const int, const int, const double, const double, double *, double *, double *, std::ostream &);

   void epsi_linesearch_update(const int, double, double, double *, double *, double *, double *);
   void epsi_sort_virtual(const int, double *, double *);

 
   /* write psi solid */
   void writepsi(char *output_filename, std::ostream &coutput) {
      cpsi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne,&nbrillouin,psi1,output_filename,coutput);
   }

   void writepsi_excited(char *output_filename, std::ostream &coutput) {
      cpsi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne,&nbrillouin,psi1_excited,output_filename,coutput);
   }

 
   /* solid energy */
   double energy() {
      myelectron->run(psi1, rho1, dng1, rho1_all);
      E[0] = (myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy());
      return E[0];
   }
 
   double psi2_energy() {
      myelectron->run(psi2, rho2, dng2, rho2_all);
      E[0] = (myelectron->energy(psi2, rho2, dng2, rho2_all) + myewald->energy());
      return E[0];
   }
 
   /* solid energy and eigenvalues */
   double energy_eigenvalues() {
      myelectron->run(psi1, rho1, dng1, rho1_all);
      E[0] = (myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy());
      
      /* generate eigenvalues */
      myelectron->gen_hml(psi1, hml);
      //mygrid->m_diagonalize(hml, eig);
      mygrid->w_diagonalize(hml, eig);
      
      return E[0];
   }
 
   /* solid energy and eigenvalues and other energies and en */
   double gen_all_energies() 
   {
      myelectron->run(psi1, rho1, dng1, rho1_all);
      myelectron->gen_energies_en(psi1, rho1, dng1, rho1_all, E, en);
      
      /*  ion-ion energy */
      E[4] = myewald->energy();
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
      //mygrid->m_diagonalize(hml, eig);
      mygrid->w_diagonalize(hml, eig);
      
      /* generate dipole */
      //mypsp->mydipole->gen_dipole(rho1);
      
      return E[0];
   }
 
   /* various solid energies */
   double eorbit() { return myelectron->eorbit(psi1); }
   double psi2_eorbit() { return myelectron->eorbit(psi2); }
   double ehartree() { return myelectron->ehartree(dng1); }
   double exc() { return myelectron->exc(rho1_all); }
   double pxc() { return myelectron->pxc(rho1); }
   double eke() { return myelectron->eke(psi1); }
   double vl_ave() { return myelectron->vl_ave(dng1); }
   double vnl_ave() { return myelectron->vnl_ave(psi1); }
   double eion() {
     double ee = 0.0;
     ee = myewald->energy();
     return ee;
   }
 
   /* solid - generate current hamiltonian */
   void gen_hml() { myelectron->gen_hml(psi1, hml); }

   void gen_vall() { myelectron->gen_vall(); }
   void get_vall(double *vall_out) { myelectron->get_vall(vall_out); }
   void set_vall(const double *vall_in) { myelectron->set_vall(vall_in); }
 
   /* solid - diagonalize the current hamiltonian */
   void diagonalize() { mygrid->w_diagonalize(hml, eig); }
 
   /* solid - call phafacs and gen_vl_potential and semicore */
   void phafacs_vl_potential_semicore() {
      mystrfac->phafac();
      myewald->phafac();
      myelectron->gen_vl_potential();
      myelectron->semicore_density_update();
   }
 
   /* apply psi2 = psi1 - dte*Hpsi1 + lmbda*psi1*/
   void sd_update(double dte) {
 
      /* apply psi2 = psi1 + dte*Hpsi1 */
      myelectron->run(psi1, rho1, dng1, rho1_all);
      
      // myelectron->add_dteHpsi((-dte),psi1,psi2);
      myelectron->add_dteHpsi((dte), psi1, psi2);
      
      /* lagrange multiplier - Expensive */
      mygrid->ggw_lambda(dte, psi1, psi2, lmbda);
      
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }
 
   /* apply psi2 = psi1 - dte*Hpsi1 + lmbda*psi1*/
   void sd_update_sic(double dte) {
 
      /* apply psi2 = psi1 + dte*Hpsi1 */
      myelectron->run(psi1, rho1, dng1, rho1_all);
      
      // myelectron->add_dteHpsi((-dte),psi1,psi2);
      myelectron->add_dteHpsi((dte), psi1, psi2);
      
      /* lagrange multiplier - Expensive */
      //mygrid->ggm_lambda_sic(dte, psi1, psi2, lmbda);
      
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }
 
   double psi_1get_Tgradient(double *G1) {
      double total_energy;
      myelectron->run(psi1, rho1, dng1, rho1_all);
      total_energy = myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy();
      myelectron->gen_hml(psi1, hml);
      myelectron->get_Tgradient(psi1, hml, G1);
     
      return total_energy;
   }
 
   double psi_1get_TSgradient(double *G1) {
      double total_energy;
      myelectron->run(psi1, rho1, dng1, rho1_all);
      total_energy = myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy();
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
   double rho_error() 
   {
      double x;
      double sumxx = 0.0;
      for (int i=0; i<nfft3d; ++i) 
      {
         x = (rho2[i] - rho1[i]);
         x += (rho2[i+(ispin-1)*nfft3d] - rho1[i+(ispin-1)*nfft3d]);
         sumxx += x * x;
      }
      return mygrid->c3db::parall->SumAll(1,sumxx) * dv;
   }
 
   std::vector<double> eig_vector() {
      return std::vector<double>(eig, &eig[neall]);
   }
 
   void psi_1local_force(double *grad_ion) {
      myelectron->vl_force(dng1, grad_ion);
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

      myion->add_contraint_force(grad_ion);
   }

   std::string print_virtual() 
   {
      std::stringstream stream;
 
      std::ios init(NULL);
      init.copyfmt(stream);
      std::string eoln = "\n";

      for (auto nbq=0; nbq<nbrillouin; ++nbq)
      {
         stream << eoln;
         stream << " Brillouin zone point: " << nbq + 1 << eoln << eoln;
         stream << " virtual orbital energies:" <<  eoln;
         int nn = ne_excited[0] - ne_excited[1];
         double ev = 27.2116;
      
         // Print the first set of excited states in reverse order without symmetry considerations
         for (int i = ne_excited[0]-1; i>=ne_excited[1]; --i) 
            stream << eig1stream(eig_excited[i + nbq*neall], eig_excited[i + nbq*neall]*ev);
      
         // Print the second set of excited states in reverse order without symmetry considerations
         for (int i = ne_excited[1]-1; i>=0; --i) 
         {
            stream << eig2stream(
                        eig_excited[i + nn + nbq*neall], eig_excited[i + nn + nbq*neall]*ev,
                        eig_excited[i + (ispin - 1) * ne_excited[0] + nbq*neall],
                        eig_excited[i + (ispin - 1) * ne_excited[0] + nbq*neall]*ev);
         }
         stream << eoln;
      }

      return stream.str();
   }


 

 
   friend std::ostream &operator<<(std::ostream &os, const Solid &mysolid) {
      /* using old style c++ formatting */
      std::ios init(NULL);
      init.copyfmt(os);
      std::string eoln = "\n";
      os << "     =============  energy results (Solid object)  =============" << eoln;
      os << eoln << eoln;
     
      os << std::fixed << " number of electrons: spin up= " << std::setw(11)
         << std::setprecision(5) << mysolid.en[0]
         << "  down= " << std::setw(11) << std::setprecision(5)
         << mysolid.en[mysolid.ispin - 1] << " (real space)";
      os << eoln << eoln;
      os << eoln;
      os << ionstream(" total     energy    : ", mysolid.E[0],mysolid.E[0]/mysolid.myion->nion);
      os << elcstream(" total orbital energy: ", mysolid.E[1],mysolid.E[1]/mysolid.neall);
      os << elcstream(" hartree energy      : ", mysolid.E[2],mysolid.E[2]/mysolid.neall);
      os << elcstream(" exc-corr energy     : ", mysolid.E[3],mysolid.E[3]/mysolid.neall);

      os << ionstream(" ion-ion energy      : ", mysolid.E[4], mysolid.E[4]/mysolid.myion->nion);
     
      os << eoln;
      os << elcstream(" kinetic (planewave) : ", mysolid.E[5],mysolid.E[5]/mysolid.neall);
      os << elcstream(" V_local (planewave) : ", mysolid.E[6],mysolid.E[6]/mysolid.neall);
      os << elcstream(" V_nl    (planewave) : ", mysolid.E[7],mysolid.E[7]/mysolid.neall);
      os << elcstream(" V_Coul  (planewave) : ", mysolid.E[8],mysolid.E[8]/mysolid.neall);
      os << elcstream(" V_xc    (planewave) : ", mysolid.E[9],mysolid.E[9]/mysolid.neall);


      //if (mysolid.myelectron->is_v_apc_on())
      //   os << ionstream(" K.S. V_APC energy   : ",mysolid.E[52],mysolid.E[52]/mysolid.myion->nion);
     
      os << " Viral Coefficient   : " << std::setw(19) << std::setprecision(10)
         << (mysolid.E[9]+mysolid.E[8]+mysolid.E[7]+mysolid.E[6])/mysolid.E[5];

      if (mysolid.myion->has_ion_constraints())
      {
         os << std::endl;
         if (mysolid.myion->has_ion_bond_constraints())
            os << " spring bond         : " << Efmt(19,10) << mysolid.E[70] << " ("
                                             << Efmt(15,5)  << mysolid.E[70]/mysolid.myion->nion << " /ion)" << std::endl;
         if (mysolid.myion->has_ion_bondings_constraints())
            os << " spring bondings     : " << Efmt(19,10) << mysolid.E[71] << " ("
                                            << Efmt(15,5)  << mysolid.E[71]/mysolid.myion->nion << " /ion)" << std::endl;
      }
      for (auto nb=0; nb<mysolid.nbrillouin; ++nb)
      {
         int nbq = mysolid.mygrid->ktoindex(nb);
         int pk = mysolid.mygrid->ktop(nb);
       
         int n = mysolid.ne[0] + mysolid.ne[1];
         double tmpeig[n];
         std::memset(tmpeig,0,n*sizeof(double));
         if (pk==mysolid.mygrid->c3db::parall->taskid_k())
            std::memcpy(tmpeig,mysolid.eig+nbq*n,n*sizeof(double));
         mysolid.mygrid->c3db::parall->Vector_SumAll(3,n,tmpeig);
         //std::memcpy(tmpeig,mysolid.eig+nbq*n,n*sizeof(double));

         os << eoln;
         os << eoln;
         os << mysolid.mygrid->mybrillouin->print_zone_point(nb);
         os << eoln;
         os << " orbital energies:" << eoln;
         int nn = mysolid.ne[0] - mysolid.ne[1];
         double ev = 27.2116;
         for (int i=0; i<nn; ++i)
            os << eig1stream(tmpeig[mysolid.ne[0]-1-i], tmpeig[mysolid.ne[0]-1-i] * ev);
         for (int i=0; i<mysolid.ne[1]; ++i)
            os << eig2stream(tmpeig[i+nn], 
                             tmpeig[i+nn]*ev,
                             tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]],
                             tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]]*ev);
         os << eoln;
      }
     
      // write dipoles
      //os << mysolid.mypsp->mydipole->shortprint_dipole();
     
      os.copyfmt(init);
     
      return os;
   }
};

} // namespace pwdft

#endif
