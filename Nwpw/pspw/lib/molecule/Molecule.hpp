#ifndef	_MOLECULE_HPP_
#define _MOLECULE_HPP_


#include        <iostream>
#include 	<iomanip> 
#include	<vector>
using namespace std;

#include	"Control2.hpp"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Ewald.hpp"
#include        "Strfac.hpp"
#include	"Pseudopotential.hpp"
#include	"Electron.hpp"
#include	"psi.hpp"

namespace pwdft {

#define	ionstream(A,B,C)	(A) << scientific << setw(19) << setprecision(10) << (B) << setw(0) <<  " (" << setw(15) << setprecision(5)  << (C) << setw(0) << " /ion)" << "\n"
#define	elcstream(A,B,C)	(A) << scientific << setw(19) << setprecision(10) << (B) << setw(0) <<  " (" << setw(15) << setprecision(5)  << (C) << setw(0) << " /electron)" << "\n"

#define	eig1stream(A,B)		scientific << setw(18) << setprecision(7) << (A) << setw(0) <<  " (" << fixed << setw(8) << setprecision(3)  << (B) << setw(0) << "eV)\n"
#define	eig2stream(A,B,C,D)	scientific << setw(18) << setprecision(7) << (A) << setw(0) <<  " (" << fixed << setw(8) << setprecision(3)  << (B) << setw(0) << "eV) " << scientific << setw(18) << setprecision(7) << (C) << setw(0) <<  " (" << fixed << setw(8) << setprecision(3)  << (D) << setw(0) << " eV)\n" 

class	Molecule {

   double omega,scal2,scal1,dv;

   int ispin,ne[2],neall,n2ft3d,shift1,shift2;
   int nfft[3];
   int version=3;
   

public:
   Pneb   *mygrid;
   Ion    *myion;
   Strfac *mystrfac;
   Ewald  *myewald;
   Electron_Operators *myelectron;
   Pseudopotential  *mypsp;



   double *psi1, *rho1, *rho1_all, *dng1;
   double *psi2, *rho2, *rho2_all, *dng2;
   double *lmbda, *hml, *eig;

   double E[60],en[2];

   bool newpsi;

   /* Constructors */
   Molecule(char *,  bool, Pneb *, Ion *, Strfac *, Ewald *, Electron_Operators *, Pseudopotential *);

   /* destructor */
   ~Molecule() {
         delete [] psi1;
         delete [] rho1;
         delete [] rho1_all;
         delete [] dng1;

         delete [] psi2;
         delete [] rho2;
         delete [] rho2_all;
         delete [] dng2;

         delete [] lmbda;
         delete [] hml;
         delete [] eig;
    }


   /* write psi molecule */
   void writepsi(char *output_filename) {
      psi_write(mygrid, &version, nfft, mygrid->lattice->unita_ptr(),
                &ispin, ne, psi1, output_filename);
   }


   /* molecule energy */
   double energy() {
      myelectron->run(psi1,rho1,dng1,rho1_all);
      E[0] = (myelectron->energy(psi1,rho1,dng1,rho1_all) +  myewald->energy());
      return E[0];
   }

   double psi2_energy() {
      myelectron->run(psi2,rho2,dng2,rho2_all);
      E[0] = (myelectron->energy(psi2,rho2,dng2,rho2_all) +  myewald->energy());
      return E[0];
   }


   /* molecule energy and eigenvalues */
   double energy_eigenvalues() {
      myelectron->run(psi1,rho1,dng1,rho1_all);
      E[0] = (myelectron->energy(psi1,rho1,dng1,rho1_all) +  myewald->energy());

      /* generate eigenvalues */
      myelectron->gen_hml(psi1,hml);
      mygrid->m_diagonalize(hml,eig);

      return E[0];
   }


   /* molecule energy and eigenvalues and other energies and en */
   double gen_all_energies() {
      myelectron->run(psi1,rho1,dng1,rho1_all);
      myelectron->gen_energies_en(psi1,rho1,dng1,rho1_all,E,en);
      E[4] = myewald->energy();
      E[0] += E[4];

      /* generate eigenvalues */
      myelectron->gen_hml(psi1,hml);
      mygrid->m_diagonalize(hml,eig);

      return E[0];
   }

   /* various molecule energies */
   double eorbit()   { return myelectron->eorbit(psi1); }
   double psi2_eorbit()  { return myelectron->eorbit(psi2); }
   double ehartree() { return myelectron->ehartree(dng1); }
   double exc()      { return myelectron->exc(rho1_all); }
   double pxc()      { return myelectron->pxc(rho1); }
   double eke()      { return myelectron->eke(psi1); }
   double vl_ave()   { return myelectron->vl_ave(dng1); }
   double vnl_ave()  { return myelectron->vnl_ave(psi1); }
   double eion()     { return myewald->energy(); }

   /* molecule - generate current hamiltonian */
   void gen_hml() { myelectron->gen_hml(psi1,hml); }

   /* molecule - diagonalize the current hamiltonian */
   void diagonalize() { mygrid->m_diagonalize(hml,eig); }

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
      myelectron->run(psi1,rho1,dng1,rho1_all);
      myelectron->add_dteHpsi((-dte),psi1,psi2);

      /* lagrange multiplier - Expensive */
      mygrid->ggm_lambda(dte,psi1,psi2,lmbda);

      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }

   double psi_1get_Tgradient(double *G1)
   {
      double total_energy;
      myelectron->run(psi1,rho1,dng1,rho1_all);
      total_energy = myelectron->energy(psi1,rho1,dng1,rho1_all) + myewald->energy();
      myelectron->gen_hml(psi1,hml);
      myelectron->get_Tgradient(psi1,hml,G1);

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
      //mygrid->r_zero_ends(rho2);
      //mygrid->r_zero_ends(rho1);
      //if (ispin==2)
      //{
      //   mygrid->r_zero_ends(&rho1[n2ft3d]);
      //   mygrid->r_zero_ends(&rho2[n2ft3d]);
      //}
      for (int i=0; i<n2ft3d; ++i)
      {
         x  = (rho2[i]-rho1[i]);
         x += (rho2[i+(ispin-1)*n2ft3d]-rho1[i+(ispin-1)*n2ft3d]);
         sumxx += x*x;
      }
      return mygrid->d3db::parall->SumAll(1,sumxx)*dv;
   }

   std::vector<double> eig_vector() { return std::vector<double>(eig,&eig[neall]); }


   void psi_1local_force(double *grad_ion)    { myelectron->vl_force(dng1,grad_ion); }
   void psi_1nonlocal_force(double *grad_ion) { myelectron->vnl_force(psi1,grad_ion); }
   void semicore_force(double *grad_ion)      { myelectron->semicore_force(grad_ion); }
   void ewald_fion(double *grad_ion)          { myewald->force(grad_ion); }

   void psi_1apc_force(double *grad_ion)      { myelectron->apc_force(dng1,grad_ion); }

  

   friend ostream& operator<<(ostream& os, const Molecule& mymolecule) {
      /* using old style c++ formatting */
      ios init(NULL);
      init.copyfmt(os);
      string eoln = "\n";
      os << "     ==============  energy results (Molecule object)  ==============" << eoln;
      os << eoln << eoln;

      os << fixed << " number of electrons: spin up= " << setw(11) << setprecision(5) << mymolecule.en[0] 
                                         << "  down= " << setw(11) << setprecision(5) << mymolecule.en[mymolecule.ispin] 
                                         << " (real space)";
      os << eoln << eoln;
      if (mymolecule.mypsp->myapc->v_apc_on)
         os << mymolecule.mypsp->myapc->shortprint_APC();
      os << eoln;
      os << ionstream(" total     energy    : ",mymolecule.E[0],mymolecule.E[0]/mymolecule.myion->nion);
      os << elcstream(" total orbital energy: ",mymolecule.E[1],mymolecule.E[1]/mymolecule.neall);
      os << elcstream(" hartree energy      : ",mymolecule.E[2],mymolecule.E[2]/mymolecule.neall);
      os << elcstream(" exc-corr energy     : ",mymolecule.E[3],mymolecule.E[3]/mymolecule.neall);
      if (mymolecule.myelectron->is_v_apc_on())
         os << ionstream(" APC energy          : ",mymolecule.E[51],mymolecule.E[51]/mymolecule.myion->nion);
      os << ionstream(" ion-ion energy      : ",mymolecule.E[4],mymolecule.E[4]/mymolecule.myion->nion);

      os << eoln;
      os << elcstream(" kinetic (planewave) : ",mymolecule.E[5],mymolecule.E[5]/mymolecule.neall);
      os << elcstream(" V_local (planewave) : ",mymolecule.E[6],mymolecule.E[6]/mymolecule.neall);
      os << elcstream(" V_nl    (planewave) : ",mymolecule.E[7],mymolecule.E[7]/mymolecule.neall);
      os << elcstream(" V_Coul  (planewave) : ",mymolecule.E[8],mymolecule.E[8]/mymolecule.neall);
      os << elcstream(" V_xc    (planewave) : ",mymolecule.E[9],mymolecule.E[9]/mymolecule.neall);
      if (mymolecule.myelectron->is_v_apc_on())
         os << ionstream(" K.S. V_APC energy   : ",mymolecule.E[52],mymolecule.E[52]/mymolecule.myion->nion);

      os << " Viral Coefficient   : "
         << setw(19) << setprecision(10) << (mymolecule.E[9] + mymolecule.E[8] + mymolecule.E[7] + mymolecule.E[6])/mymolecule.E[5];
      os << eoln;
      os << eoln;
      os << " orbital energy:" << eoln;
      int nn = mymolecule.ne[0] - mymolecule.ne[1];
      double ev = 27.2116;
      for (int i=0; i<nn; ++i)
         os << eig1stream(mymolecule.eig[i],mymolecule.eig[i]*ev);
      for (int i=0; i<mymolecule.ne[1]; ++i)
         os << eig2stream(mymolecule.eig[i+nn], 
                          mymolecule.eig[i+nn]*ev,
                          mymolecule.eig[i+(mymolecule.ispin-1)*mymolecule.ne[0]], 
                          mymolecule.eig[i+(mymolecule.ispin-1)*mymolecule.ne[0]]*ev);
      os << eoln;


      os.copyfmt(init);

      return os;
   }
   

};

}

#endif
