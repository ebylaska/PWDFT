#ifndef	_MOLECULE_HPP_
#define _MOLECULE_HPP_


using namespace std;

#include	"Control2.hpp"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Ewald.hpp"
#include        "Strfac.hpp"
#include	"Pseudopotential.hpp"
#include	"Electron.hpp"
#include	"psi.hpp"

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

   double *psi1, *rho1, *rho1_all, *dng1;
   double *psi2, *rho2, *rho2_all, *dng2;
   double *hml, *eig;

   /* Constructors */
   Molecule(char *,  Pneb *, Ion *, Strfac *, Ewald *, Electron_Operators *);

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
      return (myelectron->energy(psi1,rho1,dng1,rho1_all) +  myewald->energy());
   }


   /* molecule energy and eigenvalues */
   double energy_eigenvalues() {
      myelectron->run(psi1,rho1,dng1,rho1_all);
      double total_energy = (myelectron->energy(psi1,rho1,dng1,rho1_all) +  myewald->energy());

      myelectron->gen_hml(psi1,hml);
      mygrid->m_diagonalize(hml,eig);

      return total_energy;
   }

   /* various molecule energies */
   double eorbit()   { return myelectron->eorbit(psi1); }
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
  
   

};

#endif
