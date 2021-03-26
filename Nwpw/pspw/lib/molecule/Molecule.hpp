#ifndef	_MOLECULE_HPP_
#define _MOLECULE_HPP_


using namespace std;

#include	"Control2.hpp"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Ewald.hpp"
#include        "Strfac.hpp"
#include	"Pseudopotential.hpp"

class	Molecule {

   Control2  *mycontrol;
   Pneb   *mygrid;
   Ion    *myion;
   Strfac *mystrfac;
   Ewald  *myewald;
   Electron_Operators *myelectron;

   double *psi1, *rho1, *rho1_all;
   double *psi2, *rho2, *rho2_all;
   double *hml, *eig;

   double omega,scal2,scal1,dv;

   int ispin,neall,n2ft3d,shift1,shift2;

public:

   /* Constructors */
   Molecule(Control2 *, Pneb *, Ion *, Strfac *, Ewald *, Electron_Operators *);

   /* destructor */
   ~Molecule() {
         delete [] psi1;
         delete [] rho1;
         delete [] rho1_all;

         delete [] psi2;
         delete [] rho2;
         delete [] rho2_all;

         delete [] hml;
         delete [] eig;
    }

};

#endif
