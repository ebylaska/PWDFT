
#include	"Control2.hpp"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Ewald.hpp"
#include        "Strfac.hpp"
#include	"Electron.hpp"
#include	"psi.hpp"

#include	"Molecule.hpp"

/********************************************
 *                                          *
 *           Molecule::Molecule             *
 *                                          *
 ********************************************/
Molecule::Molecule(Control2 *mycontrol0, Pneb *mygrid0, Ion *myion0, Strfac *mystrfac0, Ewald *myewald0, Electron_Operators *myelectron0)
{
   mycontrol  = mycontrol0;
   mygrid     = mygrid0;
   myion      = myion0;
   mystrfac   = mystrfac0;
   myewald    = myewald0;
   myelectron = myelectron0;

   ispin = mygrid->ispin;
   neall = mygrid->neq[0] + mygrid->neq[1];

   psi1  = mygrid->g_allocate(1);
   psi2  = mygrid->g_allocate(1);
   rho1  = mygrid->r_nalloc(ispin);
   rho2  = mygrid->r_nalloc(ispin);
   rho1_all  = mygrid->r_nalloc(ispin);
   rho2_all  = mygrid->r_nalloc(ispin);

   hml   = mygrid->m_allocate(-1,1);
   eig   = new double[mygrid->ne[0]+mygrid->ne[1]];

   omega = mygrid->lattice->omega();
   scal1 = 1.0/((double) ((mygrid->nx)*(mygrid->ny)*(mygrid->nz)));
   scal2 = 1.0/omega;
   dv = omega*scal1;

   n2ft3d = (mygrid->n2ft3d);
   shift1 = 2*(mygrid->npack(1));
   shift2 = (mygrid->n2ft3d);

   psi_read(mygrid,mycontrol->input_movecs_filename(),psi2);
}

