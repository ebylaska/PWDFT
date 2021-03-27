
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
Molecule::Molecule(char *infilename, 
                   Pneb *mygrid0, Ion *myion0, Strfac *mystrfac0, Ewald *myewald0, Electron_Operators *myelectron0)
{
   mygrid     = mygrid0;
   myion      = myion0;
   mystrfac   = mystrfac0;
   myewald    = myewald0;
   myelectron = myelectron0;

   ispin = mygrid->ispin;
   neall = mygrid->neq[0] + mygrid->neq[1];
   ne[0] = mygrid->ne[0]; 
   ne[1] = mygrid->ne[1];
   nfft[0] = mygrid->nx; nfft[1] = mygrid->ny; nfft[2] = mygrid->nz;

   psi1  = mygrid->g_allocate(1);
   psi2  = mygrid->g_allocate(1);
   rho1  = mygrid->r_nalloc(ispin);
   rho2  = mygrid->r_nalloc(ispin);
   rho1_all  = mygrid->r_nalloc(ispin);
   rho2_all  = mygrid->r_nalloc(ispin);
   dng1 = mygrid->c_pack_allocate(0);
   dng2 = mygrid->c_pack_allocate(0);

   hml   = mygrid->m_allocate(-1,1);
   eig   = new double[mygrid->ne[0]+mygrid->ne[1]];

   omega = mygrid->lattice->omega();
   scal1 = 1.0/((double) ((mygrid->nx)*(mygrid->ny)*(mygrid->nz)));
   scal2 = 1.0/omega;
   dv = omega*scal1;

   n2ft3d = (mygrid->n2ft3d);
   shift1 = 2*(mygrid->npack(1));
   shift2 = (mygrid->n2ft3d);

   psi_read(mygrid,infilename,psi1);

   myelectron->gen_vl_potential();

/*---------------------- testing Electron Operators ---------------------- */
   double sum1;
   std::cout << "Molecule: Testing Electron Operators " << std::endl;
   //myelectron->gen_psi_r(psi1);
   //myelectron->gen_density(rho1);


   //double sum1 = mygrid->r_dsum(rho1)*dv;
   //std::cout << "integrate dn dv = " << sum1 << " dv=" << dv << std::endl;

   myelectron->run(psi1,rho1,dng1,rho1_all);
   sum1 = myelectron->energy(psi1,rho1,dng1,rho1_all);
   myelectron->gen_hml(psi1,hml);
   std::cout << "electron energy = " << sum1 <<  std::endl;
   std::cout << "electron eorbit = " << myelectron->eorbit(psi1) <<  std::endl;
   std::cout << "electron ehartr = " << myelectron->ehartree(dng1) <<  std::endl;
   std::cout << "electron exc    = " << myelectron->exc(rho1_all) <<  std::endl;
   std::cout << "electron pxc    = " << myelectron->pxc(rho1) <<  std::endl;
   std::cout << "B electron ehartr = " << myelectron->ehartree(dng1) <<  std::endl;
   std::cout << "electron vl_ave  = " << myelectron->vl_ave(dng1) <<  std::endl;
   std::cout << "electron vnl ave = " << myelectron->vnl_ave(psi1) <<  std::endl;

   //std::cout << "Z electron ehartr = " << myelectron->ehartree(dng1) <<  std::endl;

   //myelectron->run(psi1,rho1,dng1,rho1_all);
   //sum1 = myelectron->energy(psi1,rho1,dng1,rho1_all);
   //std::cout << "Z electron energy = " << sum1 <<  std::endl;


/*---------------------- testing Electron Operators ---------------------- */


}

