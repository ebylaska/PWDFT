
#include "Control2.hpp"
#include "cElectron.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"
#include "CStrfac.hpp"
#include "cpsi.hpp"

#include "Solid.hpp"

namespace pwdft {

/********************************************
 *                                          *
 *              Solid::Solid                *
 *                                          *
 ********************************************/
Solid::Solid(char *infilename, bool wvfnc_initialize, Cneb *mygrid0,
                   Ion *myion0, CStrfac *mystrfac0, Ewald *myewald0,
                   cElectron_Operators *myelectron0, CPseudopotential *mypsp0,
                   std::ostream &coutput) {
  mygrid = mygrid0;
  myion = myion0;
  mystrfac = mystrfac0;
  myewald = myewald0;
  myelectron = myelectron0;
  mypsp = mypsp0;

  nbrillouin = mygrid->nbrillouin;
  nbrillq = mygrid->nbrillq;
  ispin = mygrid->ispin;
  neall = mygrid->neq[0] + mygrid->neq[1];
  ne[0] = mygrid->ne[0];
  ne[1] = mygrid->ne[1];
  nfft[0] = mygrid->nx;
  nfft[1] = mygrid->ny;
  nfft[2] = mygrid->nz;
  for (int i = 0; i < 60; ++i)
    E[i] = 0.0;

  //psi1 = mygrid->g_allocate(1);
  //psi2 = mygrid->g_allocate(1);
  psi1 = mygrid->g_allocate_nbrillq_all();
  psi2 = mygrid->g_allocate_nbrillq_all();
  rho1 = mygrid->r_nalloc(ispin);
  rho2 = mygrid->r_nalloc(ispin);
  rho1_all = mygrid->r_nalloc(ispin);
  rho2_all = mygrid->r_nalloc(ispin);
  dng1 = mygrid->c_pack_allocate(0);
  dng2 = mygrid->c_pack_allocate(0);

  //lmbda = mygrid->m_allocate(-1, 1);
  //hml = mygrid->m_allocate(-1, 1);
  hml   = mygrid->w_allocate_nbrillq_all();
  lmbda = mygrid->w_allocate_nbrillq_all();
  eig   = new double[nbrillq*(ne[0]+ne[1])];

  omega = mygrid->lattice->omega();
  scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
  scal2 = 1.0 / omega;
  dv = omega * scal1;

  n2ft3d = (mygrid->n2ft3d);
  nfft3d = (mygrid->nfft3d);
  //shift1 = 2 * (mygrid->npack(1));
  shift1 = 2 * (mygrid->npack1_max());
  shift2 = (mygrid->n2ft3d);

  newpsi = cpsi_read(mygrid, infilename, wvfnc_initialize, psi1, coutput);

  myelectron->gen_vl_potential();

  /*---------------------- testing Electron Operators ---------------------- */
    /*  
     double sum1;
     std::cout << "Solid: Testing Electron Operators " << std::endl;

     myelectron->run(psi1,rho1,dng1,rho1_all);
     sum1 = myelectron->energy(psi1,rho1,dng1,rho1_all);
     myelectron->gen_hml(psi1,hml);
     std::cout << "electron energy = " << sum1 <<  std::endl;
     std::cout << "electron eorbit = " << myelectron->eorbit(psi1) << std::endl;
     std::cout << "electron ehartr = " << myelectron->ehartree(dng1) <<
     std::endl; std::cout << "electron exc    = " << myelectron->exc(rho1_all)
     <<  std::endl; std::cout << "electron pxc    = " << myelectron->pxc(rho1)
     <<  std::endl; 
     std::cout << "B electron ehartr = " << myelectron->ehartree(dng1) <<  std::endl; 
     std::cout << "electron vl_ave  = " << myelectron->vl_ave(dng1) <<  std::endl; 
     std::cout << "electron vnl ave ="  << myelectron->vnl_ave(psi1) <<  std::endl;
    */
  /*---------------------- testing Electron Operators ---------------------- */
}

} // namespace pwdft
