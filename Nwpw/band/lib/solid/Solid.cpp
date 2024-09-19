
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
                   Control2 &control, std::ostream &coutput) {
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

  ep = control.Eprecondition();
  sp = control.Sprecondition();
  tole = control.tolerances(0);

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

/********************************************
 *                                          *
 *           Solid::epsi_initialize         *
 *                                          *
 ********************************************/
void Solid::epsi_initialize(char *infilename, bool wvfnc_initialize, const int *nex, std::ostream &coutput) 
{
   ne_excited[0] = nex[0];
   ne_excited[1] = nex[1];

   psi1_excited = mygrid->g_allocate_excited_nbrillq_all(nex);
   psi2_excited = mygrid->g_allocate_excited_nbrillq_all(nex);

   hml_excited = mygrid->w_nex_allocate_nbrillq_all(nex);
   eig_excited = new double[nbrillq*(nex[0] + nex[1])];


   /* read psi from file if psi_exist and not forcing wavefunction initialization */
  bool  newpsi = epsi_read(mygrid, infilename, wvfnc_initialize, nex, psi1_excited, coutput);

  bool lprint = (mygrid->c3db::parall->is_master());
  if (lprint) coutput << " input epsi filename:" << infilename << std::endl;
  // to determing ne_excited look at wavefunction or look at control
  //mygrid->g_set_ne_excited(ne_excited);
  // psi1_excited = mygrid->g_allocate_excited(1);
  // psi2_excited = mygrid->g_allocate_excited(1);
  //newpsi = psi_read(mygrid, infilename, wvfnc_initialize, psi1, coutput);
}

/********************************************
 *                                          *
 *           Solid::epsi_finalize           *
 *                                          *
 ********************************************/
void Solid::epsi_finalize(char *outfilename, std::ostream &coutput)
{ 
   epsi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne_excited,&nbrillouin,
             psi1_excited,outfilename,coutput);
} 

/********************************************
 *                                          *
 *           Solid::epsi_minimize           *
 *                                          *
 ********************************************/
void Solid::epsi_minimize(double *vall, std::ostream &coutput)
{ 
   int nshift0= 2*(mygrid->neq[0]+mygrid->neq[1])*mygrid->CGrid::npack1_max();
   int nshift1= 2*(ne_excited[0]+ ne_excited[1])*mygrid->CGrid::npack1_max();

   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {
      auto nbq1 = nbq+1;
      mygrid->g_ortho_excited(nbq1,psi1+nbq*nshift0, ne_excited, psi1_excited+nbq*nshift1);
      double error_out,eorb0;
   }
/*
   mygrid->g_ortho_excited(psi1, ne_excited, psi1_excited);
   double error_out,eorb0;
  
   for (auto ms=0; ms<ispin; ++ms) 
   {
      int kshift = ms*ne_excited[0]*2*mygrid->CGrid::npack1_max();
      for (auto k=0; k<ne_excited[ms]; ++k)
      {
         int indxk = 2*mygrid->CGrid::npack1_max()*k + kshift;
         double *orb = psi1_excited + indxk;

         bool continue_outer_loop = true;
         for (int l2=1; continue_outer_loop && l2<=2; ++l2)
         {
            // Project out filled and lower virtual spaces
            mygrid->g_project_out_filled(psi1, ms, orb);
            mygrid->g_project_out_virtual(ms, ne_excited, k, psi1_excited, orb);

            // normalize 
            mygrid->g_norm(orb);

            // minimize orbital 
            bool continue_inner_loop = true;
            for (int l=0; continue_inner_loop && l<=(1+(l2-1)*3); ++l)
            {
               eorb0 = epsi_KS_update_virtual(ms,k,120,tole,0.001,vall,orb,&error_out,coutput);
               if (error_out <= tole)
                  continue_inner_loop = false; // Exit the inner loop
            }

            // If error is still too large or energy is too high, retry orthogonalization
            if (error_out > tole || eorb0 > 4.0)
            {
                if (l2 <= 1)
                {
                   //std::cout << "retry orthogonalization" << std::endl;
                   mygrid->c_pack_zero(1, orb);
                   mygrid->c_pack_addzero(1, 1.0, orb);
                   int nne[2] = {1,0};
                   //std::cout << "INTO exited_random nne=" << nne[0] << " " << nne[1] <<  std::endl;
                   //mygrid->g_generate_excited_random(nne,orb);
                   mygrid->g_project_out_filled(psi1, ms, orb);
                   mygrid->g_project_out_virtual(ms, ne_excited, k, psi1_excited, orb);
                   mygrid->g_norm(orb);
                }
                else
                   continue_outer_loop = false; // Exit the outer loop
            }
            else
               continue_outer_loop = false; // Exit the outer loop
         }

         // Store the minimized orbital energy
         eig_excited[ms*ne_excited[0] + k] = eorb0;
      } //k
   } //ms

  epsi_sort_virtual();
*/
}


} // namespace pwdft
