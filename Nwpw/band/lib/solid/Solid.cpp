
#include "Control2.hpp"
#include "cElectron.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"
#include "CStrfac.hpp"
#include "cpsi.hpp"
#include "cpsi_H.hpp"

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
                   Control2 &control, std::ostream &coutput) 
{
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
   for (int i=0; i<60; ++i)
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

   int nshift0 = 2*(mygrid->neq[0]+mygrid->neq[1])*mygrid->CGrid::npack1_max();
   int nshift1 = 2*(ne_excited[0]+ ne_excited[1])*mygrid->CGrid::npack1_max();
   int nshift2 = (ne_excited[0]+ ne_excited[1]);
   bool lprint = (mygrid->c3db::parall->is_master());

   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {

      if (lprint) coutput << std::endl << std::setw(12) << " Brillouin zone point:" << std::setw(4) << nbq+1 
                                                        << " of " << nbrillq <<  std::endl;

      double error_out,eorb0;

      auto nbq1 = nbq+1;

      double *psi_f = psi1 + nbq*nshift0;
      double *psi_v = psi1_excited + nbq*nshift1;
      double *eig_v = eig_excited + nbq*nshift2;

      mygrid->g_ortho_excited(nbq1,psi_f,ne_excited,psi_v);
  
      for (auto ms=0; ms<ispin; ++ms) 
      {
         int kshift = ms*ne_excited[0]*2*mygrid->CGrid::npack1_max();
         for (auto k=0; k<ne_excited[ms]; ++k)
         {
            int indxk = 2*mygrid->CGrid::npack1_max()*k + kshift;
            double *orb = psi_v + indxk;

            //double tum = mygrid->cc_pack_dot(nbq1,orb,orb);
            //std::cout << "k=" << k << " tum=" << tum << std::endl;

            bool continue_outer_loop = true;
            for (int l2=1; continue_outer_loop && l2<=2; ++l2)
            {
               // Project out filled and lower virtual spaces
               mygrid->g_project_out_filled(nbq1,psi_f, ms, orb);
               mygrid->g_project_out_virtual(nbq1,ms,ne_excited,k,psi_v,orb);

               // normalize 
               mygrid->g_norm(nbq1,orb);


              //for (auto kk=k; kk<ne_excited[ms]; ++kk)
             // {
              //   int indxkk = 2*mygrid->CGrid::npack1_max()*kk + kshift;
               //  double *orbkk = psi_v + indxkk;
                // double sum = mygrid->cc_pack_dot(nbq1,orbkk,orb);
              //   std::cout << "k=" << k << " kk=" << kk << " ms=" << ms << " sum=" << sum << std::endl;

               //  double rum = mygrid->rr_dot(vall,vall);;
                // std::cout << "rum=" << rum << std::endl;
             // }

               // minimize orbital 
               bool continue_inner_loop = true;
               for (int l=0; continue_inner_loop && l<=(1+(l2-1)*3); ++l)
               {
                  eorb0 = epsi_KS_update_virtual(nbq1,ms,k,120,tole,0.001,vall,orb,&error_out,coutput);
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
                      mygrid->g_project_out_filled(nbq1,psi_f,ms,orb);
                      mygrid->g_project_out_virtual(nbq1,ms,ne_excited,k,psi_v,orb);
                      mygrid->g_norm(nbq1,orb);
                   }
                   else
                      continue_outer_loop = false; // Exit the outer loop
               }
               else
                  continue_outer_loop = false; // Exit the outer loop
            }

            // Store the minimized orbital energy
            eig_v[ms*ne_excited[0] + k] = eorb0;
         } //k
      } //ms

      epsi_sort_virtual(nbq1,eig_v,psi_v);
   }
}


/********************************************
 *                                          *
 *           Solid::epsi_sort_virtual       *
 *                                          *
 ********************************************/
void Solid::epsi_sort_virtual(const int nbq1, double *eig_v, double *psi_v)
{
   double *torb = new (std::nothrow) double[2*mygrid->CGrid::npack1_max()]();
   for (auto ms=0; ms<ispin; ++ms)
   {
      int msshift = ms*ne_excited[0]*2*mygrid->CGrid::npack1_max();
      for (auto ii=0;    ii<ne_excited[ms]; ++ii)
      for (auto jj=ii+1; jj<ne_excited[ms]; ++jj)
      {
         int indxii = 2*mygrid->CGrid::npack1_max()*ii + msshift;
         int indxjj = 2*mygrid->CGrid::npack1_max()*jj + msshift;
         double *orbii = psi_v + indxii;
         double *orbjj = psi_v + indxjj;

         int i = ii + ms*ne_excited[0];
         int j = jj + ms*ne_excited[0];
         double ei = eig_v[i];
         double ej = eig_v[j];

         if (ej<ei)
         {
            std::swap(eig_v[i], eig_v[j]);
            mygrid->cc_pack_copy(nbq1,orbii,torb);
            mygrid->cc_pack_copy(nbq1,orbjj,orbii);
            mygrid->cc_pack_copy(nbq1,torb,orbjj);
         }
      }
   }
   mygrid->c_pack_deallocate(torb);

}

/********************************************
 *                                          *
 *      Solid::epsi_KS_update_virtual       *
 *                                          *
 ********************************************/
double Solid::epsi_KS_update_virtual(const int nbq1, const int ms, const int k,
                                     const int maxit_orb, const double maxerror,
                                     const double perror, double *vall, double *orb,
                                     double *error_out, std::ostream &coutput)

{           
   double *t0 = new (std::nothrow) double[2*mygrid->CGrid::npack1_max()]();
   double *r1 = new (std::nothrow) double[2*mygrid->CGrid::npack1_max()]();
   double *g  = new (std::nothrow) double[2*mygrid->CGrid::npack1_max()]();
   double *t  = new (std::nothrow) double[2*mygrid->CGrid::npack1_max()]();

   bool precondition = true;
   bool done = false;
   double error0 = 0.0;
   double e0 = 0.0;
   double eold = 0.0;
   double de0 = 0.0;
   double theta = -3.14159/600.0;
   double    lmbda_r0 = 1.0;
   int it = 0;
   int pit = 0;
   while (!done)
   {
      ++it;
      eold = e0;

      //calculate residual (steepest descent) direction for a single band
      epsi_get_gradient(nbq1,orb, vall+ms*n2ft3d, g);
      e0 = mygrid->cc_pack_dot(nbq1,orb,g);

      e0 = -e0;

      double percent_error=0.0;
      if(error0>1.0e-11) percent_error = std::abs(e0-eold)/error0;


      precondition = (std::abs(e0-eold)>(sp*maxerror));

      done = ((it > maxit_orb) || (std::abs(e0-eold)<maxerror));

      mygrid->cc_pack_copy(nbq1,g,r1);
      mygrid->cc_pack_daxpy(nbq1,(e0),orb,r1);

      //preconditioning 
      if (precondition)
      {
         ++pit;
         myelectron->get_myke()->ke_precondition(nbq1,ep,1,orb,r1);
      }

      //determine conjuagate direction ***
      double lmbda_r1 = mygrid->cc_pack_dot(nbq1,r1,r1);

      mygrid->cc_pack_copy(nbq1,r1,t);


      if (it>1) mygrid->cc_pack_daxpy(nbq1,(lmbda_r1/lmbda_r0),t0,t);
      lmbda_r0 = lmbda_r1;
      bool oneloop = true;
      bool repeat_loop = true;
      while (repeat_loop)
      {
        mygrid->cc_pack_copy(nbq1,t,t0);

         //normalize search direction, t ****
         // project out lower virtual space
         //   call psi_project_out_virtual(ii,dcpl_mb(t(1)))
         // project out filled space
         mygrid->g_project_out_filled(nbq1, psi1, ms, t);

         mygrid->g_project_out_virtual(nbq1, ms, ne_excited, k+1, psi1_excited, t);

         de0 = mygrid->cc_pack_dot(nbq1,t,t);
         de0 = 1.0/std::sqrt(de0);
         mygrid->c_pack_SMul(nbq1,de0,t);
         de0 = mygrid->cc_pack_dot(nbq1,t,g);

         //bad direction;
         if ((de0<0.0) && oneloop)
         {
             mygrid->cc_pack_copy(nbq1,g,t);
             oneloop = false;
         }
         else
            repeat_loop = false;
      }
      de0 = -2.0*de0;

      epsi_linesearch_update(nbq1,e0,de0,&theta,vall+ms*n2ft3d,orb,t);

      done = ((it > maxit_orb) ||  (std::abs(e0-eold) < maxerror));
      //done = true;
   }

   mygrid->c_pack_deallocate(t);
   mygrid->c_pack_deallocate(g);
   mygrid->c_pack_deallocate(r1);
   mygrid->c_pack_deallocate(t0);

   *error_out = std::abs(e0-eold);
   e0         = -e0;

   bool lprint = (mygrid->c3db::parall->is_master());
   if (lprint) coutput << std::setw(18) << "virtual orbital" << std::setw(4) << k+1
           << " current e=" << std::setw(10) << std::scientific << std::setprecision(3) << e0
           << " (error=" << std::setw(9) << std::scientific << std::setprecision(3) << (*error_out) << ")"
           << " iterations" << std::setw(4) << it << "(" << std::setw(4) << pit
           << " preconditioned, Ep,Sp=" << std::fixed << std::setprecision(1) << std::setw(5) << ep
           << "," << std::setw(7) << sp << ")" << std::endl;


   return e0;
}

/********************************************
 *                                          *
 *        epsi_linesearch_update            *
 *                                          *
 ********************************************/
void Solid::epsi_linesearch_update(const int nbq1, 
                                   double e0, double de0, double *theta, double *vall, double *orb, double *t)
{
   double *torb = new (std::nothrow) double[2*mygrid->CGrid::npack1_max()]();
   double *g    = new (std::nothrow) double[2*mygrid->CGrid::npack1_max()]();
   double theta0 = *theta;

   mygrid->cc_pack_copy(nbq1,orb, torb);

   // orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
   double x = std::cos(theta0);
   double y = std::sin(theta0);
   mygrid->cc_pack_SMul(nbq1,x,torb,orb);
   mygrid->cc_pack_daxpy(nbq1,y,t,orb);


   // determine theta ***
   epsi_get_gradient(nbq1,orb, vall, g);
   double e1 = mygrid->cc_pack_dot(nbq1,orb,g);
   e1 = -e1;

   x = (e0 - e1 + 0.5*de0*std::sin(2.0*theta0))/(1.0-std::cos(2*theta0));
   //x = (e1 - e0 + 0.5*de0*std::sin(2.0*theta0))/(1.0-std::cos(2*theta0));
   double theta1 = 0.5*std::atan(0.50*de0/x);

   // orb2 = orb*cos(theta) + t*sin(theta) ****
   x = std::cos(theta1);
   y = std::sin(theta1);

   double sum = mygrid->cc_pack_dot(nbq1,torb,t);
   mygrid->cc_pack_SMul(nbq1,x,torb,orb);
   mygrid->cc_pack_daxpy(nbq1,y,t,orb);
   //std::cout << "theta,x,y=" << std::setprecision(6) << theta1 << " " << x << " " << y <<  " " << e0 << " " << e1 << " " << sum << std::endl;

   mygrid->c_pack_deallocate(torb);
   mygrid->c_pack_deallocate(g);

   *theta = theta1;
}

/********************************************
 *                                          *
 *        Solid::epsi_get_gradient          *
 *                                          *
 ********************************************/
void Solid::epsi_get_gradient(const int nbq1, double *orb, double *vall, double *Horb)
{
   double *orb_r = mygrid->c_alloc();

   // fourier transform orb_r
   mygrid->c_zero(orb_r);
   mygrid->cc_pack_copy(nbq1,orb,orb_r);
   mygrid->c_unpack(nbq1,orb_r);
   mygrid->cr_fft3d(orb_r);

   mygrid->c_pack_zero(nbq1,Horb);
   cpsi_H_orb(nbq1,mygrid,myelectron->get_myke(),mypsp,orb,orb_r,vall,Horb);
   mygrid->c_pack_SMul(nbq1,-1.0,Horb);

   mygrid->c_dealloc(orb_r);
}


      


} // namespace pwdft
