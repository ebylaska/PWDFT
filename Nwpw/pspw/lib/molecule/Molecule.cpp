
#include "Control2.hpp"
#include "Electron.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "psi.hpp"

#include "Molecule.hpp"

namespace pwdft {

/********************************************
 *                                          *
 *           Molecule::Molecule             *
 *                                          *
 ********************************************/
Molecule::Molecule(char *infilename, bool wvfnc_initialize, Pneb *mygrid0,
                   Ion *myion0, Strfac *mystrfac0, Ewald *myewald0,
                   Electron_Operators *myelectron0, Pseudopotential *mypsp0,  
                   Control2 &control, std::ostream &coutput) 
{
  mygrid = mygrid0;
  myion = myion0;
  mystrfac = mystrfac0;
  myewald = myewald0;
  myelectron = myelectron0;
  mypsp = mypsp0;


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

  psi1 = mygrid->g_allocate(1);
  psi2 = mygrid->g_allocate(1);
  rho1 = mygrid->r_nalloc(ispin);
  rho2 = mygrid->r_nalloc(ispin);
  rho1_all = mygrid->r_nalloc(ispin);
  rho2_all = mygrid->r_nalloc(ispin);
  dng1 = mygrid->c_pack_allocate(0);
  dng2 = mygrid->c_pack_allocate(0);

  lmbda = mygrid->m_allocate(-1, 1);
  hml = mygrid->m_allocate(-1, 1);
  eig = new double[mygrid->ne[0] + mygrid->ne[1]];

  omega = mygrid->lattice->omega();
  scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
  scal2 = 1.0 / omega;
  dv = omega * scal1;

  n2ft3d = (mygrid->n2ft3d);
  shift1 = 2 * (mygrid->npack(1));
  shift2 = (mygrid->n2ft3d);

  newpsi = psi_read(mygrid, infilename, wvfnc_initialize, psi1, coutput);

  myelectron->gen_vl_potential();

  /*---------------------- testing Electron Operators ---------------------- */
  /*
     double sum1;
     std::cout << "Molecule: Testing Electron Operators " << std::endl;

     myelectron->run(psi1,rho1,dng1,rho1_all);
     sum1 = myelectron->energy(psi1,rho1,dng1,rho1_all);
     myelectron->gen_hml(psi1,hml);
     std::cout << "electron energy = " << sum1 <<  std::endl;
     std::cout << "electron eorbit = " << myelectron->eorbit(psi1) << std::endl;
     std::cout << "electron ehartr = " << myelectron->ehartree(dng1) <<
     std::endl; std::cout << "electron exc    = " << myelectron->exc(rho1_all)
     <<  std::endl; std::cout << "electron pxc    = " << myelectron->pxc(rho1)
     <<  std::endl; std::cout << "B electron ehartr = " <<
     myelectron->ehartree(dng1) <<  std::endl; std::cout << "electron vl_ave  =
     " << myelectron->vl_ave(dng1) <<  std::endl; std::cout << "electron vnl ave
     = " << myelectron->vnl_ave(psi1) <<  std::endl;

  */
  /*---------------------- testing Electron Operators ---------------------- */
}


/********************************************
 *                                          *
 *           Molecule::epsi_initialize      *
 *                                          *
 ********************************************/
void Molecule::epsi_initialize(char *infilename, bool wvfnc_initialize, const int *nex, std::ostream &coutput) 
{
   ne_excited[0] = nex[0];
   ne_excited[1] = nex[1];

   psi1_excited = mygrid->g_allocate_excited(nex,1);
   psi2_excited = mygrid->g_allocate_excited(nex,1);

   hml_excited = mygrid->m_nex_allocate(-1, nex);
   eig_excited = new double[nex[0] + nex[1]];
 



   /* read psi from file if psi_exist and not forcing wavefunction initialization */
  bool  newpsi = epsi_read(mygrid, infilename, wvfnc_initialize, nex, psi1_excited, coutput);

  coutput << " input epsi filename:" << infilename << std::endl;
  // to determing ne_excited look at wavefunction or look at control
  //mygrid->g_set_ne_excited(ne_excited);
 // psi1_excited = mygrid->g_allocate_excited(1);
 // psi2_excited = mygrid->g_allocate_excited(1);
  //newpsi = psi_read(mygrid, infilename, wvfnc_initialize, psi1, coutput);
}

/********************************************
 *                                          *
 *           Molecule::epsi_finalize        *
 *                                          *
 ********************************************/
void Molecule::epsi_finalize(char *outfilename, std::ostream &coutput) 
{
   epsi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne_excited,
             psi1_excited,outfilename,coutput);
}

/********************************************
 *                                          *
 *           Molecule::epsi_minimize        *
 *                                          *
 ********************************************/
/**
 * @brief Minimize the energy of excited state orbitals.
 *
 * This function performs a minimization process on the excited state orbitals of a molecule.
 * It iteratively orthogonalizes, normalizes, and minimizes each orbital, ensuring that the
 * final orbitals are properly conditioned and their energies are minimized within a specified
 * tolerance. The process involves projecting out contributions from filled and lower virtual
 * spaces and adjusting the orbital based on convergence criteria.
 *
 * @param vall Pointer to an array of potential values used during the minimization process.
 * @param coutput Output stream for logging or reporting the progress of the minimization.
 *
 * The minimization process includes the following steps:
 * - Orthogonalization against filled and lower virtual spaces.
 * - Normalization of the orbital.
 * - Iterative minimization using the `epsi_KS_update_virtual` function, checking for convergence.
 * - If the minimization does not converge within the allowed number of attempts, the orbital
 *   is reset, and the process is retried.
 * - The final minimized energy of each orbital is stored in the `eig_excited` array.
 */
void Molecule::epsi_minimize(double *vall, std::ostream &coutput)
{
   mygrid->g_ortho_excited(psi1, ne_excited, psi1_excited);
   double error_out,eorb0;

   for (auto ms=0; ms<ispin; ++ms)
   {
      int kshift = ms*ne_excited[0]*2*mygrid->PGrid::npack(1);
      for (auto k=0; k<ne_excited[ms]; ++k)
      {
         int indxk = 2*mygrid->PGrid::npack(1)*k + kshift;
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
                   mygrid->c_pack_zero(1, orb);
                   mygrid->c_pack_addzero(1, 1.0, orb);
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

}

/********************************************
 *                                          *
 *     Molecule::epsi_sort_virtual          *
 *                                          *
 ********************************************/
/**
 * @brief Sorts the excited state orbitals and corresponding eigenvalues in ascending order.
 *
 * This function sorts the excited state orbitals (`psi1_excited`) and their corresponding
 * eigenvalues (`eig_excited`) in ascending order of eigenvalues. It iterates over each pair
 * of orbitals and swaps them if they are out of order based on their eigenvalues.
 *
 * The sorting is performed separately for each spin state. If a swap is needed, both the
 * eigenvalues and the associated orbitals are swapped to maintain consistency.
 *
 * The function allocates a temporary buffer for swapping orbitals and deallocates it at
 * the end of the sorting process.
 */
void Molecule::epsi_sort_virtual()
{
   double *torb = mygrid->c_pack_allocate(1);
   for (auto ms=0; ms<ispin; ++ms)
   {
      int msshift = ms*ne_excited[0]*2*mygrid->PGrid::npack(1);
      for (auto ii=0;    ii<ne_excited[ms]; ++ii)
      for (auto jj=ii+1; jj<ne_excited[ms]; ++jj)
      {
         int indxii = 2*mygrid->PGrid::npack(1)*ii + msshift;
         int indxjj = 2*mygrid->PGrid::npack(1)*jj + msshift;
         double *orbii = psi1_excited + indxii;
         double *orbjj = psi1_excited + indxjj;

         int i = ii + ms*ne_excited[0];
         int j = jj + ms*ne_excited[0];
         double ei = eig_excited[i];
         double ej = eig_excited[j];

         if (ej<ei)
         {
            std::swap(eig_excited[i], eig_excited[j]);
            mygrid->cc_pack_copy(1,orbii,torb);
            mygrid->cc_pack_copy(1,orbjj,orbii);
            mygrid->cc_pack_copy(1,torb,orbjj);
         }
      }
   }
   mygrid->c_pack_deallocate(torb);
}

/********************************************
 *                                          *
 *     Molecule::epsi_KS_update_virtual     *
 *                                          *
 ********************************************/
//    This routine performs a KS update on orbital i
double Molecule::epsi_KS_update_virtual(const int ms, const int k, const int maxit_orb, const double maxerror, 
                                        const double perror, double *vall, double *orb, 
                                         double *error_out, std::ostream &coutput)
{

   double *t0 = mygrid->c_pack_allocate(1);
   double *r1 = mygrid->c_pack_allocate(1);
   double *g  = mygrid->c_pack_allocate(1);
   double *t  = mygrid->c_pack_allocate(1);
  
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
      epsi_get_gradient(orb, vall+ms*n2ft3d, g);
      e0 = mygrid->cc_pack_dot(1,orb,g);

      e0 = -e0;


      double percent_error=0.0;
      if(error0>1.0e-11) percent_error = std::abs(e0-eold)/error0;
     

      precondition = (std::abs(e0-eold)>(sp*maxerror));

      done = ((it > maxit_orb) || (std::abs(e0-eold)<maxerror));
     
      mygrid->cc_pack_copy(1,g,r1);
      mygrid->cc_pack_daxpy(1,(e0),orb,r1);

      //preconditioning 
      if (precondition)
      { 
         ++pit;
         myelectron->get_myke()->ke_precondition(ep,1,orb,r1);
      }

      //determine conjuagate direction ***
      double lmbda_r1 = mygrid->cc_pack_dot(1,r1,r1);

      mygrid->cc_pack_copy(1,r1,t);


      if (it>1) mygrid->cc_pack_daxpy(1,(lmbda_r1/lmbda_r0),t0,t);
      lmbda_r0 = lmbda_r1;
      bool oneloop = true;
      bool repeat_loop = true;
      while (repeat_loop)
      {
        mygrid->cc_pack_copy(1,t,t0);

         //normalize search direction, t ****
         //   call psi_project_out_virtual(ii,dcpl_mb(t(1)))
         // project out filled space
         mygrid->g_project_out_filled(psi1, ms, t);
       
         // project out lower virtual space
         mygrid->g_project_out_virtual(ms, ne_excited, k, psi1_excited, t);

         de0 = mygrid->cc_pack_dot(1,t,t);
         de0 = 1.0/std::sqrt(de0);
         mygrid->c_pack_SMul(1,de0,t);
         de0 = mygrid->cc_pack_dot(1,t,g);

         //bad direction;
         if ((de0<0.0) && oneloop)
         {
             mygrid->cc_pack_copy(1,g,t);
             oneloop = false;
         }
         else
            repeat_loop = false;
      }

      de0 = -2.0*de0;

      epsi_linesearch_update(e0,de0,&theta,vall+ms*n2ft3d,orb,t);

      done = ((it > maxit_orb) ||  (std::abs(e0-eold) < maxerror));
      //done = true;
   }
   mygrid->c_pack_deallocate(t0);
   mygrid->c_pack_deallocate(r1);
   mygrid->c_pack_deallocate(g);
   mygrid->c_pack_deallocate(t);

   *error_out = std::abs(e0-eold);
   e0         = -e0;

   coutput << std::setw(12) << "orbital" << std::setw(4) << k 
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
/**
 * @brief Updates orbital parameters using a line search method based on energy gradients.
 *
 * This function modifies the orbital parameters by applying a trigonometric update 
 * and gradient information to minimize the energy. The update involves calculating 
 * the new angle (`theta`) and rotating the orbitals (`orb`) accordingly.
 *
 * @param e0 Initial energy value.
 * @param de0 Initial energy derivative.
 * @param theta Pointer to the angle parameter.
 * @param orb Pointer to the orbital parameters.
 * @param t Pointer to the search direction.
 *
 * @return void
 */
void Molecule::epsi_linesearch_update(double e0, double de0, double *theta, double *vall, double *orb, double *t)
{
   double *torb = mygrid->c_pack_allocate(1);
   double *g    = mygrid->c_pack_allocate(1);
   double theta0 = *theta;

   mygrid->cc_pack_copy(1,orb, torb);

   // orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
   double x = std::cos(theta0);
   double y = std::sin(theta0);
   mygrid->cc_pack_SMul(1,x,torb,orb);
   mygrid->cc_pack_daxpy(1,y,t,orb);

   // determine theta ***
   epsi_get_gradient(orb, vall, g);
   double e1 = mygrid->cc_pack_dot(1,orb,g);
   e1 = -e1;

   x = (e0 - e1 + 0.5*de0*std::sin(2.0*theta0))/(1.0-std::cos(2*theta0));
   double theta1 = 0.5*std::atan(0.50*de0/x);

   // orb2 = orb*cos(theta) + t*sin(theta) ****
   x = std::cos(theta1);
   y = std::sin(theta1);
   mygrid->cc_pack_SMul(1,x,torb,orb);
   mygrid->cc_pack_daxpy(1,y,t,orb);

   mygrid->c_pack_deallocate(torb);
   mygrid->c_pack_deallocate(g);

   *theta = theta1;
}


/********************************************
 *                                          *
 *     Molecule::epsi_get_gradient          *
 *                                          *
 ********************************************/
void Molecule::epsi_get_gradient(double *orb, double *vall, double *Horb)
{
   double *orb_r = mygrid->r_alloc();

   // fourier transform orb_r
   mygrid->r_zero(orb_r);
   mygrid->cc_pack_copy(1,orb,orb_r);
   mygrid->c_unpack(1,orb_r);
   mygrid->cr_fft3d(orb_r);

   mygrid->c_pack_zero(1,Horb);
   psi_H_orb(mygrid,myelectron->get_myke(),mypsp,orb,orb_r,vall,Horb);
   mygrid->c_pack_SMul(1,-1.0,Horb);

   mygrid->r_dealloc(orb_r);
}



} // namespace pwdft

