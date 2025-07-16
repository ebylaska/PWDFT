
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

   fractional = control.fractional();
   if (fractional)
   {
      nextra[0] = control.fractional_orbitals(0);
      if (control.ispin()==2)
         nextra[1] = control.fractional_orbitals(1);
      else
         nextra[1] = 0;
      smearcorrection = 0.0;
      smeartype = control.fractional_smeartype();
      smearkT   = control.fractional_kT();

      occ1 = mygrid->initialize_occupations_with_allocation(nextra);
      occ2 = mygrid->initialize_occupations_with_allocation(nextra);

      fractional_frozen = control.fractional_frozen();
      fractional_alpha = control.fractional_alpha();
      fractional_alpha_min = control.fractional_alpha_min();
      fractional_alpha_max = control.fractional_alpha_max();
      fractional_beta = control.fractional_beta();
      fractional_gamma = control.fractional_gamma();
      fractional_rmsd_threshold = control.fractional_rmsd_threshold();
      //if (fractional_alpha < fractional_alpha_min) fractional_alpha_min =  fractional_alpha;
      //if (fractional_alpha > fractional_alpha_max) fractional_alpha_max =  fractional_alpha;
      occupation_update = fractional && !control.fractional_frozen();
   }
   else
   {
      nextra[0] = 0;
      nextra[1] = 0;
      fractional_frozen = false;
      fractional_alpha = 0.0;
      fractional_alpha_min = 0.0;
      fractional_alpha_max = 0.0;
      fractional_beta  = 0.0;
      fractional_gamma = 0.0;
      occupation_update = false;
   }
   multiplicity = control.multiplicity();
   total_charge = control.total_charge();


 
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
   hml  = mygrid->w_allocate_nbrillq_all();
   eig  = new double[nbrillq*(ne[0]+ne[1])];
   eig_prev = new double[nbrillq*(ne[0]+ne[1])];
   //hml2  = mygrid->w_allocate_nbrillq_all();
   //eig2  = new double[nbrillq*(ne[0]+ne[1])];
   lmbda = mygrid->w_allocate_nbrillq_all();
 
   omega = mygrid->lattice->omega();
   scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
   scal2 = 1.0 / omega;
   dv = omega * scal1;
 
   n2ft3d = (mygrid->CGrid::n2ft3d);
   nfft3d = (mygrid->CGrid::nfft3d);
   //shift1 = 2 * (mygrid->npack(1));
   shift1 = 2 * (mygrid->CGrid::npack1_max());
   shift2 = (mygrid->CGrid::n2ft3d);
   mshift = 2*(ne[0]*ne[0]+ne[1]*ne[1]);
 
   // Instead of always reading from file, check the force_reinit_flag
   if (force_reinit_flag) {
      if (mygrid->c3db::parall->base_stdio_print)
         coutput << "[PWDFT] Forcing reinitialization of wavefunction due to previous NaN/Inf failure." << std::endl;
      
      // Enhanced fallback strategy with multiple options
      std::string guess = control.initial_wavefunction_guess();
      
      if (guess == "atomic") {
         if (mygrid->c3db::parall->base_stdio_print)
            coutput << "[PWDFT] Using atomic guess for wavefunction reinitialization." << std::endl;
         mygrid->g_generate_atomic_guess(psi1);
      } else if (guess == "superposition") {
         if (mygrid->c3db::parall->base_stdio_print)
            coutput << "[PWDFT] Using superposition of atomic orbitals for wavefunction reinitialization." << std::endl;
         // Generate superposition of atomic orbitals
         mygrid->g_generate_superposition_guess(psi1);
      } else if (guess == "gaussian") {
         if (mygrid->c3db::parall->base_stdio_print)
            coutput << "[PWDFT] Using Gaussian basis guess for wavefunction reinitialization." << std::endl;
         // Generate Gaussian basis functions
         mygrid->g_generate_gaussian_guess(psi1);
      } else if (guess == "mixed") {
         if (mygrid->c3db::parall->base_stdio_print)
            coutput << "[PWDFT] Using mixed strategy (atomic + random) for wavefunction reinitialization." << std::endl;
         // Generate mixed guess: atomic for occupied, random for virtual
         mygrid->g_generate_mixed_guess(psi1);
      } else {
         if (mygrid->c3db::parall->base_stdio_print)
            coutput << "[PWDFT] Using random initialization for wavefunction reinitialization." << std::endl;
         mygrid->g_generate_random(psi1);
      }
      clear_force_reinit_wavefunction();
   } else {
      newpsi = cpsi_read(mygrid, infilename, wvfnc_initialize, psi1, &smearoccupation, occ2, coutput);
   }
   smearoccupation = 0;
   if (fractional)
   {
      if (newpsi) replace_excited_psi1(control,coutput);
      smearoccupation = 1;
      std::vector<double> filling = control.fractional_filling();
      if (filling.size() > 0)
      {
         int sz = filling.size();
         if (sz > (ne[0]+ne[1])) sz = ne[0]+ne[1];
         for (auto nbq=0; nbq<nbrillq; ++nbq)
         {
            int ishift = nbq*(ne[0]+ne[1]);
            std::memcpy(occ2+ishift,filling.data(),sz*sizeof(double));
         }
      }
  
      std::memcpy(occ1,occ2,nbrillq*(ne[0]+ne[1])*sizeof(double));
   }


 
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
 *         Solid::replace_excited_psi1      *
 *                                          *
 ********************************************/
/**
 * @brief Inserts excited-state orbitals from file into the current psi1 array.
 *
 * This routine replaces the top-most (virtual) portion of the `psi1` array
 * with excited-state orbitals read from an external file (e.g., `.emovecs`),
 * using the layout encoded in `Control2::fractional_orbitals()`.
 *
 * Steps performed:
 * 1. Check for existence of input excited-state file via `cpsi_filefind`.
 * 2. Initialize the `psi1_excited` data using `ecpsi_initialize(...)`.
 * 3. Copy excited orbitals into the top of `psi1` for each k-point and spin.
 * 4. Then renames the original `.emovecs` file to `.bak` via `ecpsi_rename`.
 *
 * Assumptions:
 * - The number of excited orbitals per spin is known and passed as `nex[2]`.
 * - Destination `psi1` and source `psi1_excited` are both allocated and valid.
 * - The memory layout uses packed real and imaginary components in interleaved
 *   format (2× `npack`) for each wavefunction.
 * - Indexing accounts for spin (ms), Brillouin zone index (nb), and orbital number (n).
 *
 * @param[in]  control   Reference to the control object containing I/O parameters.
 * @param[out] coutput   Stream used for logging messages (e.g., std::cout).
 */
void Solid::replace_excited_psi1(Control2 &control, std::ostream &coutput)
{
   int nex[2] = { control.fractional_orbitals(0), control.fractional_orbitals(1)};

   if (cpsi_filefind(mygrid,control.input_e_movecs_filename()))
   {
      ecpsi_initialize(control.input_e_movecs_filename(),false,nex,coutput);

      for (auto nb=0; nb<nbrillq; ++nb)
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (nex[ms]>0)
         {
            for (auto n=0; n<nex[ms]; ++n)
            {
               int sz = 2*mygrid->npack(1);
               int indxf = sz*(n + ms*mygrid->neq[0])  + nb*(mygrid->neq[0]+mygrid->neq[1]);
               int indxe = sz*(n + ms*nex[0]) + nb*(nex[0]+nex[1]);;
               std::memcpy(psi1 + indxf, psi1_excited + indxe, sz);
            }
         }
          //epsi_finalize(control.input_e_movecs_filename(),coutput);

      }

      // rename the ecpsi_filename to *.bak 
      ecpsi_rename(control.input_e_movecs_filename(),coutput);
   }
}




/********************************************
 *                                          *
 *         Solid::cpsi_KS_update             *
 *                                          *
 ********************************************/
/**
 * @brief Updates Kohn-Sham orbitals using iterative minimization over all Brillouin points and spins.
 *
 * For each Brillouin zone point and spin channel, this function:
 *   - Projects the current orbital orthogonal to higher occupied states
 *   - Normalizes the orbital
 *   - Applies an iterative update via `cpsi_KS_update_orb`
 *   - Accumulates the energy contribution weighted by Brillouin zone weight
 *
 * This is part of the orbital-by-orbital Kohn-Sham optimization loop.
 *
 * @param[in]  maxit_orb   Maximum iterations per orbital
 * @param[in]  maxerror    Maximum orbital gradient norm allowed
 * @param[in]  perror      Residual tolerance
 * @param[in]  vall        Potential array (must be complex-allocated)
 * @param[in]  ispin       Number of spin channels (1 or 2)
 * @param[in]  neq         Array of orbital counts per spin [n_alpha, n_beta]
 * @param[in]  nbrillq     Number of Brillouin zone points
 * @param[in]  psi         Input/output orbitals, packed complex format
 * @param[out] error_out   Maximum gradient norm encountered
 * @param[in,out] coutput  Output stream for logging
 *
 * @return The total energy summed over all updated orbitals and k-points.
 *
 * @note `vall` must be allocated using `c_alloc` as it's treated as a complex array.
 * @warning Orbitals (`psi`) must be packed consistently across spin and Brillouin index.
 */
double Solid::cpsi_KS_update(const int maxit_orb, const double maxerror,
                             const double perror, double *vall,
                             const int ispin, const int *neq, const int nbrillq, 
                             double *psi, double *error_out, std::ostream &coutput)
{        
   double esum = 0.0;

   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {
      int nbq1 = nbq+1;
      int kshift = nbq*(neq[0]+neq[1])*2*mygrid->CGrid::npack1_max();
      double weight =  mygrid->pbrill_weight(nbq);
      double *psi1k = psi + kshift;


      for (auto ms=0; ms<ispin; ++ms)
      {    
         int ishift = ms*neq[0]*2*mygrid->CGrid::npack1_max();
         for (auto i=neq[ms]-1; i>=0; --i)
         {
            int indx = 2*mygrid->CGrid::npack1_max()*i + ishift;
            double *orb = psi1k + indx;
 
            // orthogonalize to lower orbitals
            //mygrid->g_project_out_filled_below(psi1, ms, i, orb);

            mygrid->g_project_out_filled_above(nbq1,psi1k, ms, i, orb);
 
            // normalize
            double norm = mygrid->cc_pack_dot(nbq1,orb,orb);
            norm = 1.0/std::sqrt(norm);
            mygrid->c_pack_SMul(nbq1,norm,orb);
         
            double e0 = cpsi_KS_update_orb(nbq1,ms,i,maxit_orb,maxerror,perror,vall,psi1k,orb,
                                          error_out, coutput);
 
            esum += e0*weight;
         }
      }
   }

   // --- FULL ORTHOGONALIZATION CLEANUP ---
   mygrid->g_ortho(psi);

   return esum;
}



/********************************************
 *                                          *
 *         Solid::cpsi_KS_update_orb        *
 *                                          *
 ********************************************/
double Solid::cpsi_KS_update_orb(const int nbq1, const int ms, const int k, const int maxit_orb, const double maxerror,
                                 const double perror, double *vall, double *psi1k, double *orb,
                                 double *error_out, std::ostream &coutput)
{

   double *t0 = mygrid->c_pack1_max_allocate();
   double *r1 = mygrid->c_pack1_max_allocate();
   double *g  = mygrid->c_pack1_max_allocate();
   double *t  = mygrid->c_pack1_max_allocate();


   bool precondition = false;
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
      error0 = std::abs(e0-eold);
      eold = e0;

      //calculate residual (steepest descent) direction for a single band
      compute_Horb_for_cg(nbq1,orb,vall+ms*n2ft3d,g);
      e0 = mygrid->cc_pack_dot(nbq1,orb,g);

      e0 = -e0;
     

      double percent_error = 0.0;
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

      mygrid->cc_pack_copy(1,r1,t);

      if (it>1)
         mygrid->cc_pack_daxpy(nbq1,(lmbda_r1/lmbda_r0),t0,t);
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
         //mygrid->g_project_out_filled_below(psi1, ms, k, t);
         mygrid->g_project_out_filled_from_k_up(nbq1,psi1k, ms, k, t);

         de0 = mygrid->cc_pack_dot(nbq1,t,t);
         de0 = 1.0/std::sqrt(de0);
         //if (std::isnan(de0)) de0=0.0;
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

      cpsi_linesearch_update(nbq1, e0,de0,&theta,vall+ms*nfft3d,orb,t);



      done = ((it > maxit_orb) ||  (std::abs(e0-eold) < maxerror));
      //done = true;
   }

   mygrid->c_pack_deallocate(t0);
   mygrid->c_pack_deallocate(r1);
   mygrid->c_pack_deallocate(g);
   mygrid->c_pack_deallocate(t);

   e0         = -e0;
   *error_out = (e0-eold);

   bool lprint = (mygrid->c3db::parall->is_master());
   lprint = false;
   if (lprint) coutput << std::setw(12) << "orbital" << std::setw(4) << k+1
           << " current e=" << std::setw(10) << std::scientific << std::setprecision(3) << e0
           << " (error=" << std::setw(9) << std::scientific << std::setprecision(3) << (*error_out) << ")"
           << " iterations" << std::setw(4) << it << "(" << std::setw(4) << pit
           << " preconditioned, Ep,Sp=" << std::fixed << std::setprecision(1) << std::setw(5) << ep
           << "," << std::setw(7) << sp << ")" << std::endl;


   return e0;
}



void Solid::cpsi_linesearch_update(const int nbq1, double e0, double de0, double *theta, double *vall, double *orb, double *t)
{
   //double *torb = mygrid->c_pack_allocate(nbq1);
   //double *g    = mygrid->c_pack_allocate(nbq1);
   double *torb = mygrid->c_pack1_max_allocate();
   double *g    = mygrid->c_pack1_max_allocate();

   double theta0 = *theta;

   mygrid->cc_pack_copy(nbq1,orb, torb);

   // orb2 = orb*cos(pi/300) + t*sin(pi/300) ****
   double x = std::cos(theta0);
   double y = std::sin(theta0);
   mygrid->cc_pack_SMul(nbq1,x,torb,orb);
   mygrid->cc_pack_daxpy(nbq1,y,t,orb);

   // determine theta ***
   compute_Horb_for_cg(nbq1,orb,vall,g);
   double e1 = mygrid->cc_pack_dot(nbq1,orb,g);
   e1 = -e1;


   x = (e0 - e1 + 0.5*de0*std::sin(2.0*theta0))/(1.0-std::cos(2*theta0));
   //x = (e1 - e0 + 0.5*de0*std::sin(2.0*theta0))/(1.0-std::cos(2*theta0));
   double theta1 = 0.5*std::atan(0.50*de0/x);
   if (std::isnan(theta1)) theta1 =0.0;


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
 *           Solid::ecpsi_initialize        *
 *                                          *
 ********************************************/
void Solid::ecpsi_initialize(char *infilename, bool wvfnc_initialize, const int *nex, std::ostream &coutput) 
{
   ne_excited[0] = nex[0];
   ne_excited[1] = nex[1];

   psi1_excited = mygrid->g_allocate_excited_nbrillq_all(nex);
   psi2_excited = mygrid->g_allocate_excited_nbrillq_all(nex);

   hml_excited = mygrid->w_nex_allocate_nbrillq_all(nex);
   eig_excited = new double[nbrillq*(nex[0] + nex[1])];


   /* read psi from file if psi_exist and not forcing wavefunction initialization */
  bool  newpsi = ecpsi_read(mygrid, infilename, wvfnc_initialize, nex, psi1_excited, coutput);

  bool lprint = (mygrid->c3db::parall->is_master());
  if (lprint) coutput << " input ecpsi filename:" << infilename << std::endl;
  // to determing ne_excited look at wavefunction or look at control
  //mygrid->g_set_ne_excited(ne_excited);
  // psi1_excited = mygrid->g_allocate_excited(1);
  // psi2_excited = mygrid->g_allocate_excited(1);
  //newpsi = psi_read(mygrid, infilename, wvfnc_initialize, psi1, coutput);
}

/********************************************
 *                                          *
 *           Solid::ecpsi_rename            *
 *                                          *
 ********************************************/
/**
 * @brief Renames an excited-state wavefunction file to a backup file.
 * 
 * This function appends the `.bak` extension to the provided filename to preserve
 * the original wavefunction file after it has been used. It only executes the rename
 * operation on the master process to avoid concurrent access in parallel environments.
 * 
 * @param infilename The name of the original wavefunction file to rename.
 * @param coutput    Output stream for logging status messages (e.g., std::cout).
 */
void Solid::ecpsi_rename(const char *infilename, std::ostream &coutput)
{  
   std::string original_filename = infilename;
   std::string backup_filename   = original_filename + ".bak";
   
   // Rename the file to .bak (only master process should do this)
   if (mygrid->c3db::parall->is_master())
   {
      if (std::rename(original_filename.c_str(), backup_filename.c_str()) != 0) 
      {
         coutput << " Warning: Unable to rename " << original_filename
                 << " to " << backup_filename << std::endl;
      } 
      else 
      {
         coutput << " Moved " << original_filename << " to " << backup_filename << std::endl;
      }
   }
}



/********************************************
 *                                          *
 *           Solid::ecpsi_finalize          *
 *                                          *
 ********************************************/
void Solid::ecpsi_finalize(char *outfilename, std::ostream &coutput)
{ 
   ecpsi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne_excited,&nbrillouin,
             psi1_excited,outfilename,coutput);
} 

/********************************************
 *                                          *
 *           Solid::ecpsi_minimize          *
 *                                          *
 ********************************************/
void Solid::ecpsi_minimize(double *vall, std::ostream &coutput)
{ 

   int nshift0 = 2*(mygrid->neq[0]+mygrid->neq[1])*mygrid->CGrid::npack1_max();
   int nshift1 = 2*(ne_excited[0]+ne_excited[1])  *mygrid->CGrid::npack1_max();
   int nshift2 = (ne_excited[0]+ne_excited[1]);
   bool lprint = (mygrid->c3db::parall->is_master());

   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {
      auto nbq1 = nbq+1;

      if (lprint) coutput << std::endl << std::setw(12) << " Brillouin zone point:" << std::setw(4) << nbq+1 
                                                        << " of " << nbrillq <<  std::endl;

      double error_out,eorb0;


      double *psi_f = psi1         + nbq*nshift0;
      double *psi_v = psi1_excited + nbq*nshift1;
      double *eig_v = eig_excited  + nbq*nshift2;

      mygrid->g_ortho_excited(nbq1,psi_f,ne_excited,psi_v);
  
      for (auto ms=0; ms<ispin; ++ms) 
      {
         int kshift = ms*ne_excited[0]*2*mygrid->CGrid::npack1_max();
         for (auto k=0; k<ne_excited[ms]; ++k)
         {
            int indxk = 2*mygrid->CGrid::npack1_max()*k + kshift;
            double *orb = psi_v + indxk;

            //double tum = mygrid->cc_pack_dot(nbq1,orb,orb);
            //std::cout << "NORM k=" << k << " tum=" << tum << std::endl;

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
                  eorb0 = ecpsi_KS_update_virtual(nbq1,ms,k,120,tole,0.001,vall,orb,&error_out, psi_f,psi_v,coutput);
                  //eorb0 = ecpsi_KS_update_virtual(nbq1,ms,k,1,tole,0.001,vall,orb,&error_out, psi_f,psi_v,coutput);
                  if (error_out <= tole)
                     continue_inner_loop = false; // Exit the inner loop
               }

               // If error is still too large or energy is too high, retry orthogonalization
               if (error_out > tole || eorb0 > 4.0)
               {
                   if (l2 <= 1)
                   {
                      //mygrid->c_pack_zero(nbq1, orb);
                      //mygrid->c_pack_addzero(nbq1, 1.0, orb);
                      //int nne[2] = {1,0};
                      //std::cout << "INTO exited_random nne=" << nne[0] << " " << nne[1] <<  std::endl;
                      //mygrid->g_generate_excited_random(nne,orb);
                      //mygrid->g_project_out_filled(nbq1,psi_f,ms,orb);
                      //mygrid->g_project_out_virtual(nbq1,ms,ne_excited,k,psi_v,orb);
                      //mygrid->g_norm(nbq1,orb);
                      mygrid->c_corrector_orb(nbq1,orb);
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

      ecpsi_sort_virtual(nbq1,eig_v,psi_v);
   }
}


/********************************************
 *                                          *
 *           Solid::ecpsi_sort_virtual       *
 *                                          *
 ********************************************/
void Solid::ecpsi_sort_virtual(const int nbq1, double *eig_v, double *psi_v)
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
 *      Solid::ecpsi_KS_update_virtual       *
 *                                          *
 ********************************************/
double Solid::ecpsi_KS_update_virtual(const int nbq1, const int ms, const int k,
                                     const int maxit_orb, const double maxerror,
                                     const double perror, double *vall, double *orb,
                                     double *error_out, 
                                     double *psi_f, double *psi_v, std::ostream &coutput)

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
      compute_Horb_for_cg(nbq1,orb,vall+ms*n2ft3d,g);
      e0 = mygrid->cc_pack_dot(nbq1,orb,g);

      e0 = -e0;

      double percent_error=0.0;
      if(error0>1.0e-11) percent_error = std::abs(e0-eold)/error0;


      precondition = (std::abs(e0-eold)>(sp*maxerror));

      done = ((it > maxit_orb) || (std::abs(e0-eold)<maxerror));

      if (!done)
      {
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
            mygrid->g_project_out_filled(nbq1, psi_f, ms, t);
 
            mygrid->g_project_out_virtual(nbq1, ms, ne_excited, k+1, psi_v, t);
 
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
 
         ecpsi_linesearch_update(nbq1,e0,de0,&theta,vall+ms*nfft3d,orb,t);
 
         done = ((it > maxit_orb) ||  (std::abs(e0-eold) < maxerror));
         //done = true;
      }
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
 *        ecpsi_linesearch_update           *
 *                                          *
 ********************************************/
void Solid::ecpsi_linesearch_update(const int nbq1, 
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
   compute_Horb_for_cg(nbq1,orb, vall, g);
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

   mygrid->c_pack_deallocate(torb);
   mygrid->c_pack_deallocate(g);

   *theta = theta1;
}

/********************************************
 *                                          *
 *        Solid::compute_Horb_for_cg        *
 *                                          *
 ********************************************/
/**
 * @brief Applies the Hamiltonian to a complex orbital and computes -H * orb.
 *
 * This function transforms the input orbital to reciprocal space,
 * applies the Hamiltonian operator, and stores the result in `Horb`,
 * scaled by -1.0. It is used in conjugate gradient (CG) steps for
 * electronic structure optimization.
 *
 * @param[in]  nbq1   Index of the current Brillouin zone point
 * @param[in]  orb    Input complex orbital (packed format, complex-allocated)
 * @param[in]  vall   Local potential (complex grid, complex-allocated)
 * @param[out] Horb   Output buffer for -H * orb (complex-allocated)
 *
 * @note All grid pointers must be allocated using `mygrid->c_alloc()`.
 */
void Solid::compute_Horb_for_cg(const int nbq1, double *orb, double *vall, double *Horb)
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
