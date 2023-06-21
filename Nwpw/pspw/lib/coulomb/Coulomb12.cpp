/* Coulomb12.cpp -
   Author - Eric Bylaska
*/

#include "Coulomb12.hpp"
#include "util.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *  Coulomb12_Operator::Coulomb12_Operator *
 *                                         *
 *******************************************/
/*  This constructor initializes the Coulomb12 operator object based on the provided 
 *  Pneb object and Control2 object.                
 *                                               
 *  Parameters:                                 
 *  - mygrid: pointer to Pneb object           
 *  - control: reference to Control2 object   
 */                                           
Coulomb12_Operator::Coulomb12_Operator(Pneb *mygrid, Control2 &control)
{
   mypneb = mygrid;

   if (control.version == 3) 
   {
      has_coulomb1 = true;
      mycoulomb1 = new (std::nothrow) Coulomb_Operator(mygrid, control);
   }

   if (control.version == 4)
   {
      has_coulomb2 = true;
      mycoulomb2 = new (std::nothrow) Coulomb2_Operator(mygrid, control);
   }

   if (control.gpoisson_on())
   {
      has_dielec = true;
      relax_dielec  = control.gpoisson_relax_dielec();
      cube_dielec  = control.gpoisson_cube_dielec();

      dielec        = control.gpoisson_dielec();
      filter_dielec = control.gpoisson_filter();


      rho0       = control.gpoisson_rho0();
      beta       = control.gpoisson_beta();
      rhomin = control.gpoisson_rhomin();
      rhomax = control.gpoisson_rhomax();

      maxit_pol = control.gpoisson_maxit();
      model_pol = control.gpoisson_model();
      alpha_pol = control.gpoisson_alpha();
      tole_pol  = control.tolerances(0);
      rcut_ion  = control.gpoisson_rcut_ion();

      rcenter_pol[0] = control.gpoisson_rcenter(0);
      rcenter_pol[1] = control.gpoisson_rcenter(1);
      rcenter_pol[2] = control.gpoisson_rcenter(2);
      rmin_pol = control.gpoisson_rmin();
      rmax_pol = control.gpoisson_rmax();

      epsilon = mypneb->r_alloc();
      depsilon = mypneb->r_alloc();
      sw = mypneb->r_alloc();
      p = mypneb->r_alloc();

      epsilon_screen = mypneb->r_alloc();
      epsilon_x = mypneb->r_alloc();
      epsilon_y = mypneb->r_alloc();
      epsilon_z = mypneb->r_alloc();
      w_x = mypneb->r_alloc();
      w_y = mypneb->r_alloc();
      w_z = mypneb->r_alloc();
      rho_x = mypneb->r_alloc();
      rho_y = mypneb->r_alloc();
      rho_z = mypneb->r_alloc();

      rho_ind0 = mypneb->r_alloc();
      rho_ind1 = mypneb->r_alloc();

      rho_ion  = mypneb->r_alloc();
      dng_ion  = mypneb->r_alloc();
      v_ion    = mypneb->r_alloc();
      vdielec0 = mypneb->r_alloc();
      vks0     = mypneb->r_alloc();

      //initialize dplot capabilties
      if (cube_dielec)
         tmpcontrol = &control;

      rho_ion_set  = false;
      vdielec0_set  = false;
   }
}


/****************************************************
 *                                                  *
 *    Coulomb12_Operator::initialize_dielectric     *
 *                                                  *
 ****************************************************/
 /*
 *  Initializes the dielectric constant for the Coulomb operator. Sets the ion 
 *  charge density,calculates the electrostatic potential and energy due to the 
 *  ion, and updates the ion-dependent parts of the operator.              
 *                                                
 *  Arguments:                                    
 *  - myion0: pointer to Ion object               
 *  - mystrfac0: pointer to Strfac object        
 *                                              
 *  Returns: none                              
 */
void Coulomb12_Operator::initialize_dielectric(Ion *myion0, Strfac *mystrfac0)
{
   myion    = myion0;
   mystrfac = mystrfac0;

   if ((has_dielec)  && (!rho_ion_set))
   {
      int n2ft3d = mypneb->n2ft3d;
      this->generate_dng_ion(dng_ion);
      std::memcpy(rho_ion,dng_ion,n2ft3d*sizeof(double));
      mypneb->c_unpack(0,rho_ion);
      mypneb->cr_pfft3b(0,rho_ion);
      mypneb->r_zero_ends(rho_ion);
      rho_ion_set = true;

      if (has_coulomb1) mycoulomb1->vcoulomb(dng_ion,v_ion);
      if (has_coulomb2) mycoulomb2->vcoulomb(rho_ion,v_ion);

   }

}


/*******************************************
 *                                         *
 *    Coulomb12_Operator::v_dielectric     *
 *                                         *
 *******************************************/
double Coulomb12_Operator::v_dielectric(const double *rho, const double *dng,
                                        const double *vh, const double *vloc,
                                        double *vdielec) {
  int n2ft3d = mypneb->n2ft3d;
  int npack0 = mypneb->npack(0);
  double *Gx = mypneb->Gpackxyz(0, 0);
  double *Gy = mypneb->Gpackxyz(0, 1);
  double *Gz = mypneb->Gpackxyz(0, 2);
  double omega = mypneb->lattice->omega();
  double scal1 = 1.0 / ((double)((mypneb->nx) * (mypneb->ny) * (mypneb->nz)));
  double scal2 = 1.0 / omega;
  double dv = omega * scal1;
  double fourpi = 16.0 * std::atan(1.0);
  double overfourpi = 1.0 / fourpi;
  double energy = 0.0;

  util_andreussi_dielec(n2ft3d, dielec, rhomin, rhomax, rho, epsilon);
  util_dandreussi_dielec(n2ft3d, dielec, rhomin, rhomax, rho, depsilon);
  mypneb->rr_Divide(epsilon, depsilon);
  mypneb->r_zero_ends(depsilon);

  /* calculate gr = grad n */
  mypneb->tcc_pack_iMul(0, Gx, dng, epsilon_x);
  mypneb->tcc_pack_iMul(0, Gy, dng, epsilon_y);
  mypneb->tcc_pack_iMul(0, Gz, dng, epsilon_z);
  mypneb->c_unpack(0, epsilon_x);
  mypneb->c_unpack(0, epsilon_y);
  mypneb->c_unpack(0, epsilon_z);
  mypneb->cr_pfft3b(0,epsilon_x);
  mypneb->cr_pfft3b(0,epsilon_y);
  mypneb->cr_pfft3b(0,epsilon_z);
  mypneb->rr_Mul(depsilon, epsilon_x);
  mypneb->rr_Mul(depsilon, epsilon_y);
  mypneb->rr_Mul(depsilon, epsilon_z);

  for (auto k = 0; k < npack0; ++k) {
    auto gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
    // p[k]  = scal2*(gg)*vloc[k];
    // p[k]  = (gg)*(vh[0] - vloc[k])/fourpi;
    p[2 * k] = (gg) * (vh[2 * k] - scal2 * vloc[2 * k]) / fourpi;
    p[2 * k + 1] = (gg) * (vh[2 * k + 1] - scal2 * vloc[2 * k + 1]) / fourpi;
    sw[2 * k] = (vh[2 * k] - scal2 * vloc[2 * k]) / fourpi;
    sw[2 * k + 1] = (vh[2 * k + 1] - scal2 * vloc[2 * k + 1]) / fourpi;
  }
  // p[0] = 0.0;
  sw[0] = 0.0;
  std::cout << "p[0]=" << p[0] << " p[1]=" << p[1] << std::endl;
  std::cout << "vloc[0]=" << vloc[0] * scal2 << " vloc[1]=" << vloc[1] * scal2
            << std::endl;
  std::cout << "vh[0]=" << vh[0] << " vh[1]=" << vh[1] << std::endl;
  std::cout << "sw[0]=" << sw[0] << " sw[1]=" << sw[1] << std::endl;

  mypneb->tcc_pack_iMul(0, Gx, sw, w_x);
  mypneb->tcc_pack_iMul(0, Gy, sw, w_y);
  mypneb->tcc_pack_iMul(0, Gz, sw, w_z);
  std::cout << "w_x[0]=" << w_x[0] << " w_x[1]=" << w_x[1] << std::endl;
  std::cout << "w_y[0]=" << w_y[0] << " w_y[1]=" << w_y[1] << std::endl;
  std::cout << "w_z[0]=" << w_z[0] << " w_z[1]=" << w_z[1] << std::endl;

  std::cout << "epsilon_x[0]=" << epsilon_x[0]
            << " epsilon_x[1]=" << epsilon_x[1] << std::endl;
  std::cout << "epsilon_y[0]=" << epsilon_y[0]
            << " epsilon_y[1]=" << epsilon_y[1] << std::endl;
  std::cout << "epsilon_z[0]=" << epsilon_z[0]
            << " epsilon_z[1]=" << epsilon_z[1] << std::endl;

  mypneb->c_unpack(0, w_x);
  mypneb->c_unpack(0, w_y);
  mypneb->c_unpack(0, w_z);
  mypneb->c_unpack(0, p);
  mypneb->cr_fft3d(w_x);
  mypneb->cr_fft3d(w_y);
  mypneb->cr_fft3d(w_z);
  mypneb->cr_fft3d(p);
  mypneb->r_zero_ends(p);
  mypneb->r_zero_ends(w_x);
  mypneb->r_zero_ends(w_y);
  mypneb->r_zero_ends(w_z);

  // calculate rho_ind0
  for (auto i = 0; i < n2ft3d; ++i) {
    // rho_ind0[i] = (1.0-1.0/epsilon[i])*(p[i]) ;
    //             + overfourpi*( epsilon_x[i]*w_x[i]
    //                          + epsilon_y[i]*w_y[i]
    //                          + epsilon_z[i]*w_z[i]);
    // rho_ind0[i] = (1.0-1.0/epsilon[i])*(p[i]) + (epsilon_x[i]*w_x[i] +
    // epsilon_y[i]*w_y[i] + epsilon_z[i]*w_z[i]);

    rho_ind0[i] = (1.0 - 1.0 / epsilon[i]) * (p[i]);
    // rho_ind0[i] = (epsilon_x[i]*w_x[i] + epsilon_y[i]*w_y[i] +
    // epsilon_z[i]*w_z[i]);

    // rho_ind0[i] = (1.0-1.0/epsilon[i])*(rho[i]);
  }
  mypneb->r_zero_ends(rho_ind0);
  energy = mypneb->r_dsum(rho_ind0) * dv;
  std::cout << "rho_ind0*dv=" << energy << std::endl;

  energy = mypneb->r_dsum(p) * dv;
  std::cout << "p*dv=" << energy << std::endl;

  /*
      mypneb->r_SMul(scal1,rho_ind0);
      mypneb->rc_fft3d(rho_ind0);
      mypneb->c_pack(0,rho_ind0);
      mycoulomb1->vcoulomb(rho_ind0,vdielec);

      mypneb->c_unpack(0,vdielec);
      mypneb->cr_fft3d(vdielec);
      mypneb->r_zero_ends(vdielec);

      energy = mypneb->r_dsum(p)*dv;
      std::cout << "p*dv=" << energy << std::endl;
      */

  // std::memcpy(vdielec,rho_ind0,n2ft3d*sizeof(double));
  std::memcpy(vdielec, p, n2ft3d * sizeof(double));

  return energy;
}

/*****************************************************
 *                                                   *
 *    Coulomb12_Operator::v_dielectric_aperiodic     *
 *                                                   *
 *****************************************************/
/**
 * Calculates the aperiodic dielectric contribution to the Coulomb operator
 *
 * The aperiodic dielectric contribution is computed using a Poisson solver,
 * and it depends on the charge density rho and the electronic polarizability
 * tensor dng. The Poisson solver also requires the values of the Hartree
 * potential vh and the ion-induced charge density rho_ion, which are
 * generated or updated by this function if `move` is true or if they have
 * not been previously set.
 *
 * The result is stored in the `vdielec` array, which should have the same size
 * as the `rho` array. The optional `fion` array is used to store the ion-induced
 * forces, if provided and `move` is true.
 *
 * @param rho        Pointer to the charge density array of size `n2ft3d`.
 * @param dng        Pointer to the fourier transform of the density of size `npack0`.
 * @param vh         Pointer to the Hartree potential array of size `n2ft3d`.
 * @param vdielec    Pointer to the output array of size `n2ft3d`.
 * @param move       Whether to re-generate the rho_ion, dng_ion, and v_ion
 *                   arrays. If true, these arrays will be generated or updated
 *                   even if they have been previously set. Otherwise, they will
 *                   be reused if available.
 * @param fion       Optional pointer to the array used to store the ion-induced
 *                   forces. If not null and `move` is true, the forces will be
 *                   computed and stored in this array. Otherwise, the array
 *                   will not be modified.
 */
void Coulomb12_Operator::v_dielectric_aperiodic(const double *rho, const double *dng, const double *vh, 
                                                double *vdielec, bool move, double *fion)
{
   int n2ft3d = mypneb->n2ft3d;
   int npack0 = mypneb->npack(0);
   double *Gx = mypneb->Gpackxyz(0,0);
   double *Gy = mypneb->Gpackxyz(0,1);
   double *Gz = mypneb->Gpackxyz(0,2);
   double omega = mypneb->lattice->omega();
   double scal1 = 1.0/((double)((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
   double scal2 = 1.0/omega;
   double dv = omega*scal1;
   double fourpi = 16.0*std::atan(1.0);
   double overfourpi = 1.0/fourpi;
   double energy = 0.0;

   /* re-calcuate rho_ion, dng_ion, v_ion */
   if ((move) || (!rho_ion_set))
   {
      this->generate_dng_ion(dng_ion);
      std::memcpy(rho_ion,dng_ion,n2ft3d*sizeof(double));
      mypneb->c_unpack(0,rho_ion);
      mypneb->cr_pfft3b(0,rho_ion);
      mycoulomb2->vcoulomb(rho_ion,v_ion);
 
      rho_ion_set = true;
   }

   /* generate caclulate dielectric and dielectric gradients */
   this->dielectric_generate(rho,dng);

/*
   mypneb->rr_Divide(epsilon, depsilon);
   mypneb->r_zero_ends(depsilon);

   mypneb->tcr_pack_iMul_unpack_fft(0,Gx,dng,epsilon_x);
   mypneb->tcr_pack_iMul_unpack_fft(0,Gy,dng,epsilon_y);
   mypneb->tcr_pack_iMul_unpack_fft(0,Gz,dng,epsilon_z);
   mypneb->rr_Mul(depsilon,epsilon_x);
   mypneb->rr_Mul(depsilon,epsilon_y);
   mypneb->rr_Mul(depsilon,epsilon_z);
   */

   for (auto i=0; i<n2ft3d; ++i) 
   {
      //rho_ind0[i] = (1.0/epsilon[i]-1.0) * (rho[i] + rho_ion[i]);
      rho_ind0[i] = epsilon_screen[i]*(rho[i]+rho_ion[i]);
      p[i]        = vh[i]+v_ion[i];
   }
   mypneb->r_zero_ends(rho_ind0);
   mypneb->r_zero_ends(p);

   mypneb->r_SMul(scal1,p);
   mypneb->rc_pfft3f(0,p);
   mypneb->c_pack(0,p);

   mypneb->tcr_pack_iMul_unpack_fft(0,Gx,p,w_x);
   mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_y);
   mypneb->tcr_pack_iMul_unpack_fft(0,Gz,p,w_z);
   for (auto i=0; i<n2ft3d; ++i)
      rho_ind0[i] += overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]);


   std::memcpy(rho_ind1,rho_ind0,n2ft3d*sizeof(double));
   if (vdielec0_set)
   {
      mypneb->rr_SMul(scal1,vdielec0,p);
      mypneb->rc_pfft3f(0,p);
      mypneb->c_pack(0,p);

      mypneb->tcr_pack_iMul_unpack_fft(0,Gx,p,w_x);
      mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_y);
      mypneb->tcr_pack_iMul_unpack_fft(0,Gz,p,w_z);
      for (auto i=0; i<n2ft3d; ++i)
         rho_ind1[i] += overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]);
   }
   mycoulomb2->vcoulomb(rho_ind1,vdielec0);
   mypneb->r_zero_ends(vdielec0);

   /* iteration=0,1,... */
   double eold = 0.0;
   double epol = 0.5*mypneb->rr_dot(rho_ind1,vdielec0)*dv;
   int it = 0;
   while ((std::abs(epol-eold) > tole_pol) && (it<maxit_pol)) 
   {
      ++it;
      mypneb->rr_SMul(scal1,vdielec0,p);
      mypneb->rc_pfft3f(0,p);
      mypneb->c_pack(0,p);
      mypneb->tcr_pack_iMul_unpack_fft(0,Gx,p,w_x);
      mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_y);
      mypneb->tcr_pack_iMul_unpack_fft(0,Gz,p,w_z);
      for (auto i=0; i<n2ft3d; ++i) 
      {
       rho_ind1[i] = (1.0 - alpha_pol)*rho_ind1[i]
                   + alpha_pol*(rho_ind0[i] 
                           + overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]));
      }
      mycoulomb2->vcoulomb(rho_ind1,vdielec0);
      mypneb->r_zero_ends(vdielec0);

      eold = epol;
      epol = 0.5*mypneb->rr_dot(rho_ind1,vdielec0)*dv;

/*
      double energy = mypneb->r_dsum(rho_ind1) * dv;
      double sumion = mypneb->r_dsum(rho_ion) * dv;
      double sumelc = mypneb->r_dsum(rho) * dv;
      double eelc = mypneb->rr_dot(rho,vdielec0) * dv;
      double eion = mypneb->rr_dot(rho_ion,vdielec0) * dv;
      if (mypneb->d3db::parall->is_master())
      {
        std::cout << Ifmt(5) << it << Efmt(15,6) << energy << " " << sumelc
                  << " " << sumion << " " << epol << " " << eelc << " " << eion
                  << " " << epol-eold << std::endl;
      }
*/
   }

   if (relax_dielec)
   {
      mypneb->rrrr_Sum(vdielec0,vh,v_ion,p);
      mypneb->r_SMul(scal1,p);
      mypneb->rc_pfft3f(0,p);
      mypneb->c_pack(0,p);
 
      mypneb->tcr_pack_iMul_unpack_fft(0,Gx,p,w_x);
      mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_y);
      mypneb->tcr_pack_iMul_unpack_fft(0,Gz,p,w_z);
      for (auto i=0; i<n2ft3d; ++i)
         vks0[i] = -0.5*overfourpi*depsilon[i]*(w_x[i]*w_x[i] + w_y[i]*w_y[i] + w_z[i]*w_z[i]);
                 //- 0.5*(w_x[i]*rho_x[i] + w_y[i]*rho_y[i] + w_x[i]*rho_z[i]);
 
      mypneb->rrr_Sum(vdielec0,vks0,vdielec);
   }
   else
      std::memcpy(vdielec,vdielec0,n2ft3d*sizeof(double));
 
   double eelc = 0.5*mypneb->rr_dot(rho,vdielec0)*dv;
   double eion = 0.5*mypneb->rr_dot(rho_ion,vdielec0)*dv;
   double pelc = mypneb->rr_dot(rho,vdielec)*dv;


/*
   mypneb->rr_SMul(scal1,vdielec0,p);
   mypneb->r_zero_ends(p);
   mypneb->rc_pfft3f(0,p);
   mypneb->c_pack(0,p);
   double egion = 0.5*omega*mypneb->cc_pack_dot(0,p,dng_ion);
   std::cout << "eion=" << eion << " eion(g)=" << egion<< std::endl;
   std::cout << "omega*dng(0)=" << omega*dng_ion[0] << std::endl;
   */


   edielec = eelc + eion;
   pdielec = pelc;
   vdielec0_set = true;

    /* calculate force */
   if (move)
   {
      mypneb->rr_SMul(scal1,vdielec0,p);
      mypneb->rc_pfft3f(0,p);
      mypneb->c_pack(0,p);
      this->dng_ion_vdielec0_fion(p,fion);

      //std::memset(fion,0,3*12*sizeof(double));
      //for (auto ii=0; ii<12; ++ii)
      //   std::cout << "ii=" << ii << " fion= " << fion[3*ii] << " " << fion[3*ii+1] << " " << fion[3*ii+2] << std::endl;
   }

}

/*****************************************************
 *                                                   *
 *       Coulomb12_Operator::dielectric_generate     *
 *                                                   *
 *****************************************************/
/* 
 *   Description: This method calculates the dielectric constant and        
 *                dielectric gradients required for the calculation of the  
 *                Coulomb 1-2 operator. The dielectric constant is generated 
 *                on a grid of points and may be filtered with a periodic    
 *                Gaussian filter. This function receives the density and    
 *                the fourier transform of the density as input parameters.                     
 *                                                                           
 *   Parameters:  rho - An array containing the density values.              
 *                dng - An array containing the Fourier tranform density values.    
 *                                                                           
 *   Remarks:     This function is called only when the dielectric constant  
 *                is not set or when it needs to be updated. The dielectric  
 *                constant is used by the Coulomb12_Operator to calculate    
 *                electrostatic forces and energies.                         
 */                                                                          
void Coulomb12_Operator::dielectric_generate(const double *rho, const double *dng)
{
   if (notset_dielec || relax_dielec)
   {
      int n2ft3d = mypneb->n2ft3d;
      if (model_pol==3)
      {

         mypneb->initialize_r_grid();
         double *r_grid = mypneb->r_grid;

         util_sphere_dielec(n2ft3d,r_grid,rcenter_pol,dielec,rmin_pol,rmax_pol,epsilon);
         util_sphere_gradient_dielec(n2ft3d,r_grid,rcenter_pol,dielec,rmin_pol,rmax_pol,
                                     epsilon_x,epsilon_y,epsilon_z);
         mypneb->r_zero_ends(epsilon);
         mypneb->r_zero_ends(epsilon_x);
         mypneb->r_zero_ends(epsilon_y);
         mypneb->r_zero_ends(epsilon_z);
      }
      else
      {
         double *Gx = mypneb->Gpackxyz(0,0);
         double *Gy = mypneb->Gpackxyz(0,1);
         double *Gz = mypneb->Gpackxyz(0,2);
       
         /* generate caclcuate dielectric */
         if (model_pol==0)
         {
            util_andreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,epsilon);
            util_dandreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,depsilon);
         }
       
         if (model_pol==1)
         {
            util_andreussi2_dielec(n2ft3d,dielec,rhomin,rhomax,rho,epsilon);
            util_dandreussi2_dielec(n2ft3d,dielec,rhomin,rhomax,rho,depsilon);
         }
         if (model_pol==2)
         {
            util_fattebert_dielec(n2ft3d,dielec,beta,rho0,rho,epsilon);
            util_dfattebert_dielec(n2ft3d,dielec,beta,rho0,rho,depsilon);
         }
       
       
         /* generate caclulate dielectric gradients */
         mypneb->rr_Divide(epsilon, depsilon);
         mypneb->r_zero_ends(epsilon);
         mypneb->r_zero_ends(depsilon);
         if (filter_dielec>0.0)
         {
            mypneb->rr_periodic_gaussian_filter(filter_dielec,epsilon,sw);  mypneb->rr_copy(sw,epsilon);
            mypneb->rr_periodic_gaussian_filter(filter_dielec,depsilon,sw); mypneb->rr_copy(sw,depsilon);
         }
       
         mypneb->tcr_pack_iMul_unpack_fft(0,Gx,dng,rho_x);
         mypneb->tcr_pack_iMul_unpack_fft(0,Gy,dng,rho_y);
         mypneb->tcr_pack_iMul_unpack_fft(0,Gz,dng,rho_z);

         mypneb->tcr_pack_iMul_unpack_fft(0,Gx,dng,epsilon_x);
         mypneb->tcr_pack_iMul_unpack_fft(0,Gy,dng,epsilon_y);
         mypneb->tcr_pack_iMul_unpack_fft(0,Gz,dng,epsilon_z);
         mypneb->rr_Mul(depsilon,epsilon_x);
         mypneb->rr_Mul(depsilon,epsilon_y);
         mypneb->rr_Mul(depsilon,epsilon_z);

      }


      //screen_epsilon[i] = (1.0/epsilon[i]-1.0)
      mypneb->rr_screen0(epsilon,epsilon_screen);

      // Gaussian filter
      /*
      if (filter_dielec>0.0)
      {
         double *tmp = mypneb->r_alloc();
         mypneb->rr_periodic_gaussian_filter(filter_dielec,epsilon,tmp);   mypneb->rr_copy(tmp,epsilon);
         mypneb->rr_periodic_gaussian_filter(filter_dielec,depsilon,tmp);  mypneb->rr_copy(tmp,depsilon);
         mypneb->rr_periodic_gaussian_filter(filter_dielec,epsilon_x,tmp); mypneb->rr_copy(tmp,epsilon_x); 
         mypneb->rr_periodic_gaussian_filter(filter_dielec,epsilon_y,tmp); mypneb->rr_copy(tmp,epsilon_y);
         mypneb->rr_periodic_gaussian_filter(filter_dielec,epsilon_z,tmp); mypneb->rr_copy(tmp,epsilon_z);
         mypneb->r_dealloc(tmp);
      }
      */
      

      /*
      mypneb->rr_periodic_gaussian_filter(1.0,epsilon_x,tmp); mypneb->rr_copy(tmp,epsilon_x);
      mypneb->rr_periodic_gaussian_filter(1.0,epsilon_y,tmp); mypneb->rr_copy(tmp,epsilon_y);
      mypneb->rr_periodic_gaussian_filter(1.0,epsilon_z,tmp); mypneb->rr_copy(tmp,epsilon_z);
      mypneb->r_dealloc(tmp);
      mypneb->r_zero_ends(epsilon);
      mypneb->r_zero_ends(epsilon_x);
      mypneb->r_zero_ends(epsilon_y);
      mypneb->r_zero_ends(epsilon_z);
      mydplot.gcube_write("smooth_epsilon_x.cube", -1, "SCF dielec function",epsilon_x);
      */
     
   }
}

/*****************************************************
 *                                                   *
 *       Coulomb12_Operator::dielectric_fion         *
 *                                                   *
 *****************************************************/
 /*
 * Function: Coulomb12_Operator::dielectric_fion
 * ----------------------------------------------
 * This function calculates the force on ions due to the dielectric field. The input argument is a pointer to an
 * array fion that will contain the calculated force values. The function first calculates the electrostatic potential
 * due to the dielectric field and then calculates the force by taking the negative gradient of the potential.
 * The calculated force values are stored in the array fion.
 *
 * Parameters:
 *    - fion: pointer to an array of double values that will contain the calculated force values
 *
 * Returns: void
 */

void Coulomb12_Operator::dielectric_fion(double *fion)
{
   double scal1 = 1.0/((double)((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));

   mypneb->rr_SMul(scal1,vdielec0,p);
   mypneb->rc_pfft3f(0,p);
   mypneb->c_pack(0,p);
   this->dng_ion_vdielec0_fion(p,fion);
}

/*****************************************************
 *                                                   *
 *    Coulomb12_Operator::rho_ion_vdielec_fion       *
 *                                                   *
 *****************************************************/
 /*
                  /
                  |
 dE/dR =d/dR 0.5* | rho_ion(r)*vdielec(r) dr
                  |
                  /

 * Computes the force on the ions due to the potential from the dielectric
 * using the charge density rho_ion and the vdielec operator. The force is
 * stored in fion. The calculation is based on the equation:
 *
 * dE/dR = d/dR 0.5 * integral(rho_ion(r) * vdielec(r) dr)
 *
 * The energy E is given by:
 *
 * E = 0.5 * omega * <dng_ion(G) | vdielec(G)>
 *
 * where omega is the cell volume, G is a wavevector, and dng_ion(G) is the
 * Fourier transform of rho_ion.
 *
 * Inputs:
 *   - vg: the vdielec operator evaluated on the potential (double array)
 *   - fion: the force on the ions (double array)
 *
 * Outputs:
 *   - fion: the force on the ions (double array)
 */

void Coulomb12_Operator::dng_ion_vdielec0_fion(const double *vg, double *fion)
{
   int nion  = myion->nion;
   int *katm = myion->katm;
   int npack0 = mypneb->npack(0);

   double *Gx = mypneb->Gpackxyz(0,0);
   double *Gy = mypneb->Gpackxyz(0,1);
   double *Gz = mypneb->Gpackxyz(0,2);

   double *zv_psp  = myion->zv_psp;
   double *gauss   = mypneb->t_pack_allocate(0);
   double *xtmp    = mypneb->t_pack_allocate(0);
   //double *xtmp2   = mypneb->t_pack_allocate(0);
   double *exi     = mypneb->c_pack_allocate(0);
   double *dng_tmp = mypneb->c_pack_allocate(0);

   mypneb->t_pack_gaussian(0,rcut_ion,gauss);
   for (auto ii=0; ii<nion; ++ii)
   {
      mystrfac->strfac_pack(0,ii,exi);
      mypneb->tcc_pack_aMul(0,zv_psp[katm[ii]],gauss,exi,dng_tmp);

      //mypneb->cct_pack_iconjgMulb(0,vg,dng_tmp,xtmp);
      mypneb->cct_pack_iconjgMulb(0,dng_tmp,vg,xtmp);
      //mypneb->cct_pack_iconjgMul(0,dng_tmp,vg,xtmp2);
      //for (auto k=0; k<npack0; ++k)
      //    xtmp[k] = 0.5*(xtmp[k] - xtmp2[k]);

      fion[3*ii]   += mypneb->tt_pack_dot(0,Gx,xtmp);
      fion[3*ii+1] += mypneb->tt_pack_dot(0,Gy,xtmp);
      fion[3*ii+2] += mypneb->tt_pack_dot(0,Gz,xtmp);
   }

   mypneb->c_pack_deallocate(dng_tmp);
   mypneb->c_pack_deallocate(exi);
   mypneb->t_pack_deallocate(xtmp);
   //mypneb->t_pack_deallocate(xtmp2);
   mypneb->t_pack_deallocate(gauss);
}

/*****************************************************
 *                                                   *
 *    Coulomb12_Operator::generate_dng_ion           *
 *                                                   *
 *****************************************************/
 /*
 * Function name: Coulomb12_Operator::generate_dng_ion
 *
 * Purpose: Generates the dng_ion array, which contains the ionic density
 *          on the simulation grid. The ionic density is computed as a sum
 *          of the individual Gaussian densities of the ions in the system.
 *
 * Input:
 *   - none
 *
 * Output:
 *   - dng_ion: double pointer to the array containing the ionic density on
 *              the simulation grid
 *
 * Description:
 *   This function calculates the ionic density on the simulation grid by
 *   summing the individual Gaussian densities of the ions in the system. The
 *   Gaussian density for each ion is calculated using the strfac_pack function
 *   of the StructureFactor class and then added to the dng_ion array using the
 *   tcc_pack_aMulAdd function of the Tensor class. The Gaussian density is
 *   evaluated up to the cutoff radius rcut_ion. The ionic density is then scaled
 *   by the reciprocal of the volume of the simulation cell and stored in the
 *   dng_ion array.
 */

void Coulomb12_Operator::generate_dng_ion(double *dng_ion) 
{
   int nion = myion->nion;
   int *katm = myion->katm;
   double *gauss = mypneb->t_pack_allocate(0);
   double *exi   = mypneb->c_pack_allocate(0);
   double *zv_psp = myion->zv_psp;
   double omega = mypneb->lattice->omega();
   double scal2 = 1.0 / omega;
 
   mypneb->c_pack_zero(0,dng_ion);
   mypneb->t_pack_gaussian(0,rcut_ion,gauss);
   for (auto ii=0; ii<nion; ++ii)
   {
      mystrfac->strfac_pack(0,ii,exi);
      mypneb->tcc_pack_aMulAdd(0,zv_psp[katm[ii]],gauss,exi,dng_ion);
   }
   mypneb->c_pack_SMul(0,-scal2,dng_ion);

   mypneb->c_pack_deallocate(exi);
   mypneb->t_pack_deallocate(gauss);
}

/*****************************************************
 *                                                   *
 *    Coulomb12_Operator::v_dielectric2_aperiodic    *
 *                                                   *
 *****************************************************/
/*
let phi = phi_elc + phi_ion + phi_ind
    rho = rho_elc + rho_ion

    \nabla dot (eps\nabla phi) = -4*pi*rho

    eps*\nabla2 phi = -4*pi*rho - \nabla eps dot \nabla phi

    \nabla2 phi = -4*pi*rho/dps - \nabla eps/eps dot \nabla phi

    \nabla2 phi = -4*pi*rho + 4*pi*(1-1/eps)*rho - \nabla eps/eps dot \nabla phi

    \nabla2 phi = -4*pi*rho + 4*pi*(1-1/eps)*rho - \nabla eps/eps dot \nabla
phi_elc
                                                 - \nabla eps/eps dot \nabla
phi_ion
                                                 - \nabla eps/eps dot \nabla
phi_ind

    \nabla2 phi_ind = 4*pi*(1-1/eps)*rho - \nabla eps/eps dot \nabla phi_elc
                                         - \nabla eps/eps dot \nabla phi_ion
                                         - \nabla eps/eps dot \nabla phi_ind

    \nabla2 phi_ind = -4*pi*(1/eps-1)*rho - 4*pi*1/4pi*\nabla eps/eps dot \nabla
phi_elc
                                          - 4*pi*1/4pi*\nabla eps/eps dot \nabla
phi_ion
                                          - 4*pi*1/4pi*\nabla eps/eps dot \nabla
phi_ind

  therefore,

    \nabla2 phi_ind = -4*pi*(rho_ind0 +rho_ind1)

     where rho_ind0 = (1/eps-1)*(rho_ele+rho_ion) + 1/4pi * \nabla eps/eps dot
(\nabla phi_elc)
                                                  + 1/4pi * \nabla eps/eps dot
(\nabla phi_ion)

           rho_ind1 = 1/4pi * \nabla eps/eps dot \nabla phi_ind

*/
double Coulomb12_Operator::v_dielectric2_aperiodic(
    const double *rho, const double *dng, const double *rho_ion,
    const double *dng_ion, const double *vh, const double *vloc,
    const bool restart, double *vdielec, double *vks_dielec,
    nwpw_dplot *mydplot)
{

  bool debug = true;
  int n2ft3d = mypneb->n2ft3d;
  int npack0 = mypneb->npack(0);
  double *Gx = mypneb->Gpackxyz(0, 0);
  double *Gy = mypneb->Gpackxyz(0, 1);
  double *Gz = mypneb->Gpackxyz(0, 2);
  double omega = mypneb->lattice->omega();
  double scal1 = 1.0 / ((double)((mypneb->nx) * (mypneb->ny) * (mypneb->nz)));
  double scal2 = 1.0 / omega;
  double dv = omega * scal1;
  double fourpi = 16.0 * std::atan(1.0);
  double overfourpi = 1.0 / fourpi;
  double energy = 0.0;

  util_andreussi_dielec(n2ft3d, dielec, rhomin, rhomax, rho, epsilon);
  util_dandreussi_dielec(n2ft3d, dielec, rhomin, rhomax, rho, depsilon);
  mypneb->rr_Divide(epsilon, depsilon);
  mypneb->r_zero_ends(depsilon);

  /* calculate gr = grad epsilonn */
  mypneb->tcc_pack_iMul(0, Gx, dng, epsilon_x);
  mypneb->tcc_pack_iMul(0, Gy, dng, epsilon_y);
  mypneb->tcc_pack_iMul(0, Gz, dng, epsilon_z);
  mypneb->c_unpack(0, epsilon_x);
  mypneb->c_unpack(0, epsilon_y);
  mypneb->c_unpack(0, epsilon_z);
  mypneb->cr_fft3d(epsilon_x);
  mypneb->cr_fft3d(epsilon_y);
  mypneb->cr_fft3d(epsilon_z);
  mypneb->cr_fft3d(sw);
  mypneb->rr_Mul(depsilon, epsilon_x);
  mypneb->rr_Mul(depsilon, epsilon_y);
  mypneb->rr_Mul(depsilon, epsilon_z);

  for (auto i = 0; i < n2ft3d; ++i) {
    rho_ind0[i] = (1.0 / epsilon[i] - 1.0) * (rho[i] + rho_ion[i]);
    sw[i] = vh[i] + vloc[i];
  }
  mypneb->r_zero_ends(rho_ind0);
  mypneb->r_zero_ends(sw);

  mypneb->r_SMul(scal1, sw);
  mypneb->rc_fft3d(sw);
  mypneb->c_pack(0, sw);

  mypneb->tcc_pack_iMul(0, Gx, sw, w_x);
  mypneb->tcc_pack_iMul(0, Gy, sw, w_y);
  mypneb->tcc_pack_iMul(0, Gz, sw, w_z);
  mypneb->c_unpack(0, w_x);
  mypneb->c_unpack(0, w_y);
  mypneb->c_unpack(0, w_z);
  mypneb->cr_fft3d(w_x);
  mypneb->cr_fft3d(w_y);
  mypneb->cr_fft3d(w_z);
  mypneb->r_zero_ends(w_x);
  mypneb->r_zero_ends(w_y);
  mypneb->r_zero_ends(w_z);
  for (auto i = 0; i < n2ft3d; ++i)
    rho_ind0[i] += overfourpi * (w_x[i] * epsilon_x[i] + w_y[i] * epsilon_y[i] +
                                 w_z[i] * epsilon_z[i]);

  /* step A */
  std::memcpy(rho_ind1, rho_ind0, n2ft3d * sizeof(double));
  if (restart) {
    mypneb->rr_SMul(scal1, vdielec, p);
    mypneb->rc_fft3d(p);
    mypneb->c_pack(0, p);

    mypneb->tcr_pack_iMul_unpack_fft(0, Gx, p, w_x);
    mypneb->tcr_pack_iMul_unpack_fft(0, Gy, p, w_y);
    mypneb->tcr_pack_iMul_unpack_fft(0, Gz, p, w_z);
    for (auto i = 0; i < n2ft3d; ++i)
      rho_ind1[i] +=
          overfourpi * (w_x[i] * epsilon_x[i] + w_y[i] * epsilon_y[i] +
                        w_z[i] * epsilon_z[i]);
  }
  mycoulomb2->vcoulomb(rho_ind1, vdielec);
  mypneb->r_zero_ends(vdielec);
  mydplot->gcube_write("AP_rho_ind00.cube", -1, "SCF dielec function",
                       rho_ind1);
  mydplot->gcube_write("AP_vdielec00.cube", -1, "SCF dielec function", vdielec);

  if (mypneb->d3db::parall->is_master() && debug)
    std::cout << "   it    rho_ind1*dv   rho_elc*dv    rho_ion*dv         epol "
                 "    eelc-ind      eion-ind     epol-eold"
              << std::endl;

  /* iteration=0,1,... */
  double eold = 0.0;
  double epol = 0.5 * mypneb->rr_dot(rho_ind1, vdielec) * dv;
  double alpha = 0.41;
  int it = 0;
  while (std::abs(epol - eold) > tole_pol) {
    ++it;
    /* step A */
    mypneb->rr_SMul(scal1, vdielec, p);
    mypneb->rc_fft3d(p);
    mypneb->c_pack(0, p);
    mypneb->tcr_pack_iMul_unpack_fft(0, Gx, p, w_x);
    mypneb->tcr_pack_iMul_unpack_fft(0, Gy, p, w_y);
    mypneb->tcr_pack_iMul_unpack_fft(0, Gz, p, w_z);
    /*
           mypneb->tcc_pack_iMul(0,Gx,p,w_x);
           mypneb->tcc_pack_iMul(0,Gy,p,w_y);
           mypneb->tcc_pack_iMul(0,Gz,p,w_z);
           mypneb->c_unpack(0,w_x);
           mypneb->c_unpack(0,w_y);
           mypneb->c_unpack(0,w_z);
           mypneb->cr_fft3d(w_x);
           mypneb->cr_fft3d(w_y);
           mypneb->cr_fft3d(w_z);
           mypneb->r_zero_ends(w_x);
           mypneb->r_zero_ends(w_y);
           mypneb->r_zero_ends(w_z);
           */
    for (auto i=0; i<n2ft3d; ++i)
    {
       rho_ind1[i] = (1.0 - alpha_pol)*rho_ind1[i]
                   + alpha_pol*(rho_ind0[i] 
                           + overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]));
    }
    mycoulomb2->vcoulomb(rho_ind1, vdielec);
    mypneb->r_zero_ends(vdielec);

    /* step B */

    eold = epol;
    epol = 0.5 * mypneb->rr_dot(rho_ind1, vdielec) * dv;

    energy = mypneb->r_dsum(rho_ind1) * dv;
    double sumion = mypneb->r_dsum(rho_ion) * dv;
    double sumelc = mypneb->r_dsum(rho) * dv;
    double eelc = mypneb->rr_dot(rho, vdielec) * dv;
    double eion = mypneb->rr_dot(rho_ion, vdielec) * dv;
    if (mypneb->d3db::parall->is_master() && debug) {
      std::cout << Ifmt(5) << it << Efmt(15, 6) << energy << " " << sumelc
                << " " << sumion << " " << epol << " " << eelc << " " << eion
                << " " << epol - eold << std::endl;
    }
  }
  mydplot->gcube_write("AP_vdielecff.cube", -1, "SCF dielec function", vdielec);
  mydplot->gcube_write("AP_rho_indff.cube", -1, "SCF dielec function",
                       rho_ind1);
  mydplot->gcube_write("AP_epsilon.cube", -1, "SCF dielec function", epsilon);
  mydplot->gcube_write("AP_epsilon_x.cube", -1, "SCF dielec function",
                       epsilon_x);
  mydplot->gcube_write("AP_epsilon_y.cube", -1, "SCF dielec function",
                       epsilon_y);
  mydplot->gcube_write("AP_epsilon_z.cube", -1, "SCF dielec function",
                       epsilon_z);

  std::memcpy(sw, rho_ion, n2ft3d * sizeof(double));
  mydplot->gcube_write("AP_rho_ion.cube", -1, "SCF dielec function", sw);
  std::memcpy(sw, rho, n2ft3d * sizeof(double));
  mydplot->gcube_write("AP_rho_elc.cube", -1, "SCF dielec function", sw);

  /* */
  mypneb->rrr_Sum(vdielec, vloc, p);
  mypneb->rr_Sum(vh, p);
  mypneb->r_SMul(scal1, p);
  mypneb->rc_fft3d(p);
  mypneb->c_pack(0, p);
  mypneb->tcr_pack_iMul_unpack_fft(0, Gx, p, w_x);
  mypneb->tcr_pack_iMul_unpack_fft(0, Gy, p, w_y);
  mypneb->tcr_pack_iMul_unpack_fft(0, Gy, p, w_z);
  for (auto i = 0; i < n2ft3d; ++i)
    vks_dielec[i] = -0.5 * overfourpi * depsilon[i] *
                    (w_x[i] * w_x[i] + w_y[i] * w_y[i] + w_z[i] * w_z[i]);

  mydplot->gcube_write("AP_vks_dielec.cube", -1, "SCF dielec function",
                       vks_dielec);

  energy = mypneb->rr_dot(rho, vdielec);
  std::cout << "KS potential    rho*vdielec*dv = " << Efmt(15, 6) << energy
            << std::endl;
  energy = mypneb->rr_dot(rho, vks_dielec);
  std::cout << "KS potential rho*vks_dielec*dv = " << Efmt(15, 6) << energy
            << std::endl;

  return epol;
}

/****************************************************
 *                                                  *
 *    Coulomb12_Operator::shortprint_dielectric     *
 *                                                  *
 ****************************************************/
 /*  The function `std::string Coulomb12_Operator::shortprint_dielectric()` is a 
     member function of the `Coulomb12_Operator` class. It returns a string that 
     contains a short summary of the parameters and settings related to the dielectric calculation. 

     The function first checks if the `has_dielec` flag is set to true. If it is, then the function 
     creates a `stringstream` object named `stream` and writes a header for the dielectric field section.

     It then prints several key parameters related to the dielectric calculation, including whether or 
     not to relax the dielectric field, whether or not to cube the dielectric field, the dielectric 
     constant, and if applicable, the Gaussian filter parameters. The function then prints the polarization 
     model, maximum number of iterations, alpha value, and other related parameters.

     If the `model_pol` is set to 3, then the function also prints additional parameters related to the 
     sphere model, such as the minimum and maximum radii and the center of the sphere. If `model_pol` is 
     set to 2, then the function prints the rho0 and beta parameters. If `model_pol` is not equal to 2 or 3, 
     then the function prints the minimum and maximum rho values.

     Finally, the function returns the `stream` object as a string if `has_dielec` is true, otherwise it 
     returns an empty string.
*/

std::string Coulomb12_Operator::shortprint_dielectric()
{
   if (has_dielec){

      std::stringstream stream;

      stream << std::endl;
      stream << " dielectric field:" << std::endl;
      if (relax_dielec)
         stream << "      relax dielectric          = true"  << std::endl;
      else
         stream << "      relax dielectric          = false" << std::endl;
      if (cube_dielec)
         stream << "      cube dielectric           = true"  << std::endl;
      else
         stream << "      cube dielectric           = false" << std::endl;
      stream << "      dielectric constant -eps- = " << LFfmt(10,3) << dielec << std::endl;
      if (filter_dielec>0.0)
         stream << "      Gaussian filter -sigma-   = " << LFfmt(10,3)   << filter_dielec << std::endl;
      if (model_pol==0) stream << "      model    =  Andreussi" << std::endl;
      if (model_pol==1) stream << "      model    = Andreussi2" << std::endl;
      if (model_pol==2) stream << "      model    =  Fattebert" << std::endl;
      if (model_pol==3) stream << "      model    =  sphere"    << std::endl;
      stream << "      maxit    = " << Ifmt(10)   << maxit_pol << std::endl;
      stream << "      alpha    = " << Ffmt(10,3) << alpha_pol << std::endl;
      stream << "      rcut_ion = " << Ffmt(10,3) << rcut_ion << " au" << std::endl;
      if (model_pol==3) 
      {
         stream << "      rmin     = " << Efmt(10,3) << rmin_pol << " au" << std::endl;
         stream << "      rmax     = " << Efmt(10,3) << rmax_pol << " au" << std::endl;
         stream << "      rcenter  = <" << Ffmt(10,3) << rcenter_pol[0] << " " 
                                        << Ffmt(10,3) << rcenter_pol[1] << " " 
                                        << Ffmt(10,3) << rcenter_pol[2] << "> au" << std::endl;
      }
      else if (model_pol==2) 
      {
         stream << "      rho0     = " << Efmt(10,3) << rho0   << " au" << std::endl;
         stream << "      beta     = " << Efmt(10,3) << beta   << " au" << std::endl;
      }
      else
      {
         stream << "      rhomin   = " << Efmt(10,3) << rhomin   << " au" << std::endl;
         stream << "      rhomax   = " << Efmt(10,3) << rhomax   << " au" << std::endl;
      }
      stream << "      tole     = " << Efmt(10,3) << tole_pol << std::endl;

      return stream.str();
   }
   else
      return "";
}




} // namespace pwdft
