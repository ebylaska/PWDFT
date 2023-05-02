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
      dielec    = control.gpoisson_dielec();
      rho0      = control.gpoisson_rho0();
      beta      = control.gpoisson_beta();
      rhomin = control.gpoisson_rhomin();
      rhomax = control.gpoisson_rhomax();

      maxit_pol = control.gpoisson_maxit();
      model_pol = control.gpoisson_model();
      alpha_pol = control.gpoisson_alpha();
      tole_pol  = control.tolerances(0);
      rcut_ion  = control.gpoisson_rcut_ion();

      epsilon = mypneb->r_alloc();
      depsilon = mypneb->r_alloc();
      sw = mypneb->r_alloc();
      p = mypneb->r_alloc();

      epsilon_x = mypneb->r_alloc();
      epsilon_y = mypneb->r_alloc();
      epsilon_z = mypneb->r_alloc();
      w_x = mypneb->r_alloc();
      w_y = mypneb->r_alloc();
      w_z = mypneb->r_alloc();

      rho_ind0 = mypneb->r_alloc();
      rho_ind1 = mypneb->r_alloc();

      rho_ion  = mypneb->r_alloc();
      dng_ion  = mypneb->r_alloc();
      v_ion    = mypneb->r_alloc();
      vdielec0 = mypneb->r_alloc();
      vks0     = mypneb->r_alloc();

      rho_ion_set  = false;
      vdielec0_set  = false;
   }
}


/****************************************************
 *                                                  *
 *    Coulomb12_Operator::initialize_dielectric     *
 *                                                  *
 ****************************************************/
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
   }


   mypneb->rr_Divide(epsilon, depsilon);
   mypneb->r_zero_ends(depsilon);

   mypneb->tcr_pack_iMul_unpack_fft(0,Gx,dng,epsilon_x);
   mypneb->tcr_pack_iMul_unpack_fft(0,Gy,dng,epsilon_y);
   mypneb->tcr_pack_iMul_unpack_fft(0,Gz,dng,epsilon_z);
   mypneb->rr_Mul(depsilon,epsilon_x);
   mypneb->rr_Mul(depsilon,epsilon_y);
   mypneb->rr_Mul(depsilon,epsilon_z);

   for (auto i=0; i<n2ft3d; ++i) 
   {
      rho_ind0[i] = (1.0 / epsilon[i] - 1.0) * (rho[i] + rho_ion[i]);
      p[i] = vh[i] + v_ion[i];
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

   mypneb->rrrr_Sum(vdielec0,vh,v_ion,p);
   mypneb->r_SMul(scal1,p);
   mypneb->rc_pfft3f(0,p);
   mypneb->c_pack(0,p);

   mypneb->tcr_pack_iMul_unpack_fft(0,Gx,p,w_x);
   mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_y);
   mypneb->tcr_pack_iMul_unpack_fft(0,Gz,p,w_z);
   for (auto i=0; i<n2ft3d; ++i)
      vks0[i] = 0.5*overfourpi*depsilon[i]*(w_x[i]*w_x[i] + w_y[i]*w_y[i] + w_z[i]*w_z[i]);

   mypneb->rrr_Sum(vdielec0,vks0,vdielec);
 
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
 *       Coulomb12_Operator::dielectric_fion         *
 *                                                   *
 *****************************************************/
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

 E = 0.5*omega*<dng_ion(G)|vdielec(G)>

 iG*dng_ion(ii,G)*vdielec(G)

 g*(a.x-ia*y)*(b.x+i*b.y)

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
std::string Coulomb12_Operator::shortprint_dielectric()
{

   if (has_dielec){

      std::stringstream stream;

      stream << std::endl;
      stream << " dielectric field:" << std::endl;
      stream << "      dielectric constant -eps- =" << Ffmt(10,3) << dielec << std::endl;
      if (model_pol==0) stream << "      model    =  Andreussi" << std::endl;
      if (model_pol==1) stream << "      model    = Andreussi2" << std::endl;
      if (model_pol==2) stream << "      model    =  Fattebert" << std::endl;
      stream << "      maxit    = " << Ifmt(10)   << maxit_pol << std::endl;
      stream << "      alpha    = " << Ffmt(10,3) << alpha_pol << std::endl;
      stream << "      rcut_ion = " << Ffmt(10,3) << rcut_ion << " au" << std::endl;
      stream << "      rhomin   = " << Efmt(10,3) << rhomin   << " au" << std::endl;
      stream << "      rhomax   = " << Efmt(10,3) << rhomax   << " au" << std::endl;
      stream << "      tole     = " << Efmt(10,3) << tole_pol << std::endl;

      return stream.str();
   }
   else
      return "";
}




} // namespace pwdft
