/* Coulomb12.cpp - 
   Author - Eric Bylaska
*/

#include "util.hpp"
#include	"Coulomb12.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *  Coulomb12_Operator::Coulomb12_Operator *
 *                                         *
 *******************************************/
Coulomb12_Operator::Coulomb12_Operator(Pneb *mygrid, Control2& control)
{
   mypneb = mygrid;

   if (control.version==3)
   {
      has_coulomb1 = true;
      mycoulomb1   = new (std::nothrow) Coulomb_Operator(mygrid,control);
   }
   
   if (control.version==4) 
   {
      has_coulomb2 = true;
      mycoulomb2   = new (std::nothrow) Coulomb2_Operator(mygrid,control);
   }


   if (control.gpoisson_on())
   {
      has_dielec = true;
      dielec = control.gpoisson_dielec();
      rho0   = control.gpoisson_rho0();
      beta   = control.gpoisson_beta();
      rhomin = control.gpoisson_rhomin();
      rhomax = control.gpoisson_rhomax();
      tole_pol = control.tolerances(0);

      epsilon   = mypneb->r_alloc();
      depsilon  = mypneb->r_alloc();
      ddepsilon = mypneb->r_alloc();
      sw        = mypneb->r_alloc();
      p         = mypneb->r_alloc();

      epsilon_x = mypneb->r_alloc();
      epsilon_y = mypneb->r_alloc();
      epsilon_z = mypneb->r_alloc();
      w_x = mypneb->r_alloc();
      w_y = mypneb->r_alloc();
      w_z = mypneb->r_alloc();

      epsilon_lap = mypneb->r_alloc();

      rho_ind0 = mypneb->r_alloc();
      rho_ind1 = mypneb->r_alloc();
      rho_ion  = mypneb->r_alloc();
   }
}

/*******************************************
 *                                         *
 *    Coulomb12_Operator::v_dielectric     *
 *                                         *
 *******************************************/
double Coulomb12_Operator::v_dielectric(const double *rho, const double *dng,  
                                        const double *vh,  const double *vloc, 
                                        double *vdielec) 
{
    int n2ft3d = mypneb->n2ft3d;
    int npack0 = mypneb->npack(0);
    double *Gx = mypneb->Gpackxyz(0,0);
    double *Gy = mypneb->Gpackxyz(0,1);
    double *Gz = mypneb->Gpackxyz(0,2);
    double omega = mypneb->lattice->omega();
    double scal1  = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
    double scal2  = 1.0/omega;
    double dv = omega*scal1;
    double fourpi = 16.0*std::atan(1.0);
    double overfourpi = 1.0/fourpi;
    double energy = 0.0;

    util_andreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,epsilon);
    util_dandreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,depsilon);
    mypneb->rr_Divide(epsilon,depsilon);
    mypneb->r_zero_ends(depsilon);

    /* calculate gr = grad n */
    mypneb->tcc_pack_iMul(0,Gx,dng,epsilon_x);
    mypneb->tcc_pack_iMul(0,Gy,dng,epsilon_y);
    mypneb->tcc_pack_iMul(0,Gz,dng,epsilon_z);
    mypneb->c_unpack(0,epsilon_x);
    mypneb->c_unpack(0,epsilon_y);
    mypneb->c_unpack(0,epsilon_z);
    mypneb->cr_fft3d(epsilon_x);
    mypneb->cr_fft3d(epsilon_y);
    mypneb->cr_fft3d(epsilon_z);
    mypneb->rr_Mul(depsilon,epsilon_x);
    mypneb->rr_Mul(depsilon,epsilon_y);
    mypneb->rr_Mul(depsilon,epsilon_z);

    for (auto k=0; k<npack0; ++k)
    {
       auto gg = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
       //p[k]  = scal2*(gg)*vloc[k];
       //p[k]  = (gg)*(vh[0] - vloc[k])/fourpi;
       p[2*k]    = (gg)*(vh[2*k]  - scal2*vloc[2*k])  /fourpi;
       p[2*k+1]  = (gg)*(vh[2*k+1]- scal2*vloc[2*k+1])/fourpi;
       sw[2*k]   = (vh[2*k]   - scal2*vloc[2*k])  /fourpi;
       sw[2*k+1] = (vh[2*k+1] - scal2*vloc[2*k+1])/fourpi;
    }
    //p[0] = 0.0;
    sw[0] = 0.0;
    std::cout << "p[0]=" << p[0] << " p[1]=" << p[1] <<  std::endl;
    std::cout << "vloc[0]=" << vloc[0]*scal2 << " vloc[1]=" << vloc[1]*scal2 <<  std::endl;
    std::cout << "vh[0]=" << vh[0] << " vh[1]=" << vh[1] <<  std::endl;
    std::cout << "sw[0]=" << sw[0] << " sw[1]=" << sw[1] <<  std::endl;

    mypneb->tcc_pack_iMul(0,Gx,sw,w_x);
    mypneb->tcc_pack_iMul(0,Gy,sw,w_y);
    mypneb->tcc_pack_iMul(0,Gz,sw,w_z);
    std::cout << "w_x[0]=" << w_x[0] << " w_x[1]=" << w_x[1] <<  std::endl;
    std::cout << "w_y[0]=" << w_y[0] << " w_y[1]=" << w_y[1] <<  std::endl;
    std::cout << "w_z[0]=" << w_z[0] << " w_z[1]=" << w_z[1] <<  std::endl;

    std::cout << "epsilon_x[0]=" << epsilon_x[0] << " epsilon_x[1]=" << epsilon_x[1] <<  std::endl;
    std::cout << "epsilon_y[0]=" << epsilon_y[0] << " epsilon_y[1]=" << epsilon_y[1] <<  std::endl;
    std::cout << "epsilon_z[0]=" << epsilon_z[0] << " epsilon_z[1]=" << epsilon_z[1] <<  std::endl;

    mypneb->c_unpack(0,w_x);
    mypneb->c_unpack(0,w_y);
    mypneb->c_unpack(0,w_z);
    mypneb->c_unpack(0,p);
    mypneb->cr_fft3d(w_x);
    mypneb->cr_fft3d(w_y);
    mypneb->cr_fft3d(w_z);
    mypneb->cr_fft3d(p);
    mypneb->r_zero_ends(p);
    mypneb->r_zero_ends(w_x);
    mypneb->r_zero_ends(w_y);
    mypneb->r_zero_ends(w_z);


    // calculate rho_ind0
    for (auto i=0; i<n2ft3d; ++i)
    {
        //rho_ind0[i] = (1.0-1.0/epsilon[i])*(p[i]) ;
        //            + overfourpi*( epsilon_x[i]*w_x[i]
        //                         + epsilon_y[i]*w_y[i]
        //                         + epsilon_z[i]*w_z[i]);
        //rho_ind0[i] = (1.0-1.0/epsilon[i])*(p[i]) + (epsilon_x[i]*w_x[i] + epsilon_y[i]*w_y[i] + epsilon_z[i]*w_z[i]);
                  
        rho_ind0[i] = (1.0-1.0/epsilon[i])*(p[i]);
        //rho_ind0[i] = (epsilon_x[i]*w_x[i] + epsilon_y[i]*w_y[i] + epsilon_z[i]*w_z[i]);

        //rho_ind0[i] = (1.0-1.0/epsilon[i])*(rho[i]);
    }
    mypneb->r_zero_ends(rho_ind0);
    energy = mypneb->r_dsum(rho_ind0)*dv;
    std::cout << "rho_ind0*dv=" << energy << std::endl;

    energy = mypneb->r_dsum(p)*dv;
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

    //std::memcpy(vdielec,rho_ind0,n2ft3d*sizeof(double));
    std::memcpy(vdielec,p,n2ft3d*sizeof(double));


    return energy;
}

/*****************************************************
 *                                                   *
 *    Coulomb12_Operator::v_dielectric_aperiodic     *
 *                                                   *
 *****************************************************/
double Coulomb12_Operator::v_dielectric_aperiodic(const double *rho, const double *dng,  
                                                  const double *vh,  const double *vloc, 
                                                  double *vdielec,
                                                  nwpw_dplot *mydplot) 
{

    int n2ft3d = mypneb->n2ft3d;
    int npack0 = mypneb->npack(0);
    double *Gx = mypneb->Gpackxyz(0,0);
    double *Gy = mypneb->Gpackxyz(0,1);
    double *Gz = mypneb->Gpackxyz(0,2);
    double omega = mypneb->lattice->omega();
    double scal1  = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
    double scal2  = 1.0/omega;
    double dv = omega*scal1;
    double fourpi = 16.0*std::atan(1.0);
    double overfourpi = 1.0/fourpi;
    double energy = 0.0;

    util_andreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,epsilon);
    util_dandreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,depsilon);
    mypneb->rr_Divide(epsilon,depsilon);
    mypneb->r_zero_ends(depsilon);

    /* calculate gr = grad epsilonn */
    mypneb->tcc_pack_iMul(0,Gx,dng,epsilon_x);
    mypneb->tcc_pack_iMul(0,Gy,dng,epsilon_y);
    mypneb->tcc_pack_iMul(0,Gz,dng,epsilon_z);
    mypneb->c_unpack(0,epsilon_x);
    mypneb->c_unpack(0,epsilon_y);
    mypneb->c_unpack(0,epsilon_z);
    mypneb->cr_fft3d(epsilon_x);
    mypneb->cr_fft3d(epsilon_y);
    mypneb->cr_fft3d(epsilon_z);
    mypneb->rr_Mul(depsilon,epsilon_x);
    mypneb->rr_Mul(depsilon,epsilon_y);
    mypneb->rr_Mul(depsilon,epsilon_z);
    mypneb->r_zero_ends(epsilon_x);
    mypneb->r_zero_ends(epsilon_y);
    mypneb->r_zero_ends(epsilon_z);

    mypneb->rr_SMul(scal1,vloc,rho_ion);
    mypneb->r_zero_ends(rho_ion);
    mypneb->rc_fft3d(rho_ion);
    mypneb->c_pack(0,rho_ion);

    for (auto k=0; k<npack0; ++k)
    {
        auto gg = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
        rho_ion[2*k]   *= (gg/fourpi);
        rho_ion[2*k+1] *= (gg/fourpi);
    }
    //rho_ion[0] = -rho[0];
    //rho_ion[1] = -rho[1];
    mypneb->c_unpack(0,rho_ion);
    mypneb->cr_fft3d(rho_ion);
    mypneb->r_zero_mends(rho_ion);
    mydplot->gcube_write("AP1_rho_ion0.cube",-1,"SCF dielec function",rho_ion);

    for (auto i=0; i<n2ft3d; ++i)
    {
        rho_ind0[i] = (1.0/epsilon[i]-1.0)*(rho[i]+rho_ion[i]);
        sw[i]       = vh[i] + vloc[i];
    }
    mypneb->r_zero_mends(rho_ind0);
    mypneb->r_zero_ends(sw);

    mypneb->r_SMul(scal1,sw);
    mypneb->rc_fft3d(sw);
    mypneb->c_pack(0,sw);

    mypneb->tcc_pack_iMul(0,Gx,sw,w_x);
    mypneb->tcc_pack_iMul(0,Gy,sw,w_y);
    mypneb->tcc_pack_iMul(0,Gz,sw,w_z);
    mypneb->c_unpack(0,w_x);
    mypneb->c_unpack(0,w_y);
    mypneb->c_unpack(0,w_z);
    mypneb->cr_fft3d(w_x);
    mypneb->cr_fft3d(w_y);
    mypneb->cr_fft3d(w_z);
    mypneb->r_zero_mends(w_x);
    mypneb->r_zero_mends(w_y);
    mypneb->r_zero_mends(w_z);
    for (auto i=0; i<n2ft3d; ++i)
    {
        rho_ind0[i] += overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]);
    }
    mypneb->r_zero_mends(rho_ind0);
    mycoulomb2->vcoulomb(rho_ind0,vdielec);
    mypneb->r_zero_ends(vdielec);

    mydplot->gcube_write("AP1_rho_ind0.cube",-1,"SCF dielec function",rho_ind0);


     std::cout << "   it    rho_ind1*dv    rho_ion*dv         epol     eelc-ind      eion-ind     epol-eold" << std::endl;
    /* iteration=0,1,... */
    std::memcpy(rho_ind1,rho_ind0,n2ft3d*sizeof(double));
    double eold = 0.0;
    double epol = 0.5*mypneb->rr_dot(rho_ind1,vdielec)*dv;
    double alpha = 0.41;
    int it = 0;
    while (std::abs(epol-eold)>tole_pol)
    {
       ++it;
       mypneb->rr_SMul(scal1,vdielec,p);
       mypneb->rc_fft3d(p);
       mypneb->c_pack(0,p);

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
       for (auto i=0; i<n2ft3d; ++i)
       {
           rho_ind1[i] = (1.0-alpha)*rho_ind1[i] +
                       + alpha*(rho_ind0[i] + overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]));
       }
       mypneb->r_zero_mends(rho_ind1);
       mycoulomb2->vcoulomb(rho_ind1,vdielec);
       mypneb->r_zero_ends(vdielec);

       eold = epol;
       epol = 0.5*mypneb->rr_dot(rho_ind1,vdielec)*dv;

       energy = mypneb->r_dsum(rho_ind1)*dv;
       double sumion = mypneb->r_dsum(rho_ion)*dv;
       double eelc = mypneb->rr_dot(rho,vdielec)*dv;
       double eion = mypneb->rr_dot(rho_ion,vdielec)*dv;
       std::cout << Ifmt(5) << it << Efmt(15,6) << energy << " "
                                                << sumion << " "
                                                << epol   << " "
                                                << eelc   << " "
                                                << eion   << " "
                                                << epol-eold << std::endl;
    }
    mydplot->gcube_write("AP1_vdielecff.cube",-1,"SCF dielec function",vdielec);
    mydplot->gcube_write("AP1_rho_indff.cube",-1,"SCF dielec function",rho_ind1);
    mydplot->gcube_write("AP1_epsilon.cube",-1,"SCF dielec function",epsilon);
    mydplot->gcube_write("AP1_epsilon_x.cube",-1,"SCF dielec function",epsilon_x);
    mydplot->gcube_write("AP1_epsilon_y.cube",-1,"SCF dielec function",epsilon_y);
    mydplot->gcube_write("AP1_epsilon_z.cube",-1,"SCF dielec function",epsilon_z);
    mydplot->gcube_write("AP1_epsilon_lap.cube",-1,"SCF dielec function",epsilon_lap);

    std::memcpy(sw,rho_ion,n2ft3d*sizeof(double));
    mydplot->gcube_write("AP1_rho_ion.cube",-1,"SCF dielec function",sw);
    std::memcpy(sw,rho,n2ft3d*sizeof(double));
    mydplot->gcube_write("AP1_rho_elc.cube",-1,"SCF dielec function",sw);
    std::memcpy(sw,vloc,n2ft3d*sizeof(double));
    mydplot->gcube_write("AP1_vloc.cube",-1,"SCF dielec function",sw);




    return energy;
}

/*****************************************************
 *                                                   *
 *    Coulomb12_Operator::generate_dng_ion           *
 *                                                   *
 *****************************************************/
void Coulomb12_Operator::generate_dng_ion(Pneb *mypneb, Ion *myion, Strfac *mystrfac, double rc,  double *dng_ion)
{
    int nion  = myion->nion;
    int *katm = myion->katm;
    double *gauss  = mypneb->t_pack_allocate(0);
    double *zv_psp = myion->zv_psp;
    double *exi = mypneb->r_alloc();
    double omega = mypneb->lattice->omega();
    double scal2  = 1.0/omega;

    mypneb->c_pack_zero(0,dng_ion);
    mypneb->t_pack_gaussian(0,rc,gauss);
    for (auto ii=0; ii<nion; ++ii)
    {
       mystrfac->strfac_pack(0,ii,exi);
       mypneb->tcc_pack_aMulAdd(0,zv_psp[katm[ii]],gauss,exi,dng_ion);
    }
    mypneb->c_pack_SMul(0,-scal2,dng_ion);
    mypneb->r_dealloc(exi);
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

     \nabla2 phi = -4*pi*rho + 4*pi*(1-1/eps)*rho - \nabla eps/eps dot \nabla phi_elc
                                                  - \nabla eps/eps dot \nabla phi_ion
                                                  - \nabla eps/eps dot \nabla phi_ind

     \nabla2 phi_ind = 4*pi*(1-1/eps)*rho - \nabla eps/eps dot \nabla phi_elc
                                          - \nabla eps/eps dot \nabla phi_ion
                                          - \nabla eps/eps dot \nabla phi_ind

     \nabla2 phi_ind = -4*pi*(1/eps-1)*rho - 4*pi*1/4pi*\nabla eps/eps dot \nabla phi_elc
                                           - 4*pi*1/4pi*\nabla eps/eps dot \nabla phi_ion
                                           - 4*pi*1/4pi*\nabla eps/eps dot \nabla phi_ind

   therefore,

     \nabla2 phi_ind = -4*pi*(rho_ind0 +rho_ind1) 

      where rho_ind0 = (1/eps-1)*(rho_ele+rho_ion) + 1/4pi * \nabla eps/eps dot (\nabla phi_elc)
                                                   + 1/4pi * \nabla eps/eps dot (\nabla phi_ion)

            rho_ind1 = 1/4pi * \nabla eps/eps dot \nabla phi_ind

*/
double Coulomb12_Operator::v_dielectric2_aperiodic(const double *rho,     const double *dng,  
                                                   const double *rho_ion, const double *dng_ion,  
                                                   const double *vh,      const double *vloc, 
                                                   const bool restart, double *vdielec, double *vks_dielec,
                                                   nwpw_dplot *mydplot) 
{
    bool debug = true;
    int n2ft3d = mypneb->n2ft3d;
    int npack0 = mypneb->npack(0);
    double *Gx = mypneb->Gpackxyz(0,0);
    double *Gy = mypneb->Gpackxyz(0,1);
    double *Gz = mypneb->Gpackxyz(0,2);
    double omega = mypneb->lattice->omega();
    double scal1  = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
    double scal2  = 1.0/omega;
    double dv = omega*scal1;
    double fourpi = 16.0*std::atan(1.0);
    double overfourpi = 1.0/fourpi;
    double energy = 0.0;

    util_andreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,epsilon);
    util_dandreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,depsilon);
    util_ddandreussi_dielec(n2ft3d,dielec,rhomin,rhomax,rho,ddepsilon);
    mypneb->rr_Divide(epsilon,depsilon);
    mypneb->rr_Divide(epsilon,ddepsilon);
    mypneb->r_zero_ends(depsilon);
    mypneb->r_zero_ends(ddepsilon);

    /* calculate gr = grad epsilonn */
    mypneb->tcc_pack_iMul(0,Gx,dng,epsilon_x);
    mypneb->tcc_pack_iMul(0,Gy,dng,epsilon_y);
    mypneb->tcc_pack_iMul(0,Gz,dng,epsilon_z);
    mypneb->c_unpack(0,epsilon_x);
    mypneb->c_unpack(0,epsilon_y);
    mypneb->c_unpack(0,epsilon_z);
    mypneb->cr_fft3d(epsilon_x);
    mypneb->cr_fft3d(epsilon_y);
    mypneb->cr_fft3d(epsilon_z);
    mypneb->cr_fft3d(sw);
    mypneb->rr_Mul(depsilon,epsilon_x);
    mypneb->rr_Mul(depsilon,epsilon_y);
    mypneb->rr_Mul(depsilon,epsilon_z);
    

    for (auto i=0; i<n2ft3d; ++i)
    {
        rho_ind0[i] = (1.0/epsilon[i]-1.0)*(rho[i]+rho_ion[i]);
        sw[i]       = vh[i] + vloc[i];
    }
    mypneb->r_zero_ends(rho_ind0);
    mypneb->r_zero_ends(sw);

    mypneb->r_SMul(scal1,sw);
    mypneb->rc_fft3d(sw);
    mypneb->c_pack(0,sw);

    mypneb->tcc_pack_iMul(0,Gx,sw,w_x);
    mypneb->tcc_pack_iMul(0,Gy,sw,w_y);
    mypneb->tcc_pack_iMul(0,Gz,sw,w_z);
    mypneb->c_unpack(0,w_x);
    mypneb->c_unpack(0,w_y);
    mypneb->c_unpack(0,w_z);
    mypneb->cr_fft3d(w_x);
    mypneb->cr_fft3d(w_y);
    mypneb->cr_fft3d(w_z);
    mypneb->r_zero_ends(w_x);
    mypneb->r_zero_ends(w_y);
    mypneb->r_zero_ends(w_z);
    for (auto i=0; i<n2ft3d; ++i)
        rho_ind0[i] += overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]);


    /* step A */
    std::memcpy(rho_ind1,rho_ind0,n2ft3d*sizeof(double));
    if (restart)
    {
       mypneb->rr_SMul(scal1,vdielec,p);
       mypneb->rc_fft3d(p);
       mypneb->c_pack(0,p);

       mypneb->tcr_pack_iMul_unpack_fft(0,Gx,p,w_x);
       mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_y);
       mypneb->tcr_pack_iMul_unpack_fft(0,Gz,p,w_z);
       for (auto i=0; i<n2ft3d; ++i)
          rho_ind1[i] +=  overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]);
    }
    mycoulomb2->vcoulomb(rho_ind1,vdielec);
    mypneb->r_zero_ends(vdielec);
    mydplot->gcube_write("AP_rho_ind00.cube",-1,"SCF dielec function",rho_ind1);
    mydplot->gcube_write("AP_vdielec00.cube",-1,"SCF dielec function",vdielec);

    if (mypneb->d3db::parall->is_master() && debug)
       std::cout << "   it    rho_ind1*dv   rho_elc*dv    rho_ion*dv         epol     eelc-ind      eion-ind     epol-eold" << std::endl;

    /* iteration=0,1,... */
    double eold = 0.0;
    double epol = 0.5*mypneb->rr_dot(rho_ind1,vdielec)*dv;
    double alpha = 0.41;
    int it = 0;
    while (std::abs(epol-eold)>tole_pol)
    {
       ++it;
       /* step A */
       mypneb->rr_SMul(scal1,vdielec,p);
       mypneb->rc_fft3d(p);
       mypneb->c_pack(0,p);
       mypneb->tcr_pack_iMul_unpack_fft(0,Gx,p,w_x);
       mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_y);
       mypneb->tcr_pack_iMul_unpack_fft(0,Gz,p,w_z);
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
           rho_ind1[i] = (1.0-alpha)*rho_ind1[i] +
                       + alpha*(rho_ind0[i] + overfourpi*(w_x[i]*epsilon_x[i] + w_y[i]*epsilon_y[i] + w_z[i]*epsilon_z[i]));
       }
       mycoulomb2->vcoulomb(rho_ind1,vdielec);
       mypneb->r_zero_ends(vdielec);

       /* step B */


       eold = epol;
       epol = 0.5*mypneb->rr_dot(rho_ind1,vdielec)*dv;

       energy = mypneb->r_dsum(rho_ind1)*dv;
       double sumion = mypneb->r_dsum(rho_ion)*dv;
       double sumelc = mypneb->r_dsum(rho)*dv;
       double eelc = mypneb->rr_dot(rho,vdielec)*dv;
       double eion = mypneb->rr_dot(rho_ion,vdielec)*dv;
       if (mypneb->d3db::parall->is_master() && debug) {
          std::cout << Ifmt(5) << it << Efmt(15,6) << energy << " "
                                                   << sumelc << " "
                                                   << sumion << " "
                                                   << epol   << " "
                                                   << eelc   << " "
                                                   << eion   << " "
                                                   << epol-eold << std::endl;
       }
    }
    mydplot->gcube_write("AP_vdielecff.cube",-1,"SCF dielec function",vdielec);
    mydplot->gcube_write("AP_rho_indff.cube",-1,"SCF dielec function",rho_ind1);
    mydplot->gcube_write("AP_epsilon.cube",-1,"SCF dielec function",epsilon);
    mydplot->gcube_write("AP_epsilon_x.cube",-1,"SCF dielec function",epsilon_x);
    mydplot->gcube_write("AP_epsilon_y.cube",-1,"SCF dielec function",epsilon_y);
    mydplot->gcube_write("AP_epsilon_z.cube",-1,"SCF dielec function",epsilon_z);

    std::memcpy(sw,rho_ion,n2ft3d*sizeof(double));
    mydplot->gcube_write("AP_rho_ion.cube",-1,"SCF dielec function",sw);
    std::memcpy(sw,rho,n2ft3d*sizeof(double));
    mydplot->gcube_write("AP_rho_elc.cube",-1,"SCF dielec function",sw);

    /* */
    mypneb->rrr_Sum(vdielec,vloc,p);
    mypneb->rr_Sum(vh,p);
    mypneb->r_SMul(scal1,p);
    mypneb->rc_fft3d(p);
    mypneb->c_pack(0,p);
    mypneb->tcr_pack_iMul_unpack_fft(0,Gx,p,w_x);
    mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_y);
    mypneb->tcr_pack_iMul_unpack_fft(0,Gy,p,w_z);
    for (auto i=0; i<n2ft3d; ++i)
       vks_dielec[i] = -0.5*overfourpi*depsilon[i]*(w_x[i]*w_x[i] + w_y[i]*w_y[i] + w_z[i]*w_z[i]);

    mydplot->gcube_write("AP_vks_dielec.cube",-1,"SCF dielec function",vks_dielec);

    energy = mypneb->rr_dot(rho,vdielec);
    std::cout << "KS potential    rho*vdielec*dv = " << Efmt(15,6) << energy << std::endl;
    energy = mypneb->rr_dot(rho,vks_dielec);
    std::cout << "KS potential rho*vks_dielec*dv = " << Efmt(15,6) << energy << std::endl;

    return epol;

}


}
