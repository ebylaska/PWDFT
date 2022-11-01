
#include        "iofmt.hpp"
#include	"nwpw_dipole.hpp"

namespace pwdft {


/* Constructor */ 

/*************************************
 *                                   *
 *      nwpw_dipole::nwpw_dipole     *
 *                                   *
 *************************************/
nwpw_dipole::nwpw_dipole(Ion *myion0, Pneb *mypneb0, Strfac *mystrfac0, Control2&)
{
   mypneb   = mypneb0;
   myion    = myion0;
   mystrfac = mystrfac0;

   ispin  = mypneb->ispin;
   n2ft3d = mypneb->n2ft3d;
   ne    = mypneb->ne;

   double omega = mypneb->lattice->omega();
   double scal1 = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
   dv = omega*scal1;

}

/* Object functions */

/*************************************
 *                                   *
 *      nwpw_dipole::gen_dipole      *
 *                                   *
 *************************************/
void nwpw_dipole::gen_dipole(const double *dn) 
{

   double *r_sym_grid = mypneb->r_nalloc(3);

   // center of mass 
   double GX = 0.00;
   double GY = 0.00;
   double GZ = 0.00;
   double tmass = 0.0;
   for (auto ii=0; ii<myion->nion; ++ii)
   {
      double massii = myion->mass[ii];
      tmass += massii;
      GX += massii*myion->rion(0,ii);
      GY += massii*myion->rion(1,ii);
      GZ += massii*myion->rion(2,ii);
   }
   GX /= tmass;
   GY /= tmass;
   GZ /= tmass;

   // crystal center of ionic charge
   int ncut = 20;
   int n0 = ncut-2;
   int n1 = ncut-2;
   int n2 = ncut-2;

   double x = n0*mypneb->lattice->unita(0,0);
   double y = n0*mypneb->lattice->unita(1,0);
   double z = n0*mypneb->lattice->unita(2,0);
   double rmax = std::sqrt(x*x + y*y + z*z);

   x = n1*mypneb->lattice->unita(0,1);
   y = n1*mypneb->lattice->unita(1,1);
   z = n1*mypneb->lattice->unita(2,1);
   double r = std::sqrt(x*x + y*y + z*z);
   if (r<rmax) rmax = r;

   x = n2*mypneb->lattice->unita(0,1);
   y = n2*mypneb->lattice->unita(1,1);
   z = n2*mypneb->lattice->unita(2,1);
   r = std::sqrt(x*x + y*y + z*z);
   if (r<rmax) rmax = r;

   double cqGX = 0.0;
   double cqGY = 0.0;
   double cqGZ = 0.0;
   double tcharge = 0.0;
   for (auto ii=0; ii<myion->nion; ++ii)
   {
      int ia = myion->katm[ii];
      double qion = myion->zv_psp[ia];

      for (auto n2=-ncut; n2<=ncut; ++n2)
      for (auto n1=-ncut; n1<=ncut; ++n1)
      for (auto n0=-ncut; n0<=ncut; ++n0)
      {
         x = myion->rion(0,ii)
           + n0*mypneb->lattice->unita(0,0)
           + n1*mypneb->lattice->unita(0,1)
           + n2*mypneb->lattice->unita(0,2);
         y = myion->rion(1,ii)
           + n0*mypneb->lattice->unita(1,0)
           + n1*mypneb->lattice->unita(1,1)
           + n2*mypneb->lattice->unita(1,2);
         z = myion->rion(2,ii)
           + n0*mypneb->lattice->unita(2,0)
           + n1*mypneb->lattice->unita(2,1)
           + n2*mypneb->lattice->unita(2,2);
         r = std::sqrt(x*x+y*y+z*z);
         if (r<=rmax)
         {
            cqGX += qion*x;
            cqGY += qion*y;
            cqGZ += qion*z;
            tcharge += qion;
         }
      }
   }
   cqGX /= tcharge;
   cqGY /= tcharge;
   cqGZ /= tcharge;

   // molecular center of ionic charge
   double qGX = 0.0;
   double qGY = 0.0;
   double qGZ = 0.0;
   tcharge = 0.0;
   for (auto ii=0; ii<myion->nion; ++ii)
   {
      int ia = myion->katm[ii];
      double qion = myion->zv_psp[ia];

      tcharge += qion;
      qGX += qion*myion->rion(0,ii);
      qGY += qion*myion->rion(1,ii);
      qGZ += qion*myion->rion(2,ii);
   }
   qGX /= tcharge;
   qGY /= tcharge;
   qGZ /= tcharge;


   // calculate the center of density 
   mypneb->generate_r_sym_grid(r_sym_grid);
   //double cdv1[3],cdv2[3],cdv3[3];
   mypneb->nrr_vdot(3,r_sym_grid,dn,mcdv1);
   mcdv1[0] *= dv;
   mcdv1[1] *= dv;
   mcdv1[2] *= dv;

   // check for ferromagnetic case 
   if (ne[ispin-1]!=0)
   {
      mypneb->nrr_vdot(3,r_sym_grid,dn+(ispin-1)*n2ft3d,mcdv2);
      mcdv2[0] *= dv;
      mcdv2[1] *= dv;
      mcdv2[2] *= dv;
   }
   else
   {
      mcdv2[0] = 0.0;
      mcdv2[1] = 0.0;
      mcdv2[2] = 0.0;
   }

   mcdv3[0] = mcdv1[0] + mcdv2[0];
   mcdv3[1] = mcdv1[1] + mcdv2[1];
   mcdv3[2] = mcdv1[2] + mcdv2[2];


   // calculate dipole with respect to center of mass 
   mypneb->generate_r_sym_mask(r_sym_grid);
   rmax = dv*mypneb->rr_dot(r_sym_grid,dn);
   mcdv1[0] /= rmax;
   mcdv1[1] /= rmax;
   mcdv1[2] /= rmax;
   if (ne[ispin-1]!=0)
   {
     double rmax2 = dv*mypneb->rr_dot(r_sym_grid,dn+(ispin-1)*n2ft3d);
     mcdv2[0] /= rmax2;
     mcdv2[1] /= rmax2;
     mcdv2[2] /= rmax2;
     rmax += rmax2;
   }
   mcdv3[0] /= rmax;
   mcdv3[1] /= rmax;
   mcdv3[2] /= rmax;


   // calculate dipole with respect to center of mass 
   double pcharge   = tcharge;;
   double ncharge   = ((double) (ne[0]+ne[ispin-1]));
   mdipole[0] = -ncharge*mcdv3[0] + pcharge*qGX - GX*(pcharge-ncharge);
   mdipole[1] = -ncharge*mcdv3[1] + pcharge*qGY - GY*(pcharge-ncharge);
   mdipole[2] = -ncharge*mcdv3[2] + pcharge*qGZ - GZ*(pcharge-ncharge);
   mqv1[0] = qGX;
   mqv1[1] = qGY;
   mqv1[2] = qGZ;

   mypneb->r_dealloc(r_sym_grid);
}

/***********************************************
 *                                             *
 *      nwpw_dipole::gen_molecular_dipole      *
 *                                             *
 ***********************************************/
void nwpw_dipole::gen_molecular_dipole(const double *dn, double *dipole) 
{
   double *r_sym_grid = mypneb->r_nalloc(3);

   // center of mass 
   double GX = 0.00;
   double GY = 0.00;
   double GZ = 0.00;
   double tmass = 0.0;
   for (auto ii=0; ii<myion->nion; ++ii)
   {
      double massii = myion->mass[ii];
      tmass += massii;
      GX += massii*myion->rion(0,ii);
      GY += massii*myion->rion(1,ii);
      GZ += massii*myion->rion(2,ii);
   }
   GX /= tmass;
   GY /= tmass;
   GZ /= tmass;


   // molecular center of ionic charge
   double qGX = 0.0;
   double qGY = 0.0;
   double qGZ = 0.0;
   double tcharge = 0.0;
   for (auto ii=0; ii<myion->nion; ++ii)
   {
      int ia = myion->katm[ii];
      double qion = myion->zv_psp[ia];

      tcharge += qion;
      qGX += qion*myion->rion(0,ii);
      qGY += qion*myion->rion(1,ii);
      qGZ += qion*myion->rion(2,ii);
   }
   qGX /= tcharge;
   qGY /= tcharge;
   qGZ /= tcharge;

   // calculate the center of density 
   mypneb->generate_r_sym_grid(r_sym_grid);
   double cdv1[3],cdv2[3],cdv3[3];
   mypneb->nrr_vdot(3,r_sym_grid,dn,cdv1);
   cdv1[0] *= dv;
   cdv1[1] *= dv;
   cdv1[2] *= dv;

   // check for ferromagnetic case 
   if (ne[ispin-1]>0)
   {
      mypneb->nrr_vdot(3,r_sym_grid,dn+n2ft3d,cdv2);
      cdv2[0] *= dv;
      cdv2[1] *= dv;
      cdv2[2] *= dv;
   }
   else
   {
      cdv2[0] = 0.0;
      cdv2[1] = 0.0;
      cdv2[2] = 0.0;
   }

   cdv3[0] = cdv1[0] + cdv2[0];
   cdv3[1] = cdv1[1] + cdv2[1];
   cdv3[2] = cdv1[2] + cdv2[2];

   // calculate dipole with respect to center of mass 
   mypneb->generate_r_sym_mask(r_sym_grid);
   double rmax = dv*mypneb->rr_dot(r_sym_grid,dn);
   cdv1[0] /= rmax;
   cdv1[1] /= rmax;
   cdv1[2] /= rmax;
   if (ne[ispin-1]>0)
   {
     int rmax2 = dv*mypneb->rr_dot(r_sym_grid,dn+(ispin-1)*n2ft3d);
     cdv2[0] /= rmax2;
     cdv2[1] /= rmax2;
     cdv2[2] /= rmax2;
     std::cout << "rmax2=" << rmax2 << std::endl;
     rmax += rmax2;
   }
   cdv3[0] /= rmax;
   cdv3[1] /= rmax;
   cdv3[2] /= rmax;

   // calculate dipole with respect to center of mass 
   double pcharge   = tcharge;;
   double ncharge   = ((double) (ne[0]+ne[ispin-1]));
   dipole[0] = -ncharge*cdv3[0] + pcharge*qGX - GX*(pcharge-ncharge);
   dipole[1] = -ncharge*cdv3[1] + pcharge*qGY - GY*(pcharge-ncharge);
   dipole[2] = -ncharge*cdv3[2] + pcharge*qGZ - GZ*(pcharge-ncharge);

   mypneb->r_dealloc(r_sym_grid);
}

/***********************************************
 *                                             *
 *        nwpw_dipole::gen_Resta_dipole        *
 *                                             *
 ***********************************************/
void nwpw_dipole::gen_Resta_dipole(const double *psi, double *dipole) 
{

   // allocate psi_r and translate matrices
   double *psi_r  = mypneb->h_allocate();
   double *psi_r2 = mypneb->h_allocate();


   mypneb->h_deallocate(psi_r);
   mypneb->h_deallocate(psi_r2);
}

/***********************************************
 *                                             *
 *        nwpw_dipole::shortprint_dipole       *
 *                                             *
 ***********************************************/

std::string nwpw_dipole::shortprint_dipole()
{
   double mu = std::sqrt(mdipole[0]*mdipole[0] + mdipole[1]*mdipole[1] + mdipole[2]*mdipole[2]);
   std::stringstream stream;

   stream << "== Center of Charge ==" << std::endl << std::endl;
   stream << "spin up    = (" << Ffmt(10,4) << mcdv1[0] << " " 
                              << Ffmt(10,4) << mcdv1[1] << " " 
                              << Ffmt(10,4) << mcdv1[2] << " )" << std::endl;
   stream << "spin down  = (" << Ffmt(10,4) << mcdv2[0] << " " 
                              << Ffmt(10,4) << mcdv2[1] << " " 
                              << Ffmt(10,4) << mcdv2[2] << " )" << std::endl;
   stream << "     total = (" << Ffmt(10,4) << mcdv3[0] << " " 
                              << Ffmt(10,4) << mcdv3[1] << " " 
                              << Ffmt(10,4) << mcdv3[2] << " )" << std::endl;
   stream << "ionic      = (" << Ffmt(10,4) << mqv1[0] << " " 
                              << Ffmt(10,4) << mqv1[1] << " " 
                              << Ffmt(10,4) << mqv1[2] << " )" << std::endl;
   stream << std::endl;
   stream << "== Molecular Dipole wrt Center of Mass ==" << std::endl << std::endl;
   stream << "mu   = (" << Ffmt(10,4) << mdipole[0] << " " 
                        << Ffmt(10,4) << mdipole[1] << " " 
                        << Ffmt(10,4) << mdipole[2] << " ) au" << std::endl;
   stream << "|mu| =  " << Ffmt(10,4) <<  mu           << " au ( " 
                        << Ffmt(10,4) <<  mu*autoDebye << " Debye )"  
                        << std::endl;
   stream << std::endl;

   return stream.str();
}


}
