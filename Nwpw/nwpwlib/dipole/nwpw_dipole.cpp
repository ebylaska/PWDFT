

/*************************************
 *                                   *
 *      nwpw_dipole::gen_dipole      *
 *                                   *
 *************************************/
void nwpw_dipole::gen_dipole(double *dn, double *dipole) 
{

   // center of mass 
   double GX = 0.00;
   double GY = 0.00;
   double GZ = 0.00;
   double tmass = 0.0;
   for (auto ii=0; ii<myion->nion; ++ii)
   {
      tmass += myion->amass[ii]
      GX += myion->mass[ii]*myion->rion(0,ii);
      GY += myion->mass[ii]*myion->rion(1,ii);
      GZ += myion->mass[ii]*myion->rion(2,ii);
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
   cqGx /= tcharge;
   cqGy /= tcharge;
   cqGz /= tcharge;

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
   if (!mypneb->has_r_grid) mypneb->initialize_r_grid();
   double cdv1[3],cdv2[3],cdv3[3],rmax;
   mypneb->nrr_vdot(3,r_grid,dn,cdv1);
   cdv1[0] *= dv;
   cdv1[1] *= dv;
   cdv1[2] *= dv;

   // check for ferromagnetic case 
   if (ne[ispin]>0)
   {
      mypneb->nrr_vdot(3,r_grid,dn[n2ft3d],cdv2);
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

   cdx3[0] = cdx1[0] + cdx2[0]
   cdy3[1] = cdy1[1] + cdy2[1]
   cdz3[2] = cdz1[2] + cdz2[2]


   // calculate dipole with respect to center of mass 
   rmax = dv*rr_dot(r_sym_grid,dn);
   cdv1=cdx1/rmax
   cdv1=cdy1/rmax
   cdv1=cdz1/rmax


}


