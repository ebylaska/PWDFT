/* Pseudopotential.C -
   Author - Eric Bylaska
*/

#include        <iostream>
#include	<cstring>
#include        <cmath>

#include	"nwpw_timing.hpp"
#include	"util.hpp"
#include        "Psp1d_Hamann.hpp"
#include        "gdevice.hpp"


#include        "blas.h"

#include	"compressed_io.hpp"
#include	"Pseudopotential.hpp"

/*******************************************
 *                                         *
 *              vpp_read_header            *
 *                                         *
 *******************************************/
static bool vpp_read_header(char *fname,
                            char *comment, int *psp_type, int *version, int nfft[], double unita[],
                            char *atom, double *amass, double *zv)
{
   int i,ifound;

   ifound = cfileexists(fname);


   if (ifound>0)
   {
      openfile(5,fname,"r");
      cread(5,comment,80);
      comment[79] = '\0';
      i = 78;
      while (comment[i] == ' ')
        comment[i--] = '\0';

      iread(5,psp_type,1);
      iread(5,version,1);
      iread(5,nfft,3);
      dread(5,unita,9);
      cread(5,atom,2);
      dread(5,amass,1);
      dread(5,zv,1);
      closefile(5);
   }

   return (ifound>0);
}

/*****************************************************
 *                                                   *
 *                vpp_formatter_check                *
 *                                                   *
 *****************************************************/

static bool vpp_formatter_check(PGrid *mygrid, char *fname)
{
   char comment[80],atom[2];
   int psp_type,version,nfft[3];
   double unita[9],amass,zv;
   bool reformat;
   double tol=1.0e-9;

   reformat = true;
   if (vpp_read_header(fname,
                            comment, &psp_type, &version, nfft, unita,
                            atom, &amass, &zv))
   {
      reformat = false;
      for (auto i=0; i<9; ++i)
         reformat = reformat || (fabs(mygrid->lattice->unita1d(i)-unita[i])>tol);
      reformat = reformat || (mygrid->nx!=nfft[0]);
      reformat = reformat || (mygrid->ny!=nfft[1]);
      reformat = reformat || (mygrid->nz!=nfft[2]);
      //reformat = reformat || (control.pversion!=version);
   }
   return reformat;
}






/*******************************************
 *                                         *
 *           Multiply_Gijl_sw1             *
 *                                         *
 *******************************************/
static void Multiply_Gijl_sw1(int nn,
                       const int nprj,
                       const int nmax,
                       const int lmax,
                       int *n_prj,
                       int *l_prj,
                       int *m_prj,
                       double *G,
                       double *sw1,
                       double *sw2)
{
   int a,b,na,nb;
   int one  = 1;
   int zero = 0;
   double rzero = 0.0;
   int nmax2 = nmax*nmax;
   int nnn = nn*nprj;
   int nna = nn;

   DCOPY_PWDFT(nnn,&rzero,zero,sw2,one);

   for (b=0; b<nprj; ++b)
      for (a=0; a<nprj; ++a)
         if ((l_prj[a]==l_prj[b]) && (m_prj[a]==m_prj[b]))
         {
            na = n_prj[a]-1;
            nb = n_prj[b]-1;
            //na = n_prj[a];
            //nb = n_prj[b];

           DAXPY_PWDFT(nna,
                   (G[nb + na*nmax + nmax2*l_prj[a]]),
                   &(sw1[a*nn]),one,
                   &(sw2[b*nn]),one);
          }
}

/*******************************************
 *                                         *
 *                vpp_read                 *
 *                                         *
 *******************************************/
static void vpp_read(PGrid *mygrid,
                     char *fname,
                     char *comment,
                     int *psp_type,
                     int *version,
                     int *nfft,
                     double *unita,
                     char   *atom,
                     double *amass,
                     double *zv,
                     int *lmmax,
                     int *lmax,
                     int *locp,
                     int *nmax,
                     double **rc,
                     int *nprj,
                     int **n_projector,
                     int **l_projector,
                     int **m_projector,
                     int **b_projector,
                     double **Gijl,
                     bool   *semicore,
                     double *rcore,
                     double **ncore,
                     double *vl,
                     double **vnl)
{
   int i,nn;
   double   *tmp2,*prj;
   Parallel *parall = mygrid->parall;

   if (parall->is_master())
   {
      std::cout << "reading formatted psp filename: " << fname << std::endl;
      openfile(5,fname,"r");
      cread(5,comment,80);
      comment[79] = '\0';
      i = 78;
      while (comment[i] == ' ')
        comment[i--] = '\0';

      iread(5,psp_type,1);
      iread(5,version,1);
      iread(5,nfft,3);
      dread(5,unita,9);
      cread(5,atom,2);
      dread(5,amass,1);
      dread(5,zv,1);
      iread(5,lmax,1);
      iread(5,locp,1);
      iread(5,nmax,1);
   }
   parall->Brdcst_cValues(0,0,80,comment);
   parall->Brdcst_iValue(0,0,psp_type);
   parall->Brdcst_iValue(0,0,version);
   parall->Brdcst_iValues(0,0,3,nfft);
   parall->Brdcst_Values(0,0,9,unita);
   parall->Brdcst_cValues(0,0,2,atom);
   parall->Brdcst_Values(0,0,1,amass);
   parall->Brdcst_Values(0,0,1,zv);
   parall->Brdcst_iValue(0,0,lmax);
   parall->Brdcst_iValue(0,0,locp);
   parall->Brdcst_iValue(0,0,nmax);
   *lmmax=((*lmax)+1)*((*lmax)+1) - (2*(*locp)+1);

   *rc = new double[*lmax+1];
   if (parall->is_master())
   {
      dread(5,*rc,*lmax+1);
      iread(5,nprj,1);
   }
   parall->Brdcst_Values(0,0,*lmax+1,*rc);
   parall->Brdcst_iValue(0,0,nprj);
   if (*nprj > 0)
   {
      *n_projector = new int[*nprj];
      *l_projector = new int[*nprj];
      *m_projector = new int[*nprj];
      *b_projector = new int[*nprj];
      if (parall->is_master())
      {
         iread(5,*n_projector,*nprj);
         iread(5,*l_projector,*nprj);
         iread(5,*m_projector,*nprj);
         iread(5,*b_projector,*nprj);
      }
      parall->Brdcst_iValues(0,0,*nprj,*n_projector);
      parall->Brdcst_iValues(0,0,*nprj,*l_projector);
      parall->Brdcst_iValues(0,0,*nprj,*m_projector);
      parall->Brdcst_iValues(0,0,*nprj,*b_projector);

      nn = (*nmax)*(*nmax)*(*lmax+1);
      *Gijl = new double[nn];
      if (parall->is_master())
      {
         dread(5,*Gijl,nn);
      }
      parall->Brdcst_Values(0,0,nn,*Gijl);
   }
   if (parall->is_master()) dread(5,rcore,1);
   parall->Brdcst_Values(0,0,1,rcore);
   if (*rcore > 0.0)
      *semicore = true;
   else
      *semicore = false;


   /* readin vl 3d block */
   tmp2 = new double [mygrid->nfft3d];
   mygrid->t_read(5,tmp2,-1);
   mygrid->t_pack(0,tmp2);
   mygrid->tt_pack_copy(0,tmp2,vl);

   /* reading vnl 3d block */
   if (*nprj > 0)
   {
      *vnl = new double[(*nprj)*(mygrid->npack(1))];
      prj = *vnl;
      for (i=0; i<(*nprj); ++i)
      {
         mygrid->t_read(5,tmp2,-1);
         mygrid->t_pack(1,tmp2);
         mygrid->tt_pack_copy(1,tmp2,&prj[i*mygrid->npack(1)]);
      }
   }
   if (*semicore)
   {
      nn     = 5*mygrid->npack(0);
      *ncore = new double[nn];
      prj    = *ncore;

      mygrid->t_read(5,tmp2,-1);
      mygrid->t_pack(0,tmp2);
      mygrid->tt_pack_copy(0,tmp2,prj);

      mygrid->t_read(5,tmp2,-1);
      mygrid->t_pack(0,tmp2);
      mygrid->tt_pack_copy(0,tmp2,&prj[2*mygrid->npack(0)]);

      mygrid->t_read(5,tmp2,-1);
      mygrid->t_pack(0,tmp2);
      mygrid->tt_pack_copy(0,tmp2,&prj[3*mygrid->npack(0)]);

      mygrid->t_read(5,tmp2,-1);
      mygrid->t_pack(0,tmp2);
      mygrid->tt_pack_copy(0,tmp2,&prj[4*mygrid->npack(0)]);
   }

   delete [] tmp2;

   if (parall->is_master()) closefile(5);
}



/*******************************************
 *                                         *
 *                vpp_write                *
 *                                         *
 *******************************************/
static void vpp_write(PGrid *mygrid,
                     char *fname,
                     char *comment,
                     int psp_type,
                     int version,
                     int *nfft,
                     double *unita,
                     char   *atom,
                     double amass,
                     double zv,
                     int lmmax,
                     int lmax,
                     int locp,
                     int nmax,
                     double *rc,
                     int nprj,
                     int *n_projector,
                     int *l_projector,
                     int *m_projector,
                     int *b_projector,
                     double *Gijl,
                     bool   semicore,
                     double rcore,
                     double *ncore,
                     double *vl,
                     double *vnl)
{
   int i,nn;
   double   *tmp2,*prj;
   Parallel *parall = mygrid->parall;

   if (parall->is_master())
   {
      openfile(6,fname,"w");
      comment[79] = '\0';
      cwrite(6,comment,80);

      iwrite(6,&psp_type,1);
      iwrite(6,&version,1);
      iwrite(6,nfft,3);
      dwrite(6,unita,9);
      cwrite(6,atom,2);
      dwrite(6,&amass,1);
      dwrite(6,&zv,1);
      iwrite(6,&lmax,1);
      iwrite(6,&locp,1);
      iwrite(6,&nmax,1);
   }

   if (parall->is_master())
   {
      dwrite(6,rc,lmax+1);
      iwrite(6,&nprj,1);
   }
   if (nprj > 0)
   {
      if (parall->is_master())
      {
         iwrite(6,n_projector,nprj);
         iwrite(6,l_projector,nprj);
         iwrite(6,m_projector,nprj);
         iwrite(6,b_projector,nprj);
      }
      nn = (nmax)*(nmax)*(lmax+1);
      if (parall->is_master())
      {
         dwrite(6,Gijl,nn);
      }
   }
   if (parall->is_master()) dwrite(6,&rcore,1);



   /* readin vl 3d block */
   tmp2 = new double [mygrid->nfft3d];
   mygrid->tt_pack_copy(0,vl,tmp2);
   mygrid->t_unpack(0,tmp2);
   mygrid->t_write(6,tmp2,-1);



   /* reading vnl 3d block */
   if (nprj > 0)
   {
      prj = vnl;
      for (i=0; i<(nprj); ++i)
      {
         mygrid->tt_pack_copy(1,&prj[i*mygrid->npack(1)],tmp2);
         mygrid->t_unpack(1,tmp2);
         mygrid->t_write(6,tmp2,-1);
      }
   }
   if (semicore)
   {
      prj    = ncore;

      mygrid->tt_pack_copy(0,prj,tmp2);
      mygrid->t_unpack(0,tmp2);
      mygrid->t_write(6,tmp2,-1);


      mygrid->tt_pack_copy(0,&prj[2*mygrid->npack(0)],tmp2);
      mygrid->t_unpack(0,tmp2);
      mygrid->t_write(6,tmp2,-1);

      mygrid->tt_pack_copy(0,&prj[3*mygrid->npack(0)],tmp2);
      mygrid->t_unpack(0,tmp2);
      mygrid->t_write(6,tmp2,-1);

      mygrid->tt_pack_copy(0,&prj[4*mygrid->npack(0)],tmp2);
      mygrid->t_unpack(0,tmp2);
      mygrid->t_write(6,tmp2,-1);
   }

   delete [] tmp2;

   if (parall->is_master()) closefile(6);
}



/*******************************************
 *                                         *
 *            semicore_check               *
 *                                         *
 *******************************************/
static double semicore_check(PGrid *mygrid, bool semicore, double rcore, double *ncore)
{
   double sum = 0.0;
   if (semicore)
   {
      double omega = mygrid->lattice->omega();
      double scal1 = 1.0/((double) ((mygrid->nx)*(mygrid->ny)*(mygrid->nz)));
      //double scal2 = 1.0/lattice_omega();
      //double dv    = lattice_omega()*scal1;
      double scal2 = 1.0/omega;
      double dv    = omega*scal1;
      double *tmp  = mygrid->r_alloc();

      /* put sqrt(core-density) at atom position */
      mygrid->tc_pack_copy(0,ncore,tmp);
      mygrid->c_SMul(0,scal2,tmp);

      /* Put put tmp into real space */
      mygrid->c_unpack(0,tmp);
      mygrid->cr_fft3d(tmp);

      /*  square it  */
      mygrid->r_sqr(tmp);

      /* integrate it */
      sum = mygrid->r_dsum(tmp) * dv;

      mygrid->r_dealloc(tmp);
   }
   return sum;
}

static int convert_psp_type(char *test)
{
   int psp_type = 0;
   if (test[0]=='0') psp_type = 0;
   if (test[0]=='1') psp_type = 1;
   if (test[0]=='2') psp_type = 2;
   if (test[0]=='3') psp_type = 3;
   if (test[0]=='4') psp_type = 4;
   if (test[0]=='5') psp_type = 5;
   if (test[0]=='6') psp_type = 6;
   if (test[0]=='7') psp_type = 7;
   if (test[0]=='8') psp_type = 8;
   if (test[0]=='9') psp_type = 9;

   return psp_type;
}

/*******************************************
 *                                         *
 *            vpp_get_psp_type             *
 *                                         *
 *******************************************/
static int vpp_get_psp_type(Parallel *myparall, char *pspname)
{
   int psp_type;
   char atom[2];

   if (myparall->is_master())
   {
      FILE *fp = std::fopen(pspname,"r");
      std::fscanf(fp,"%s",atom);
      psp_type = convert_psp_type(atom);
      //if (psp_type>0) std::fscanf(fp,"%s",atom);
      fclose(fp);
   }
   myparall->Brdcst_iValue(0,0,&psp_type);

   return psp_type;
}



/*******************************************
 *                                         *
 *                vpp_generate             *
 *                                         *
 *******************************************/
static void vpp_generate(PGrid *mygrid,
                     char *pspname,
                     char *fname,
                     char *comment,
                     int *psp_type,
                     int *version,
                     int *nfft,
                     double *unita,
                     char   *atom,
                     double *amass,
                     double *zv,
                     int *lmmax,
                     int *lmax,
                     int *locp,
                     int *nmax,
                     double **rc,
                     int *nprj,
                     int **n_projector,
                     int **l_projector,
                     int **m_projector,
                     int **b_projector,
                     double **Gijl,
                     bool   *semicore,
                     double *rcore,
                     double **ncore,
                     double *vl,
                     double **vnl)
{
   int i,nn;
   double   *tmp2,*prj;
   Parallel *myparall = mygrid->parall;

   *psp_type = vpp_get_psp_type(myparall, pspname);

   //std::cout << "in vpp_generate, psp_type=" << *psp_type << std::endl;
   //std::cout << "in vpp_generate" << std::endl;

   if ((*psp_type==0) || (*psp_type==9))
   {
      int nray = mygrid->n_ray();
      Psp1d_Hamann psp1d(myparall,pspname);
      //std::cout << "in vpp_generate Hamann psp1d.psp_type=" << psp1d.psp_type << std::endl;
      //std::cout << "in vpp_generate Hamann nray=" << nray << std::endl;
      //std::cout << "in vpp_generate Hamann Gmax=" << mygrid->Gmax_ray() << std::endl;
      //std::cout << "in vpp_generate Hamann Gmin=" << mygrid->Gmin_ray() << std::endl;
      //std::cout << "in vpp_generate Hamann dGmin=" << mygrid->dGmin_ray() << std::endl;
      //std::cout << "in vpp_generate Hamann lmax=" << psp1d.lmax << std::endl;
      //std::cout << "in vpp_generate Hamann locp=" << psp1d.locp << std::endl;
      //std::cout << "in vpp_generate Hamann nmax=" << psp1d.nmax << std::endl;
      //std::cout << "in vpp_generate Hamann nprj=" << psp1d.nprj << std::endl;
      //std::cout << "in vpp_generate Hamann n_extra=" << psp1d.n_extra << std::endl;

      nfft[0] = mygrid->nx;
      nfft[1] = mygrid->ny;
      nfft[2] = mygrid->nz;
      for (auto i=0; i<9; ++i)
         unita[i] = mygrid->lattice->unita1d(i);

      atom[0] = psp1d.atom[0];
      atom[1] = psp1d.atom[1];
      *amass  = psp1d.amass;
      *zv     = psp1d.zv;
      for (auto i=0; i<80; ++i)
          comment[i] = psp1d.comment[i];

      *psp_type = psp1d.psp_type;
      *version  = psp1d.version;
      *lmax     = psp1d.lmax;
      *locp     = psp1d.locp;
      *nmax     = psp1d.nmax;
      *lmmax=((*lmax)+1)*((*lmax)+1) - (2*(*locp)+1);

      *nprj     = psp1d.nprj;
      *semicore = psp1d.semicore;
      *rcore    = psp1d.rcore;

      *rc = new double[*lmax+1];
      for (auto l=0; l<=(*lmax); ++l)
         (*rc)[l] = psp1d.rc[l];


      /* allocate Gijl and copy from psp1d */
      int nn = (psp1d.nmax)*(psp1d.nmax)*(psp1d.lmax+1);
      *Gijl = new double[nn];
      for (auto l=0; l<nn; ++l)
      {
          //std::cout << "in vpp_generate Hamann l=" << l << " norm=" << psp1d.vnlnrm[l] << std::endl;
         (*Gijl)[l] = psp1d.vnlnrm[l];
      }

      /* allocate n_projector, l_projector, m_projector, and b_projector and copy from psp1d */
      if (psp1d.nprj>0)
      {
         *n_projector = new int[psp1d.nprj];
         *l_projector = new int[psp1d.nprj];
         *m_projector = new int[psp1d.nprj];
         *b_projector = new int[psp1d.nprj];

         for (auto l=0; l<psp1d.nprj; ++l)
         {
            (*n_projector)[l] = psp1d.n_prj[l];
            (*l_projector)[l] = psp1d.l_prj[l];
            (*m_projector)[l] = psp1d.m_prj[l];
            (*b_projector)[l] = psp1d.b_prj[l];
         }
      }


      /*  allocate and generate ray formatted grids */
      double *G_ray = mygrid->generate_G_ray();
      double *vl_ray  = new double [nray];
      double *vnl_ray = new double [(psp1d.lmax+1+psp1d.n_extra)*nray];
      double *rho_sc_k_ray = new double [2*nray];
      psp1d.vpp_generate_ray(myparall,nray,G_ray,vl_ray,vnl_ray,rho_sc_k_ray);


      /* filter the ray formatted grids */
      double ecut = mygrid->lattice->ecut();
      double wcut = mygrid->lattice->wcut();
      util_filter(nray,G_ray,ecut,vl_ray);
      for (auto l=0; l<(psp1d.lmax+1+psp1d.n_extra); ++l)
         util_filter(nray,G_ray,wcut,&(vnl_ray[l*nray]));
      if (*semicore)
      {
         util_filter(nray,G_ray,ecut,rho_sc_k_ray);
         util_filter(nray,G_ray,ecut,&(rho_sc_k_ray[nray]));
      }

      /* allocate vnl and ncore  generate formated grids */
      *vnl = new double[(psp1d.nprj)*(mygrid->npack(1))];
      if (*semicore)
         *ncore = new double[5*mygrid->npack(0)];

      /*  generate formatted grids using splines */
      psp1d.vpp_generate_spline(mygrid,nray,G_ray,vl_ray,vnl_ray,rho_sc_k_ray,
                                vl,*vnl,*ncore);

      /* deallocate ray formatted grids */
      delete [] rho_sc_k_ray;
      delete [] vnl_ray;
      delete [] vl_ray;
      delete [] G_ray;


   }
   else
   {
      std::cout << "in vpp_generate Not finished, psp_type = " << *psp_type <<  std::endl;
   }
}



/* Constructors */

/*******************************************
 *                                         *
 *     Pseudopotential::Pseudopotential    *
 *                                         *
 *******************************************/
Pseudopotential::Pseudopotential(Ion *myionin, Pneb *mypnebin, Strfac *mystrfacin, Control2& control)
{
   int ia,version,nfft[3];
   int *n_ptr,*l_ptr,*m_ptr,*b_ptr;
   double *rc_ptr,*G_ptr,*vnl_ptr,*ncore_ptr;
   double unita[9];
   char fname[256],pspname[256],aname[2];

   myion    = myionin;
   mypneb   = mypnebin;
   mystrfac = mystrfacin;

   npsp = myion->nkatm;
   nprj_max = 0;

   psp_type = new int[npsp];
   lmax     = new int[npsp];
   lmmax    = new int[npsp];
   locp     = new int[npsp];
   nmax     = new int[npsp];
   nprj     = new int[npsp];
   semicore = new bool[npsp+1];

   n_projector = new int* [npsp];
   l_projector = new int* [npsp];
   m_projector = new int* [npsp];
   b_projector = new int* [npsp];

   zv          = new double[npsp];
   amass       = new double[npsp];
   rcore       = new double[npsp];
   ncore_sum   = new double[npsp];
   rc          = new double* [npsp];
   vl          = new double* [npsp];
   for (ia=0; ia<npsp; ++ia)
      vl[ia] = new double [mypneb->npack(0)];
   Gijl        = new double* [npsp];
   vnl         = new double* [npsp];
   ncore_atom  = new double* [npsp];
   semicore_density = mypneb->r_alloc();

   comment  = new  char* [npsp];
   for (ia=0; ia<npsp; ++ia) comment[ia] = new char[80];

   semicore[npsp] = false;
   for (ia=0; ia<npsp; ++ia)
   {
      strcpy(fname,myion->atom(ia));
      strcat(fname,".vpp");
      control.add_permanent_dir(fname);
      if (vpp_formatter_check(mypneb,fname))
      {
         strcpy(pspname,myion->atom(ia));
         strcat(pspname,".psp");
         control.add_permanent_dir(pspname);
         vpp_generate(mypneb,
		      pspname,fname,
		      comment[ia],&psp_type[ia],&version,nfft,unita,aname,
		      &amass[ia],&zv[ia],&lmmax[ia],&lmax[ia],&locp[ia],&nmax[ia],
		      &rc_ptr,&nprj[ia],&n_ptr,&l_ptr,&m_ptr,
		      &b_ptr,&G_ptr,&semicore[ia],&rcore[ia],
		      &ncore_ptr,vl[ia],&vnl_ptr);


/* *** still a problem with vpp_write in semicore stuff ****
         vpp_write(mypneb,
                  fname,
                  comment[ia],psp_type[ia],version,nfft,unita,aname,
                  amass[ia],zv[ia],lmmax[ia],lmax[ia],locp[ia],nmax[ia],
                  rc_ptr,nprj[ia],n_ptr,l_ptr,m_ptr,b_ptr,G_ptr,semicore[ia],rcore[ia],
                  ncore_ptr,vl[ia],vnl_ptr);
*/

      }
      else
      {
         // ******** debug *********
         /*
         strcpy(pspname,myion->atom(ia));
         strcat(pspname,".psp");
         control.add_permanent_dir(pspname);
         vpp_generate(mypneb,
                  pspname,fname,
                  comment[ia],&psp_type[ia],&version,nfft,unita,aname,
                  &amass[ia],&zv[ia],&lmmax[ia],&lmax[ia],&locp[ia],&nmax[ia],
                  &rc_ptr,&nprj[ia],&n_ptr,&l_ptr,&m_ptr,&b_ptr,&G_ptr,&semicore[ia],&rcore[ia],
                  &ncore_ptr,vl[ia],&vnl_ptr);

         std::cout << "vpp_generate Gptr l= 0, norm=" << G_ptr[0] << std::endl;;
         std::cout << "vpp_generate Gptr l= 1, norm=" << G_ptr[1] << std::endl;;
         std::cout << "vpp_generate Gptr l= 2, norm=" << G_ptr[2] << std::endl;;
          std::cout << "vpp_generate vl[0] = " << vl[ia][0] << std::endl;
          std::cout << "vpp_generate vl[1] = " << vl[ia][1] << std::endl;
          std::cout << "vpp_generate vl[431] = " << vl[ia][431] << std::endl;
          std::cout << "vpp_generate vl[9431] = " << vl[ia][9431] << std::endl;
          std::cout << "vpp_generate vnl[431] = " << vnl_ptr[431] << std::endl;
          std::cout << "vpp_generate vnl[31+2*npack1] = " << vnl_ptr[31+2*mypneb->npack(1)] << std::endl;
          std::cout << std::endl;
          delete [] rc_ptr;
          delete [] n_ptr;
          delete [] l_ptr;
          delete [] m_ptr;
          delete [] b_ptr;
          delete [] G_ptr;
          delete [] vnl_ptr;
          if (semicore[ia]) delete [] ncore_ptr;
         // ******** debug *********
         std::cout << "VPP READ = " << fname << std::endl;
         */

         vpp_read(mypneb,
                  fname,
                  comment[ia],&psp_type[ia],&version,nfft,unita,aname,
                  &amass[ia],&zv[ia],&lmmax[ia],&lmax[ia],&locp[ia],&nmax[ia],
                  &rc_ptr,&nprj[ia],&n_ptr,&l_ptr,&m_ptr,
		  &b_ptr,&G_ptr,&semicore[ia],&rcore[ia],
                  &ncore_ptr,vl[ia],&vnl_ptr);
      }
       /*
         std::cout << "vpp_read Gptr l= 0, norm=" << G_ptr[0] << std::endl;;
         std::cout << "vpp_read Gptr l= 1, norm=" << G_ptr[1] << std::endl;;
         std::cout << "vpp_read Gptr l= 2, norm=" << G_ptr[2] << std::endl;;
          std::cout << "vpp_read vl[0] = " << vl[ia][0] << std::endl;
          std::cout << "vpp_read vl[1] = " << vl[ia][1] << std::endl;
          std::cout << "vpp_read vl[431] = " << vl[ia][431] << std::endl;
          std::cout << "vpp_read vl[9431] = " << vl[ia][9431] << std::endl;
          std::cout << "vpp_read vnl[0] = " << vnl_ptr[0] << std::endl;
          std::cout << "vpp_read vnl[1] = " << vnl_ptr[1] << std::endl;
          std::cout << "vpp_read vnl[431] = " << vnl_ptr[431] << std::endl;
          std::cout << "vpp_read vnl[31+2*npack1] = " << vnl_ptr[31+2*mypneb->npack(1)] << std::endl;

          std::cout << "VNLS = " << mypneb->tt_pack_dot(1,vnl_ptr,vnl_ptr) << std::endl;
          std::cout << "VNLPX = " << mypneb->tt_pack_dot(1,&(vnl_ptr[mypneb->npack(1)]), &(vnl_ptr[mypneb->npack(1)]) ) << std::endl;
      */

      rc[ia]          = rc_ptr;
      n_projector[ia] = n_ptr;
      l_projector[ia] = l_ptr;
      m_projector[ia] = m_ptr;
      b_projector[ia] = b_ptr;
      Gijl[ia]        = G_ptr;
      vnl[ia]         = vnl_ptr;
      if (nprj[ia]>nprj_max) nprj_max = nprj[ia];

      if (semicore[ia])
      {
         ncore_atom[ia] = ncore_ptr;
         ncore_sum[ia]  = semicore_check(mypneb,semicore[ia],rcore[ia],ncore_atom[ia]);
         semicore[npsp] = true;
      }
   }

   /* define the maximum number of projectors  */
   nprj_max *= 10;
   // nprj_max = 0;
   // for (auto ii=0; ii < myion->nion; ++ii)
   //     nprj_max += nprj[myion->katm[ii]];

}
/*******************************************
 *                                         *
 *     Pseudopotential::v_nonlocal         *
 *                                         *
 *******************************************/
void Pseudopotential::v_nonlocal(double *psi, double *Hpsi)
{

   nwpw_timing_function ftimer(6);
   bool done;
   int ii,ia,l,nshift0,sd_function,i;
   int jj,ll,jstart,jend,nprjall;
   double *exi;
   double *prjtmp,*sw1,*sw2,*prj,*vnlprj;
   Parallel *parall;
   double omega = mypneb->lattice->omega();
   //double scal = 1.0/lattice_omega();
   double scal = 1.0/omega;
   int one=1;
   int ntmp,nshift,nn;
   double rone  = 1.0;
   double rmone = -1.0;

   nn     = mypneb->neq[0]+mypneb->neq[1];
   nshift0 = mypneb->npack(1);
   nshift = 2*mypneb->npack(1);
   exi    = new double[nshift];
   prjtmp = new double[nprj_max*nshift];
   sw1    = new double[nn*nprj_max];
   sw2    = new double[nn*nprj_max];

   parall = mypneb->d3db::parall;

#if 1

   //Copy psi to device
   gdevice_psi_copy_host2gpu(nshift0,nn,psi);
   gdevice_hpsi_copy_host2gpu(nshift0,nn,Hpsi);
   
   ii = 0;
   while (ii<(myion->nion))
   {
      ia      = myion->katm[ii];
      nprjall = 0;
      jstart  = ii;
      done = false;
      while (!done)
      {
         //generate projectors
         if (nprj[ia]>0)
         {
            mystrfac->strfac_pack(1,ii,exi);
            for (l=0; l<nprj[ia]; ++l)
            {
               sd_function = !(l_projector[ia][l] & 1);
               prj = &(prjtmp[(l+nprjall)*nshift]);
               vnlprj = &(vnl[ia][l*nshift0]);
               if (sd_function)
                  mypneb->tcc_Mul( 1,vnlprj,exi,prj);
               else
                  mypneb->tcc_iMul(1,vnlprj,exi,prj);
            }
            nprjall += nprj[ia];
         }
         ++ii;
         if (ii<(myion->nion))
         {
            ia = myion->katm[ii];
            done = ((nprjall+nprj[ia]) > nprj_max);
         }
         else
         {
            done = true;
         }
      }
      jend = ii;
      mypneb->cc_pack_inprjdot(1,nn,nprjall,psi,prjtmp,sw1);
      parall->Vector_SumAll(1,nn*nprjall,sw1);

      /* sw2 = Gijl*sw1 */
      ll = 0;
      for (jj=jstart; jj<jend; ++jj)
      {
         ia = myion->katm[jj];
         if (nprj[ia]>0)
         {
            Multiply_Gijl_sw1(nn,nprj[ia],nmax[ia],lmax[ia],
                              n_projector[ia],l_projector[ia],m_projector[ia],
                              Gijl[ia],&(sw1[ll*nn]),&(sw2[ll*nn]));
            ll += nprj[ia];
         }
      }

      ntmp = nn*nprjall;
      DSCAL_PWDFT(ntmp,scal,sw2,one);

      // DGEMM_PWDFT((char*) "N",(char*) "T",nshift,nn,nprjall,
      // 		  rmone,
      // 		  prjtmp,nshift,
      // 		  sw2,   nn,
      // 		  rone,
      // 		  Hpsi,nshift);
      gdevice_NT_dgemm(nshift,nn,nprjall,rmone,prjtmp,sw2,rone,Hpsi);
   }
   gdevice_hpsi_copy_gpu2host(nshift0,nn,Hpsi);
#else

   for (ii=0; ii<(myion->nion); ++ii)
   {
      ia = myion->katm[ii];
      if (nprj[ia]>0)
      {
         /* structure factor */
         mystrfac->strfac_pack(1,ii,exi);

         /* generate sw1's and projectors */
         for (l=0; l<nprj[ia]; ++l)
         {
            sd_function = !(l_projector[ia][l] & 1);
            prj = &(prjtmp[l*nshift]);
            vnlprj = &(vnl[ia][l*nshift0]);
            if (sd_function)
               mypneb->tcc_Mul( 1,vnlprj,exi,prj);
            else
               mypneb->tcc_iMul(1,vnlprj,exi,prj);
            //mypneb->cc_pack_indot(1,nn,psi,prj,&(sw1[l*nn]));
         }
         ntmp = nprj[ia];
         mypneb->cc_pack_inprjdot(1,nn,ntmp,psi,prjtmp,sw1);
         parall->Vector_SumAll(1,nn*nprj[ia],sw1);


         /* sw2 = Gijl*sw1 */
         Multiply_Gijl_sw1(nn,nprj[ia],nmax[ia],lmax[ia],
                           n_projector[ia],l_projector[ia],m_projector[ia],
                           Gijl[ia],sw1,sw2);

         /* do Kleinman-Bylander Multiplication */
         ntmp = nn*nprj[ia];
         DSCAL_PWDFT(ntmp,scal,sw2,one);

         ntmp = nprj[ia];

         // DGEMM_PWDFT((char*) "N",(char*) "T",nshift,nn,ntmp,
         //       rmone,
         //       prjtmp,nshift,
         //       sw2,   nn,
         //       rone,
         //       Hpsi,nshift);
	 gdevice_NT_dgemm(nshift,nn,ntmp,rmone,prjtmp,sw2,rone,Hpsi);

      } /*if nprj>0*/
   } /*ii*/
#endif

   delete [] sw2;
   delete [] sw1;
   delete [] prjtmp;
   delete [] exi;
}



/*******************************************
 *                                         *
 *     Pseudopotential::v_nonlocal_fion    *
 *                                         *
 *******************************************/
void Pseudopotential::v_nonlocal_fion(double *psi, double *Hpsi, const bool move, double *fion)
{
   nwpw_timing_function ftimer(6);
   bool done;
   int ii,ia,l,nshift0,sd_function,i,n;
   int jj,ll,jstart,jend,nprjall;

   double *exi;
   double *prjtmp,*sw1,*sw2,*prj,*vnlprj;
   double *Gx,*Gy,*Gz,*xtmp,*sum;
   Parallel *parall;
   double omega = mypneb->lattice->omega();
   //double scal = 1.0/lattice_omega();
   double scal = 1.0/omega;
   int one=1;
   int three=3;
   int ntmp,nshift,nn,ispin;

   double rone  = 1.0;
   double rmone = -1.0;

   nn     = mypneb->neq[0]+mypneb->neq[1];
   ispin  = mypneb->ispin;
   nshift0 = mypneb->npack(1);
   nshift = 2*mypneb->npack(1);
   exi    = new double[nshift];
   prjtmp = new double[nprj_max*nshift];
   sw1    = new double[nn*nprj_max];
   sw2    = new double[nn*nprj_max];

   //Copy psi to device
   gdevice_psi_copy_host2gpu(nshift0,nn,psi);
   gdevice_hpsi_copy_host2gpu(nshift0,nn,Hpsi);

   if (move)
   {
      xtmp = new double[nshift0];
      sum  = new double[3*nn*nprj_max];
      Gx = new double [mypneb->nfft3d];
      Gy = new double [mypneb->nfft3d];
      Gz = new double [mypneb->nfft3d];
      mypneb->tt_copy(mypneb->Gxyz(0),Gx);
      mypneb->tt_copy(mypneb->Gxyz(1),Gy);
      mypneb->tt_copy(mypneb->Gxyz(2),Gz);
      mypneb->t_pack(1,Gx);
      mypneb->t_pack(1,Gy);
      mypneb->t_pack(1,Gz);

      //Gx  = mypneb->Gpackxyz(1,0);
      //Gy  = mypneb->Gpackxyz(1,1);
      //Gz  = mypneb->Gpackxyz(1,2);
   }

   parall = mypneb->d3db::parall;

#if 1

   ii = 0;
   while (ii<(myion->nion))
   {
      ia      = myion->katm[ii];
      nprjall = 0;
      jstart  = ii;
      done = false;
      while (!done)
      {
         //generate projectors
         if (nprj[ia]>0)
         {
            mystrfac->strfac_pack(1,ii,exi);
            for (l=0; l<nprj[ia]; ++l)
            {
               sd_function = !(l_projector[ia][l] & 1);
               prj = &(prjtmp[(l+nprjall)*nshift]);
               vnlprj = &(vnl[ia][l*nshift0]);
               if (sd_function)
                  mypneb->tcc_Mul( 1,vnlprj,exi,prj);
               else
                  mypneb->tcc_iMul(1,vnlprj,exi,prj);

               if (move) 
               {
                  for (n=0; n<nn; ++n)
                  {
                     mypneb->cct_iconjgMul(1,prj,&psi[n*nshift],xtmp);
                     sum[3*n  +3*nn*(l+nprjall)] = mypneb->tt_pack_idot(1,Gx,xtmp);
                     sum[3*n+1+3*nn*(l+nprjall)] = mypneb->tt_pack_idot(1,Gy,xtmp);
                     sum[3*n+2+3*nn*(l+nprjall)] = mypneb->tt_pack_idot(1,Gz,xtmp);
                   }
               }
            }
            nprjall += nprj[ia];
         }
         ++ii;
         if (ii<(myion->nion))
         {
            ia = myion->katm[ii];
            done = ((nprjall+nprj[ia]) > nprj_max);
         }
         else
         {
            done = true;
         }
      }
      jend = ii;
      mypneb->cc_pack_inprjdot(1,nn,nprjall,psi,prjtmp,sw1);
      parall->Vector_SumAll(1,nn*nprjall,sw1);
      if (move) parall->Vector_SumAll(1,3*nn*nprjall,sum);

      /* sw2 = Gijl*sw1 */
      ll = 0;
      for (jj=jstart; jj<jend; ++jj)
      {
         ia = myion->katm[jj];
         if (nprj[ia]>0)
         {
            Multiply_Gijl_sw1(nn,nprj[ia],nmax[ia],lmax[ia],
                              n_projector[ia],l_projector[ia],m_projector[ia],
                              Gijl[ia],&(sw1[ll*nn]),&(sw2[ll*nn]));
            ll += nprj[ia];
         }
      }

      ntmp = nn*nprjall;
      DSCAL_PWDFT(ntmp,scal,sw2,one);

      // DGEMM_PWDFT((char*) "N",(char*) "T",nshift,nn,nprjall,
      // 		  rmone,
      // 		  prjtmp,nshift,
      // 		  sw2,   nn,
      // 		  rone,
      // 		  Hpsi,nshift);
      gdevice_NT_dgemm(nshift,nn,nprjall,rmone,prjtmp,sw2,rone,Hpsi);

      if (move)
         for (ll=0; ll<nprjall; ++ll)
            for (n=0; n<nn; ++n)
            {
               fion[3*ii]   +=  (3-ispin)*2.0*DDOT_PWDFT(nn, &sw2[ll*nn], one, &sum[  3*nn*ll], three);
               fion[3*ii+1] +=  (3-ispin)*2.0*DDOT_PWDFT(nn, &sw2[ll*nn], one, &sum[1+3*nn*ll], three);
               fion[3*ii+2] +=  (3-ispin)*2.0*DDOT_PWDFT(nn, &sw2[ll*nn], one, &sum[2+3*nn*ll], three);
            }

   }
   gdevice_hpsi_copy_gpu2host(nshift0,nn,Hpsi);



#else

   for (ii=0; ii<(myion->nion); ++ii)
   {
      ia = myion->katm[ii];
      if (nprj[ia]>0)
      {
         /* structure factor */
         mystrfac->strfac_pack(1,ii,exi);

         /* generate sw1's and projectors */
         for (l=0; l<nprj[ia]; ++ l)
         {
            sd_function = !(l_projector[ia][l] & 1);
            prj = &prjtmp[l*nshift];
            vnlprj = &vnl[ia][l*nshift0];
            if (sd_function)
               mypneb->tcc_Mul( 1,vnlprj,exi,prj);
            else
               mypneb->tcc_iMul(1,vnlprj,exi,prj);
            mypneb->cc_pack_indot(1,nn,psi,prj,&sw1[l*nn]);
         }
         parall->Vector_SumAll(1,nn*nprj[ia],sw1);

         /* sw2 = Gijl*sw1 */
         Multiply_Gijl_sw1(nn,nprj[ia],nmax[ia],lmax[ia],
                           n_projector[ia],l_projector[ia],m_projector[ia],
                           Gijl[ia],sw1,sw2);

         /* do Kleinman-Bylander Multiplication */
         ntmp = nn*nprj[ia];
         DSCAL_PWDFT(ntmp,scal,sw2,one);

         ntmp = nprj[ia];
         DGEMM_PWDFT((char*) "N",(char*) "T",nshift,nn,ntmp,
               rmone,
               prjtmp,nshift,
               sw2,   nn,
               rone,
               Hpsi,nshift);

         if (move)
         {
            for (l=0; l<nprj[ia]; ++ l)
            {
               prj = &prjtmp[l*nshift];
               for (n=0; n<nn; ++n)
               {
                  mypneb->cct_iconjgMul(1,prj,&psi[n*nshift],xtmp);
                  sum[3*n+1] = mypneb->tt_pack_idot(1,Gy,xtmp);
                  sum[3*n+2] = mypneb->tt_pack_idot(1,Gz,xtmp);
                  sum[3*n]   = mypneb->tt_pack_idot(1,Gx,xtmp);
                  //std::cout << ii << " " << n << " " << l << " " << sum[3*n] << " " << sum[3*n+1] << " " << sum[3*n+2] << "    SW2= " << sw2[n+l*nn] << std::endl;
                }
                parall->Vector_SumAll(1,3*nn,sum);

                fion[3*ii]   +=  (3-ispin)*2.0*DDOT_PWDFT(nn, &sw2[l*nn], one, sum,     three);
                fion[3*ii+1] +=  (3-ispin)*2.0*DDOT_PWDFT(nn, &sw2[l*nn], one, &sum[1], three);
                fion[3*ii+2] +=  (3-ispin)*2.0*DDOT_PWDFT(nn, &sw2[l*nn], one, &sum[2], three);
            }
         }

      } /*if nprj>0*/
   } /*ii*/
#endif

   if (move)
   {
      delete [] xtmp;
      delete [] sum;
      delete [] Gx;
      delete [] Gy;
      delete [] Gz;
   }
   delete [] sw2;
   delete [] sw1;
   delete [] prjtmp;
   delete [] exi;
}




/*******************************************
 *                                         *
 *     Pseudopotential::v_local            *
 *                                         *
 *******************************************/
void Pseudopotential::v_local(double *vout, const bool move, double *dng, double *fion)
{
   nwpw_timing_function ftimer(5);
   int ii,ia,nshift,npack0;
   double *exi,*vtmp,*xtmp,*Gx,*Gy,*Gz;

   npack0 = mypneb->npack(0);
   if (move)
   {
      xtmp = new double[npack0];
      Gx = new double [mypneb->nfft3d];
      Gy = new double [mypneb->nfft3d];
      Gz = new double [mypneb->nfft3d];
      mypneb->tt_copy(mypneb->Gxyz(0),Gx);
      mypneb->tt_copy(mypneb->Gxyz(1),Gy);
      mypneb->tt_copy(mypneb->Gxyz(2),Gz);
      mypneb->t_pack(0,Gx);
      mypneb->t_pack(0,Gy);
      mypneb->t_pack(0,Gz);

      //Gx  = mypneb->Gpackxyz(0,0);
      //Gy  = mypneb->Gpackxyz(0,1);
      //Gz  = mypneb->Gpackxyz(0,2);
   }

   mypneb->c_zero(0,vout);
   nshift = 2*npack0;
   exi    = new double[nshift];
   vtmp   = new double[nshift];
   for (ii=0; ii<(myion->nion); ++ii)
   {
      ia = myion->katm[ii];
      mystrfac->strfac_pack(0,ii,exi);
      //mypneb->tcc_MulSum2(0,vl[ia],exi,vout);
      mypneb->tcc_Mul(0,vl[ia],exi,vtmp);
      mypneb->cc_Sum2(0,vtmp,vout);

      if (move)
      {
         double xx =  mypneb->cc_pack_dot(0,dng,dng);
         double yy =  mypneb->cc_pack_dot(0,vtmp,vtmp);

         mypneb->cct_iconjgMulb(0,dng,vtmp,xtmp);
         double zz =  mypneb->tt_pack_dot(0,xtmp,xtmp);

         fion[3*ii]   = mypneb->tt_pack_dot(0,Gx,xtmp);
         fion[3*ii+1] = mypneb->tt_pack_dot(0,Gy,xtmp);
         fion[3*ii+2] = mypneb->tt_pack_dot(0,Gz,xtmp);
      }
   }
   delete [] exi;
   delete [] vtmp;
   if (move)
   {
       delete [] xtmp;
       delete [] Gx;
       delete [] Gy;
       delete [] Gz;
   }
}


/********************************************
 *                                          *
 * Pseudopotential::semicore_density_update *
 *                                          *
 ********************************************/

void Pseudopotential::semicore_density_update()
{
   int ii,ia;
   double omega = mypneb->lattice->omega();
   double scal2 = 1.0/omega;
   //double scal1 = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
   //double dv    = omega*scal1;
   double *exi = mypneb->c_pack_allocate(0);
   double *tmp = mypneb->r_alloc();

   mypneb->r_zero(semicore_density);
   for (ii=0; ii<(myion->nion); ++ii)
   {
      ia = myion->katm[ii];
      if (semicore[ia])
      {
         mystrfac->strfac_pack(0,ii,exi);
         mypneb->tcc_Mul(0,ncore_atom[ia],exi,tmp);

         /* Put put tmp into real space */
         mypneb->c_unpack(0,tmp);
         mypneb->cr_fft3d(tmp);

         /*  square it  */
         mypneb->r_sqr(tmp);
         mypneb->rr_Sum(tmp,semicore_density);
      }
   }
   mypneb->r_SMul(scal2*scal2,semicore_density);

   mypneb->r_dealloc(tmp);
   mypneb->c_pack_deallocate(exi);
}


