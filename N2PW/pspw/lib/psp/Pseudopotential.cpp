/* Pseudopotential.C - 
   Author - Eric Bylaska
*/
/*


#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cstdlib>
*/

#include	<string.h>
#include        <cmath>

#include	"lattice.hpp"
#include	"compressed_io.hpp"
#include	"blas.h"

#include	"Pseudopotential.hpp"


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

   dcopy_(&nnn,&rzero,&zero,sw2,&one);
   for (b=0; b<nprj; ++b)
      for (a=0; a<nprj; ++a)
         if ((l_prj[a]==l_prj[b]) && (m_prj[a]==m_prj[b]))
         {
            na = n_prj[a]-1;
            nb = n_prj[b]-1;
           daxpy_(&nna,
                   &G[nb + na*nmax + nmax2*l_prj[a]],
                   &sw1[a*nn],&one,
                   &sw2[b*nn],&one);
          }
}




static void psp_read(PGrid *mygrid,
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
                     int *semicore,
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
      *semicore = 1;
   else
      *semicore = 0;


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
  

/* Constructors */

/*******************************************
 *                                         *
 *     Pseudopotential::Pseudopotential    *
 *                                         *
 *******************************************/
Pseudopotential::Pseudopotential(Ion *myionin, Pneb *mypnebin, Strfac *mystrfacin)
{
   int ia,version,nfft[3];
   int *n_ptr,*l_ptr,*m_ptr,*b_ptr;
   double *rc_ptr,*G_ptr,*vnl_ptr,*ncore_ptr;
   double unita[9];
   char fname[80],aname[2];

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
   semicore = new int[npsp];

   n_projector = new int* [npsp];
   l_projector = new int* [npsp];
   m_projector = new int* [npsp];
   b_projector = new int* [npsp];

   zv          = new double[npsp];
   amass       = new double[npsp];
   rcore       = new double[npsp];
   rc          = new double* [npsp];
   vl          = new double* [npsp];
   for (ia=0; ia<npsp; ++ia) 
      vl[ia] = new double [mypneb->npack(0)];
   Gijl        = new double* [npsp];
   vnl         = new double* [npsp];
   ncore_atom  = new double* [npsp];
   
   comment  = new  char* [npsp];
   for (ia=0; ia<npsp; ++ia) comment[ia] = new char[80];

   for (ia=0; ia<npsp; ++ia)
   {
      strcpy(fname,myion->atom(ia));
      strcat(fname,".vpp");
      psp_read(mypneb,
               fname,
               comment[ia],&psp_type[ia],&version,nfft,unita,aname,
               &amass[ia],&zv[ia],&lmmax[ia],&lmax[ia],&locp[ia],&nmax[ia],
               &rc_ptr,&nprj[ia],&n_ptr,&l_ptr,&m_ptr,&b_ptr,&G_ptr,&semicore[ia],&rcore[ia],
               &ncore_ptr,vl[ia],&vnl_ptr);

      rc[ia]          = rc_ptr;
      n_projector[ia] = n_ptr;
      l_projector[ia] = l_ptr;
      m_projector[ia] = m_ptr;
      b_projector[ia] = b_ptr;
      Gijl[ia]        = G_ptr;
      vnl[ia]         = vnl_ptr;
      if (nprj[ia]>nprj_max) nprj_max = nprj[ia];
      if (semicore[ia])
         ncore_atom[ia]  = ncore_ptr;
   }

}
/*******************************************
 *                                         *
 *     Pseudopotential::v_nonlocal         *
 *                                         *
 *******************************************/
void Pseudopotential::v_nonlocal(double *psi, double *Hpsi)
{
   int ii,ia,l,nshift0,sd_function,i;
   double *exi;
   double *prjtmp,*sw1,*sw2,*prj,*vnlprj;
   Parallel *parall;
   double scal = 1.0/lattice_omega();
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
         dscal_(&ntmp,&scal,sw2,&one);

        ntmp = nprj[ia];

        dgemm_("N","T",&nshift,&nn,&ntmp,
               &rmone,
               prjtmp,&nshift,
               sw2,   &nn,
               &rone,
               Hpsi,&nshift);


      } /*if nprj>0*/
   } /*ii*/

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
void Pseudopotential::v_local(double *vout)
{
   int ii,ia,nshift;
   double *exi;

   mypneb->c_zero(0,vout);
   nshift = 2*mypneb->npack(0);
   exi    = new double[nshift];
   for (ii=0; ii<(myion->nion); ++ii)
   {
      ia = myion->katm[ii];
      mystrfac->strfac_pack(0,ii,exi);
      mypneb->tcc_MulSum2(0,vl[ia],exi,vout);

   }
   delete [] exi;
}
