/* Psp1d_Hamaan.cpp - 
   Author - Eric Bylaska
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

#include	"Parallel.hpp"
#include	"Psp1d_Hamann.hpp"


#define FMT1    "%lf"
#define FMT2    " %lf %lf"
#define FMT10   "%10.3lf %10.3lf %10.3lf"


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
 *               read_vpwpup               *
 *                                         *
 *******************************************/
static void read_vpwpup(Parallel *myparall, FILE *fp, int nrho, int lmax0, int lmax,
                        double *rho, double *vp, double *wp, double *up)
{
   int i,l;
   double xx;

   if (myparall->is_master())
   {
      for (i=0; i<nrho; ++i)
      {
         std::fscanf(fp,FMT1,&rho[i]);
         for (l=0; l<=lmax0; ++l)
         {
            std::fscanf(fp,FMT1,&xx);
            if (l<=lmax) vp[i+l*nrho] = xx;
         }
      }
      for (i=0; i<nrho; ++i)
      {
         std::fscanf(fp,FMT1,&rho[i]);
         for (l=0; l<=lmax0; ++l)
         {
            std::fscanf(fp,FMT1,&xx);
            if (l<=lmax) wp[i+l*nrho] = xx;
         }
         for (l=0; l<=lmax0; ++l)
         {
            std::fscanf(fp,FMT1,&xx);
            if (l<=lmax) up[i+l*nrho] = xx;
         }
      }
   }
   myparall->Brdcst_Values(0,0,nrho,rho);
   myparall->Brdcst_Values(0,0,nrho*(lmax+1),vp);
   myparall->Brdcst_Values(0,0,nrho*(lmax+1),wp);
   myparall->Brdcst_Values(0,0,nrho*(lmax+1),up);
}

/*******************************************
 *                                         *
 *               read_vpwp                 *
 *                                         *
 *******************************************/
static void read_vpwp(Parallel *myparall, FILE *fp, int nrho, int lmax0, int lmax,
                      double *rho, double *vp, double *wp)
{
   int i,l;
   double xx;

   if (myparall->is_master())
   {
      for (i=0; i<nrho; ++i)
      {
         std::fscanf(fp,FMT1,&rho[i]);
         for (l=0; l<=lmax0; ++l)
         {
            std::fscanf(fp,FMT1,&xx);
            if (l<=lmax) vp[i+l*nrho] = xx;
         }
      }
      for (i=0; i<nrho; ++i)
      {
         std::fscanf(fp,FMT1,&rho[i]);
         for (l=0; l<=lmax0; ++l)
         {
            std::fscanf(fp,FMT1,&xx);
            if (l<=lmax) wp[i+l*nrho] = xx;
         }
      }
   }
   myparall->Brdcst_Values(0,0,nrho,rho);
   myparall->Brdcst_Values(0,0,nrho*(lmax+1),vp);
   myparall->Brdcst_Values(0,0,nrho*(lmax+1),wp);
}


/*******************************************
 *                                         *
 *               read_semicore             *
 *                                         *
 *******************************************/
static void read_semicore(Parallel *myparall, FILE *fp, int *isemicore, double *rcore, int nrho, double *semicore)
{
   int    i,isemicore0;
   double xx,yy,rcore0;

   isemicore0 = 0;
   rcore0 = 0.0;
   if (myparall->is_master())
   {
      if (std::fscanf(fp,FMT1,&xx)!=EOF)
      {
         rcore0 = xx;
         isemicore0 = 1;
         for (i=0; i<(2*nrho); ++i)
         {
            std::fscanf(fp,FMT2,&xx,&yy);
            semicore[i] = yy;
         }
      }
   }
   myparall->Brdcst_iValue(0,0,&isemicore0);
   myparall->Brdcst_Values(0,0,1,&rcore0);
   if (isemicore0>0)
      myparall->Brdcst_Values(0,0,2*nrho,semicore);
   *rcore = rcore0;
   *isemicore = isemicore0;
}

static double dsum(int n, double *x, int incx)
{
   double stemp = 0.0;
   for (int i=0; i<(n*incx); i+=incx) stemp += x[i];
   return stemp;
}

static double simpson(int n, double *y, double h)
{
   int ne = n/2;
   int no = ne+1;

   double s = 3.0*dsum(no,y,2) + 4.0*dsum(ne,&y[1],2) - y[0] - y[n-1];
   return (s*h/3.0);
}


/*******************************************
 *                                         *
 *             generate_r3_matrix          *
 *                                         *
 *******************************************/
/*
   Computes the matrix elements for EFG tensor calculations

   r3_matrix(li,lj) = <uli|1/r3|ulj> - <wli|1/r3|wlj>
*/
static void generate_r3_matrix(int nrho, int lmax, double drho, 
                               double *rho, double *wp, double *up, 
                               double *r3_matrix)
{
   int i,li,lj;
   int lmax2 = (lmax+1)*(lmax+1);
   double coeff;
   double fourpi = 16.0*atan(1.0);
   double *f = (double *) new double[nrho];


   for (i=0; i<lmax2; ++i) r3_matrix[i] = 0.0;

   for (lj=0; lj<=lmax; ++lj)
   for (li=lj; li<=lmax; ++li)
      if ((li+lj)>0)
      {
         for (i=0; i<nrho; ++i)
            f[i] = ( up[i+li*nrho]*up[i+lj*nrho]
                   - wp[i+li*nrho]*wp[i+lj*nrho])
                  /(rho[i]*rho[i]*rho[i]);
         coeff = fourpi*simpson(nrho,f,drho);
         r3_matrix[li+lj*lmax] = coeff;
         if (li!=lj) r3_matrix[lj+li*lmax] = coeff;
      }

   delete [] f;
}




/* Constructors */

/*******************************************
 *                                         *
 *     Psp1d_Hamann::Psp1d_Hamann          *
 *                                         *
 *******************************************/
Psp1d_Hamann::Psp1d_Hamann(Parallel *myparall, const char *psp_name)
{
   double xx;
   FILE *fp;

   if (myparall->is_master())
   {
      fp = std::fopen(psp_name,"r");
      std::fscanf(fp,"%s",atom);
      ihasae = convert_psp_type(atom);
      if (ihasae>0) std::fscanf(fp,"%s",atom);

      std::fscanf(fp,FMT1,&zv);
      std::fscanf(fp,FMT1,&amass);
      std::fscanf(fp,"%d",&lmax0);
      std::fscanf(fp,"%d",&lmax);
      std::fscanf(fp,"%d",&locp);
      std::fscanf(fp,FMT1,&rlocal);
      for (int l=0; l<=lmax0; ++l)
         std::fscanf(fp,FMT1,&rc[l]);
      std::fscanf(fp,"%d",&nrho);
      std::fscanf(fp,FMT1,&drho);
      std::fscanf(fp," %79[^\n]",comment);
      if (lmax>lmax0) lmax = lmax0;
      if (lmax<0)     lmax = lmax0;
      if (locp>lmax)  locp = lmax;
      if (locp<0)     locp = lmax;
   }
   myparall->Brdcst_cValues(0,0,2,atom);
   myparall->Brdcst_iValue(0,0,&ihasae);
   myparall->Brdcst_Values(0,0,1,&zv);
   myparall->Brdcst_Values(0,0,1,&amass);
   myparall->Brdcst_iValue(0,0,&lmax0);
   myparall->Brdcst_iValue(0,0,&lmax);
   myparall->Brdcst_iValue(0,0,&locp);
   myparall->Brdcst_Values(0,0,(lmax0+1),rc);
   myparall->Brdcst_iValue(0,0,&nrho);
   myparall->Brdcst_Values(0,0,1,&drho);
   myparall->Brdcst_cValues(0,0,80,comment);

   rho = (double *) new double[nrho];
   vp  = (double *) new double[(lmax+1)*nrho];
   wp  = (double *) new double[(lmax+1)*nrho];
   semicore = (double *) new double[2*nrho];
   if (ihasae>0) 
   {
      up  = (double *) new double[(lmax+1)*nrho];
      r3_matrix = (double *) new double[(lmax+1)*(lmax+1)];
   }
   if (ihasae>0)
   {
       read_vpwpup(myparall,fp,nrho,lmax0,lmax,rho,vp,wp,up);
       generate_r3_matrix(nrho,lmax,drho,rho,wp,up,r3_matrix);
   }
   else
       read_vpwp(myparall,fp,nrho,lmax0,lmax,rho,vp,wp);
   read_semicore(myparall,fp,&isemicore,&rcore,nrho,semicore);

}

