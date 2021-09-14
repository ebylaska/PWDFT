/* Psp1d_Hamaan.cpp - 
   Author - Eric Bylaska
*/

#include	<iostream>
#include	<cstdio>
#include	<cstdlib>
#include	<cstring>
#include	<cmath>
//#include        "blas.h"

using namespace std;


#include	"Parallel.hpp"
#include	"util.hpp"
#include	"Psp1d_Hamann.hpp"

namespace pwdft {
using namespace pwdft;


#define FMT1    "%lf"
#define FMT2    " %lf %lf"
#define FMT10   "%10.3lf %10.3lf %10.3lf"

/*******************************************
 *                                         *
 *              util_matinvert             *
 *                                         *
 *******************************************/
/* Calculates matrix inverse based on Gauss-Jordan elimination 
  method with partial pivoting.
*/
static void util_matinvert(int n, int nmax, double *a)
{
   int irow;
   double big,tmp,pivinv;
   int *indx = new int[nmax];
   for (auto i=0; i<n; ++i)
   {
      big = 0.0;
      for (auto j=i; j<n; ++j)
         if (fabs(a[j+i*nmax]) >= big)
         {
            big = fabs(a[j+j*nmax]);
            irow = j;
         }
      if (big<=1.0e-9)
      {
         printf("Failed to invert matix\n");
         exit(99);
      }
      indx[i] = irow;

      if (irow!=i)
      {
         for (auto j=0; j<n; ++j)
         {
            tmp = a[irow+j*nmax];
            a[irow+j*nmax] = a[i+j*nmax];
            a[i+j*nmax]=tmp;
         }
      }

      pivinv = 1.0/a[i+i*nmax];
      a[i+i*nmax]= 1.0;

      for (auto j=0; j<n; ++j)
         a[i+j*nmax] *= pivinv;

      for (auto l=0; l<n; ++l)
         if (l!=i)
         {
            tmp = a[l+i*nmax];
            a[l+i*nmax] = 0.0;
            for (auto j=0; j<n; ++j)
               a[l+j*nmax] -= a[i+j*nmax]*tmp;
         }
   }

   delete [] indx;
}

/*******************************************
 *                                         *
 *              util_simpson               *
 *                                         *
 *******************************************/
static double util_simpson(int n, double *y, double h)
{
   double s = -y[0]-y[n-1];
   for (auto i=0; i<n; i+=2) s += 2.0*y[i];
   for (auto i=1; i<n; i+=2) s += 4.0*y[i];

   return s*h/3.0;
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
static void read_vpwp(Parallel *myparall, FILE *fp, int nrho, int lmax0, int lmax, int n_extra,
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
         for (l=0; l<=(lmax0+n_extra); ++l)
         {
            std::fscanf(fp,FMT1,&xx);
            if (l<=lmax) wp[i+l*nrho] = xx;
            if (l>lmax0) wp[i+(l+lmax-lmax0)*nrho] = xx;
         }
      }
   }
   myparall->Brdcst_Values(0,0,nrho,rho);
   myparall->Brdcst_Values(0,0,nrho*(lmax+1),vp);
   myparall->Brdcst_Values(0,0,nrho*(lmax+1+n_extra),wp);
}


/*******************************************
 *                                         *
 *               read_semicore             *
 *                                         *
 *******************************************/
static void read_semicore(Parallel *myparall, FILE *fp, int *isemicore, double *rcore, int nrho, double *semicore_r)
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
            semicore_r[i] = yy;
         }
      }
   }

   myparall->Brdcst_iValue(0,0,&isemicore0);
   myparall->Brdcst_Values(0,0,1,&rcore0);
   if (isemicore0>0)
      myparall->Brdcst_Values(0,0,2*nrho,semicore_r);
   *rcore = rcore0;
   *isemicore = isemicore0;
}

//static double dsum(int n, double *x, int incx)
//{
//   double stemp = 0.0;
//   for (int i=0; i<(n*incx); i+=incx) stemp += x[i];
//   return stemp;
//}

//static double simpson(int n, double *y, double h)
//{
//   int ne = n/2;
//   int no = ne+1;
//
//   double s = 3.0*dsum(no,y,2) + 4.0*dsum(ne,&y[1],2) - y[0] - y[n-1];
//   return (s*h/3.0);
//}


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
         coeff = fourpi*util_simpson(nrho,f,drho);
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
   int nn;

   if (myparall->is_master())
   {
      fp = std::fopen(psp_name,"r");
      std::fscanf(fp,"%s",atom);
      psp_type = convert_psp_type(atom);
      if (psp_type>0) std::fscanf(fp,"%s",atom);

      std::fscanf(fp,FMT1,&zv);
      std::fscanf(fp,FMT1,&amass);
      std::fscanf(fp,"%d",&lmax0);
      std::fscanf(fp,"%d",&lmax);
      std::fscanf(fp,"%d",&locp);
      std::fscanf(fp,FMT1,&rlocal);
      for (int l=0; l<=lmax0; ++l)
         std::fscanf(fp,FMT1,&rc[l]);

      n_extra = 0;
      for (int l=0; l<=lmax0; ++l) n_expansion[l] = 1;
      if (psp_type==2) 
      {
         std::fscanf(fp,"%d",&n_extra);
         for (int l=0; l<=lmax0; ++l)
         {
            std::fscanf(fp,"%d",&nn);
            n_expansion[l] = nn;
         }
      }

      std::fscanf(fp,"%d",&nrho);
      std::fscanf(fp,FMT1,&drho);
      std::fscanf(fp," %79[^\n]",comment);
      if (lmax>lmax0) lmax = lmax0;
      if (lmax<0)     lmax = lmax0;
      if (locp>lmax)  locp = lmax;
      if (locp<0)     locp = lmax;
   }
   myparall->Brdcst_cValues(0,0,2,atom);
   //myparall->Brdcst_iValue(0,0,&ihasae);
   myparall->Brdcst_iValue(0,0,&psp_type);
   myparall->Brdcst_Values(0,0,1,&zv);
   myparall->Brdcst_Values(0,0,1,&amass);
   myparall->Brdcst_iValue(0,0,&lmax0);
   myparall->Brdcst_iValue(0,0,&lmax);
   myparall->Brdcst_iValue(0,0,&locp);
   myparall->Brdcst_Values(0,0,(lmax0+1),rc);
   myparall->Brdcst_iValue(0,0,&nrho);
   myparall->Brdcst_Values(0,0,1,&drho);
   myparall->Brdcst_cValues(0,0,80,comment);

   myparall->Brdcst_iValue(0,0,&n_extra);
   myparall->Brdcst_iValues(0,0,(lmax0+1),n_expansion);
   nprj = 0;
   nmax = 1;
   for (int l=0; l<=(lmax); ++l)
      if (l!=locp)
      {
         nprj += n_expansion[l]*(2*l+1);
         if (n_expansion[l]>nmax) nmax = n_expansion[l];
      }

   rho = (double *) new double[nrho];
   vp  = (double *) new double[(lmax+1)*nrho];
   wp  = (double *) new double[(lmax+1+n_extra)*nrho];
   vnlnrm   = (double *) new double[nmax*nmax*(lmax+1)];
   rho_sc_r = (double *) new double[2*nrho];
   if (psp_type==9) 
   {
      up  = (double *) new double[(lmax+1)*nrho];
      r3_matrix = (double *) new double[(lmax+1)*(lmax+1)];
      read_vpwpup(myparall,fp,nrho,lmax0,lmax,rho,vp,wp,up);
      generate_r3_matrix(nrho,lmax,drho,rho,wp,up,r3_matrix);
   }
   else
       read_vpwp(myparall,fp,nrho,lmax0,lmax,n_extra,rho,vp,wp);
   read_semicore(myparall,fp,&isemicore,&rcore,nrho,rho_sc_r);
   semicore = (isemicore==1);

   /* Define non-local pseudopotential  */
   for (auto l=0; l<=lmax; ++l) 
      if (l!=locp)
         for (auto i=0; i<nrho; ++i)
            vp[i+l*nrho] = vp[i+l*nrho] - vp[i+locp*nrho];

   /* set up indx(n,l) --> to wp */
   int indx[5*4];
   int nb = lmax+1;
   for (auto l=0; l<=lmax; ++l)
   {
      indx[l*5] = l;
      for (auto n1=1; n1<n_expansion[l]; ++n1)
      {
         indx[n1+l*5] = nb;
         ++nb;
      }
   }

   version = 3;
   /* Normarization constants */
   double a;
   double *f = new double[nrho];
   for (auto l=0; l<=lmax; ++l)
   {
      if (l!=locp)
      {
         for (auto n2=0; n2<n_expansion[l]; ++n2)
         {
            for (auto i=0; i<nrho; ++i)
               f[i] = vp[i+l*nrho]*wp[i+indx[n2+l*5]*nrho]*wp[i+indx[n2+l*5]*nrho];

            a = util_simpson(nrho,f,drho);
            vnlnrm[n2+n2*nmax+l*nmax*nmax] = a;

            for (auto n1=n2+1; n1<n_expansion[l]; ++n1)
            {
               for (auto i=0; i<nrho; ++i)
                  f[i] = vp[i+l*nrho]*wp[i+indx[n1+l*5]*nrho]*wp[i+indx[n2+l*5]*nrho];
               a = util_simpson(nrho,f,drho);
               vnlnrm[n1+n2*nmax+l*nmax*nmax] = a;
               vnlnrm[n2+n1*nmax+l*nmax*nmax] = a;
  
            }
         }
         if (n_expansion[l]==1)
         {
            vnlnrm[l*nmax*nmax] = 1/a;
         }
         else if (n_expansion[l]==2)
         {
            double d = vnlnrm[l*nmax*nmax]*vnlnrm[1+1*nmax+l*nmax*nmax]
                     - vnlnrm[0+1*nmax+l*nmax*nmax]*vnlnrm[1+0*nmax+l*nmax*nmax];
            double q = vnlnrm[l*nmax*nmax];
            vnlnrm[l*nmax*nmax]          = vnlnrm[1+1*nmax+l*nmax*nmax]/d;
            vnlnrm[1+1*nmax+l*nmax*nmax] = q/d;
            vnlnrm[0+1*nmax+l*nmax*nmax] = -vnlnrm[0+1*nmax+l*nmax*nmax]/d;
            vnlnrm[1+0*nmax+l*nmax*nmax] = -vnlnrm[1+0*nmax+l*nmax*nmax]/d;
         }
         else
         {
            util_matinvert(n_expansion[l],nmax,&(vnlnrm[l*nmax*nmax]));
         }
      }
      else
      {
         for (int n2=0; n2<nmax; ++n2)
         for (int n1=n2; n1<nmax; ++n1)
         {
            vnlnrm[n1 + n2*nmax + l*nmax*nmax] = 0.0;
            vnlnrm[n2 + n1*nmax + l*nmax*nmax] = 0.0;
         }
      }
   }
   delete [] f;

   /* define n_prj, l_prj, m_prj, b_prj */
   if (nprj>0)
   {
      n_prj = new int[nprj];
      l_prj = new int[nprj];
      m_prj = new int[nprj];
      b_prj = new int[nprj];
      //std::cout << "define n_prj nprj=" << nprj << std::endl;

      int lcount = nprj;
      /* f projectors */
      if ((locp!=3) && (lmax>2))
         for (auto n=0; n<n_expansion[3]; ++n)
         {
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 3;
            m_prj[lcount] = -3;
            b_prj[lcount] = indx[n+3*5];

            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 3;
            m_prj[lcount] = -2;
            b_prj[lcount] = indx[n+3*5];
          
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 3;
            m_prj[lcount] = -1;
            b_prj[lcount] = indx[n+3*5];
          
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 3;
            m_prj[lcount] = 0;
            b_prj[lcount] = indx[n+3*5];
          
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 3;
            m_prj[lcount] = 1;
            b_prj[lcount] = indx[n+3*5];
          
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 3;
            m_prj[lcount] = 2;
            b_prj[lcount] = indx[n+3*5];
          
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 3;
            m_prj[lcount] = 3;
            b_prj[lcount] = indx[n+3*5];
         }

      /* d projectors */
      if ((locp!=2) && (lmax>1))
         for (auto n=0; n<n_expansion[2]; ++n)
         {
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 2;
            m_prj[lcount] = -2;
            b_prj[lcount] = indx[n+2*5];

            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 2;
            m_prj[lcount] = -1;
            b_prj[lcount] = indx[n+2*5];

            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 2;
            m_prj[lcount] = 0;
            b_prj[lcount] = indx[n+2*5];

            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 2;
            m_prj[lcount] = 1;
            b_prj[lcount] = indx[n+2*5];

            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 2;
            m_prj[lcount] = 2;
            b_prj[lcount] = indx[n+2*5];
         }

      /* p projectors */
      if ((locp!=1) && (lmax>0))
         for (auto n=0; n<n_expansion[1]; ++n)
         {
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 1;
            m_prj[lcount] = -1;
            b_prj[lcount] = indx[n+1*5];

            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 1;
            m_prj[lcount] = 0;
            b_prj[lcount] = indx[n+1*5];

            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 1;
            m_prj[lcount] = 1;
            b_prj[lcount] = indx[n+1*5];
         }

      /* s projectors */
      if (locp!=0)
         for (auto n=0; n<n_expansion[0]; ++n)
         {
            --lcount;
            n_prj[lcount] = n+1;
            l_prj[lcount] = 0;
            m_prj[lcount] = 0;
            b_prj[lcount] = indx[n+0*5];
         }
   }

}



/*******************************************
 *                                         *
 *     Psp1d_Hamann::vpp_generate_ray      *
 *                                         *
 *******************************************/
void Psp1d_Hamann::vpp_generate_ray(Parallel *myparall, int nray, double *G_ray, double *vl_ray, double *vnl_ray, double *rho_sc_k_ray)
{

   /* set up indx(n,l) --> to wp */
   int indx[5*4];
   int nb = lmax+1;
   for (auto l=0; l<=lmax; ++l)
   {
      indx[l*5] = l;
      for (auto n1=1; n1<n_expansion[l]; ++n1)
      {
         indx[n1+l*5] = nb;
         ++nb;
      }
   }

   double pi    = 4.00*atan(1.0);
   double twopi = 2.0*pi;
   double forpi = 4.0*pi;

   double P0 = sqrt(forpi);
   double P1 = sqrt(3.0*forpi);
   double P2 = sqrt(15.0*forpi);
   double P3 = sqrt(105.0*forpi);

   double zero = 0.0;
   int    izero = 0;
   int    ione  = 1;
   int    nray2 = 2*nray;
   int    lmaxnray = (lmax+1+n_extra)*nray;

   double q;
   double *cs = new double[nrho];
   double *sn = new double[nrho];
   double *f  = new double[nrho];
   double a,xx;

   //std::cout << "gen vl_ray version=" << version<< " " <<  std::endl;
   //dcopy_(&nray,    &zero,&izero,vl_ray,&ione);
   //dcopy_(&lmaxnray,&zero,&izero,vnl_ray,&ione);
   //dcopy_(&nray2,   &zero,&izero,rho_sc_k_ray,&ione);
  
   memset(vl_ray,0,nray*sizeof(double));
   memset(vnl_ray,0,lmaxnray*sizeof(double));
   memset(rho_sc_k_ray,0,nray2*sizeof(double));

   for (auto k1=(1+myparall->taskid()); k1<nray; k1+=myparall->np())
   {
      q=G_ray[k1];
      for (auto i=0; i<nrho; ++i)
      {
         cs[i]=cos(q*rho[i]);
         sn[i]=sin(q*rho[i]);
      }

      /* h projectors */
      /* h projectors */
      /* f projectors */
      if ((locp!=3) && (lmax>2))
      {
         for (auto n=0; n<n_expansion[3]; ++n)
         {
            f[0] = 0.0;
            for (auto i=1; i<nrho; ++i)
            {
                xx = q*rho[i];
                a  = sn[i]/xx;
                a  = 15.0*(a-cs[i])/(xx*xx) - 6*a + cs[i];
                f[i]=a*wp[i+indx[n+3*5]*nrho]*vp[i+3*nrho];
            }
            vnl_ray[k1+indx[n+3*5]*nray]=P3*util_simpson(nrho,f,drho)/q;
         } 
      }

      /* d projectors */
      if ((locp!=2) && (lmax>1))
      {
         for (auto n=0; n<n_expansion[2]; ++n)
         {
            f[0] = 0.0;
            for (auto i=1; i<nrho; ++i)
            {
                a=3.0*(sn[i]/(q*rho[i])-cs[i])/(q*rho[i])-sn[i];
                f[i]=a*wp[i+indx[n+2*5]*nrho]*vp[i+2*nrho];
            }
            vnl_ray[k1+indx[n+2*5]*nray]=P2*util_simpson(nrho,f,drho)/q;
         } 
      }

      /* p projectors */
      if ((locp!=1) && (lmax>0))
      {
         for (auto n=0; n<n_expansion[1]; ++n)
         {
            f[0] = 0.0;
            for (auto i=1; i<nrho; ++i)
            {
                a = (sn[i]/(q*rho[i])-cs[i]);
                f[i]=a*wp[i+indx[n+1*5]*nrho]*vp[i+1*nrho];
            }
            vnl_ray[k1+indx[n+1*5]*nray]=P1*util_simpson(nrho,f,drho)/q;
         } 
      }

      /* s projectors */
      if (locp!=0)
      {
         for (auto n=0; n<n_expansion[0]; ++n)
         {
            for (auto i=0; i<nrho; ++i)
                f[i]=sn[i]*wp[i+indx[n+0*5]*nrho]*vp[i+0*nrho];
            vnl_ray[k1+indx[n+0*5]*nray]=P0*util_simpson(nrho,f,drho)/q;
         } 
      }

      /* local */
      if (version==3)
      {
         for (auto i=0; i<nrho; ++i)
            f[i]=rho[i]*vp[i+locp*nrho]*sn[i];
         vl_ray[k1]=util_simpson(nrho,f,drho)*forpi/q - zv*forpi*cs[nrho-1]/(q*q);
      }
      else if (version==4)
      {
         for (auto i=0; i<nrho; ++i)
            f[i]=(rho[i]*vp[i+locp*nrho]+zv*erf(rho[i]/rlocal))*sn[i];
         vl_ray[k1]=util_simpson(nrho,f,drho)*forpi/q;
      }

      /* semicore density */
      if (semicore)
      {
         for (auto i=0; i<nrho; ++i)
            f[i]=rho[i]*sqrt(rho_sc_r[i])*sn[i];
         rho_sc_k_ray[k1]=util_simpson(nrho,f,drho)*forpi/q;

         for (auto i=0; i<nrho; ++i)
            f[i]=(sn[i]/(q*rho[i])-cs[i])*rho_sc_r[i+nrho]*rho[i];
         rho_sc_k_ray[k1+nray]=util_simpson(nrho,f,drho)*forpi/q;
      }
   }
   myparall->Vector_SumAll(0,2*nray,rho_sc_k_ray);
   myparall->Vector_SumAll(0,nray,vl_ray);
   myparall->Vector_SumAll(0,lmaxnray,vnl_ray);

   /* G==0 local */
   if (version==3)
   {
      for (auto i=0; i<nrho; ++i)
      {
        f[i]=vp[i+locp*nrho]*rho[i]*rho[i];
      }
      vl_ray[0]=forpi*util_simpson(nrho,f,drho)+twopi*zv*rho[nrho-1]*rho[nrho-1];
   }
   else if (version==4) 
   {
      for (auto i=0; i<nrho; ++i)
        f[i]= (vp[i+locp*nrho]*rho[i]+zv*erf(rho[i]/rlocal))*rho[i];
      vl_ray[0]=forpi*util_simpson(nrho,f,drho);
   }

   /* G==0 semicore */
   if (semicore)
   {
      for (auto i=0; i<nrho; ++i)
         f[i] = sqrt(rho_sc_r[i])*rho[i]*rho[i];
      rho_sc_k_ray[0]      = forpi*util_simpson(nrho,f,drho);
      rho_sc_k_ray[0+nray] = 0.0;
   }

   /* G==0 vnl */
   for (auto l=0; l<=lmax; ++l)
      for (auto n=0; n<n_expansion[l]; ++n)
         vnl_ray[0+indx[n+l*5]*nray] = 0.0;   

   /* only j0 is non-zero at zero */
   if (locp!=0)
      for (auto n=0; n<n_expansion[0]; ++n) 
      {
         for (auto i=0; i<nrho; ++i)
            f[i]=rho[i]*wp[i+indx[n+0*5]*nrho]*vp[i+0*nrho];
         vnl_ray[0+indx[n+0*5]*nray]=P0*util_simpson(nrho,f,drho);
      }

   delete [] f;
   delete [] sn;
   delete [] cs;
}



/*******************************************
 *                                         *
 *   Psp1d_Hamann::vpp_generate_spline     *
 *                                         *
 *******************************************/
void Psp1d_Hamann::vpp_generate_spline(PGrid *mygrid, int nray, double *G_ray, double *vl_ray, double *vnl_ray, double *rho_sc_k_ray,
                                       double *vl, double *vnl, double *rho_sc_k)
{

   /* set up indx(n,l) --> to wp */
   int indx[5*4];
   int nb = lmax+1;
   for (auto l=0; l<=lmax; ++l)
   {
      indx[l*5] = l;
      for (auto n1=1; n1<n_expansion[l]; ++n1)
      {
         indx[n1+l*5] = nb;
         ++nb;
      }
   }

   double pi    = 4.00*atan(1.0);

   /* allocate spline grids */
   double *vl_splineray  = new double [nray];
   double *vnl_splineray = new double [(lmax+1+n_extra)*nray];
   double *rho_sc_k_splineray = new double [2*nray];
   double *tmp_splineray  = new double [nray];


   /* setup cubic bsplines */
   double dG = G_ray[2] - G_ray[1];

   /* five point formula */
   double yp1 = ( -50.00*vl_ray[1]
                 + 96.00*vl_ray[2]
                 - 72.00*vl_ray[3]
                 + 32.00*vl_ray[4]
                 -  6.00*vl_ray[5])/(24.00*dG);
   util_spline(&(G_ray[1]),&(vl_ray[1]),nray-1,yp1,0.00,&(vl_splineray[1]),tmp_splineray);


   for (auto l=0; l<=lmax; ++l)
      if (l!=locp)
         for (auto n=0; n<n_expansion[l]; ++n)
            util_spline(G_ray,&(vnl_ray[indx[n+5*l]*nray]),nray,0.00,0.00,&(vnl_splineray[indx[n+5*l]*nray]),tmp_splineray);
   
   if (semicore)
   {
      util_spline(G_ray,rho_sc_k_ray,nray,0.00,0.00,rho_sc_k_splineray,tmp_splineray);
      util_spline(G_ray,&(rho_sc_k_ray[nray]),nray,0.00,0.00,&(rho_sc_k_splineray[nray]),tmp_splineray);
   }


   double q, qx, qy, qz, xx;
   double *gx, *gy, *gz;
   int npack0 = mygrid->npack(0);
   int npack1 = mygrid->npack(1);
   int nx,lcount;
   mygrid->t_pack_nzero(0,1,vl);
   mygrid->t_pack_nzero(1,nprj,vnl);
   if (semicore) mygrid->t_pack_nzero(0,4,rho_sc_k);

   /* generate vl and rho_sc_k */
   gx = mygrid->Gpackxyz(0,0);
   gy = mygrid->Gpackxyz(0,1);
   gz = mygrid->Gpackxyz(0,2);
   for (auto k=0; k<npack0; ++k)
   {
      qx = gx[k]; qy = gy[k]; qz = gz[k];
      q = sqrt(qx*qx + qy*qy + qz*qz);
      nx = (int) floor(q/dG);

      if (q>1.0e-9)
      {
         qx /= q; qy /= q; qz /= q;
         vl[k] = util_splint(&(G_ray[1]),&(vl_ray[1]),&(vl_splineray[1]),nray-1,nx,q);
         if (semicore)
         {
            rho_sc_k[k] = util_splint(G_ray,rho_sc_k_ray,rho_sc_k_splineray,nray,nx,q);
            xx          = util_splint(G_ray,&(rho_sc_k_ray[nray]),&(rho_sc_k_splineray[nray]),nray,nx,q);
            rho_sc_k[k+npack0]   = xx*qx;
            rho_sc_k[k+2*npack0] = xx*qy;
            rho_sc_k[k+3*npack0] = xx*qz;
         }
      }
      else
      {

         vl[k] = vl_ray[0];
         if (semicore)
         {
            rho_sc_k[k] = rho_sc_k_ray[0];
            rho_sc_k[k+npack0]   = 0.0;
            rho_sc_k[k+2*npack0] = 0.0;
            rho_sc_k[k+3*npack0] = 0.0;
         }
      }
   }

   /* generate vnl */
   gx = mygrid->Gpackxyz(1,0);
   gy = mygrid->Gpackxyz(1,1);
   gz = mygrid->Gpackxyz(1,2);



   for (auto k=0; k<npack1; ++k)
   {
      qx = gx[k]; qy = gy[k]; qz = gz[k];
      q = sqrt(qx*qx + qy*qy + qz*qz);
      nx = (int) floor(q/dG);

      if (q>1.0e-9)
      {
         qx /= q; qy /= q; qz /= q;
         lcount = nprj;

         /* f projectors */

         if ((locp!=3) && (lmax>2))
            for (auto n=0; n<n_expansion[3]; ++n)
            {
               xx = util_splint(G_ray,&(vnl_ray[indx[n+3*5]*nray]),&(vnl_splineray[indx[n+3*5]*nray]),nray,nx,q);
               --lcount;
               vnl[k+lcount*npack1]=  xx * qy*(3.00*(1.00-qz*qz)-4.00*qy*qy)/sqrt(24.00);
               --lcount;
               vnl[k+lcount*npack1]=  xx * qx*qy*qz;
               --lcount;
               vnl[k+lcount*npack1]=  xx * qy*(5.00*qz*qz-1.00)/sqrt(40.00);
               --lcount;
               vnl[k+lcount*npack1]=  xx * qz*(5.00*qz*qz-3.00)/sqrt(60.00);
               --lcount;
               vnl[k+lcount*npack1]=  xx * qx*(5.00*qz*qz-1.00)/sqrt(40.00);
               --lcount;
               vnl[k+lcount*npack1]=  xx * qz*(qx*qx - qy*qy)/2.00;
               --lcount;
               vnl[k+lcount*npack1]=  xx * qx*(4.00*qx*qx-3.00*(1.00-qz*qz))/sqrt(24.00);

            }

         /* d projectors */
         if ((locp!=2) && (lmax>1))
            for (auto n=0; n<n_expansion[2]; ++n)
            {
               xx = util_splint(G_ray,&(vnl_ray[indx[n+2*5]*nray]),&(vnl_splineray[indx[n+2*5]*nray]),nray,nx,q);
               --lcount;
               vnl[k+lcount*npack1]=  xx * qx*qy;
               --lcount;
               vnl[k+lcount*npack1]=  xx * qy*qz;
               --lcount;
               vnl[k+lcount*npack1]=  xx * (3.00*qz*qz-1.00)/(2.00*sqrt(3.00));
               --lcount;
               vnl[k+lcount*npack1]=  xx * qz*qx;
               --lcount;
               vnl[k+lcount*npack1]=  xx * (qx*qx-qy*qy)/(2.00);

            }

         /* p projectors */
         if ((locp!=1) && (lmax>0))
            for (auto n=0; n<n_expansion[1]; ++n)
            {
               xx = util_splint(G_ray,&(vnl_ray[indx[n+1*5]*nray]),&(vnl_splineray[indx[n+1*5]*nray]),nray,nx,q);
               --lcount;
               vnl[k+lcount*npack1]=  xx * qy;
               --lcount;
               vnl[k+lcount*npack1]=  xx * qz;
               --lcount;
               vnl[k+lcount*npack1]=  xx * qx;
            }

         /* s projectors */
         if (locp!=0)
            for (auto n=0; n<n_expansion[0]; ++n)
            {
               xx = util_splint(G_ray,&(vnl_ray[indx[n+0*5]*nray]),&(vnl_splineray[indx[n+0*5]*nray]),nray,nx,q);
               --lcount;
               vnl[k+lcount*npack1]=  xx;
            }
      }
      else
      {
         for (auto l=0; l<nprj; ++l)
            vnl[k+l*npack1] = 0.0;

         /* only j0 is non-zero at zero */
         if (locp!=0)
            for (auto n=0; n<n_expansion[0]; ++n)
               vnl[k+indx[n+0*5]*npack1] =  vnl_ray[0+indx[n+0*5]*nray];
      }
   }


   /*  deallocate spineray formatted grids */
   delete [] tmp_splineray;
   delete [] rho_sc_k_splineray;
   delete [] vnl_splineray;
   delete [] vl_splineray;

}

}


