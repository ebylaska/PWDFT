/* Pneb.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/

#include        <iostream>
#include        <cstdio>
#include        <cstring> //memset()
#include        <cmath>
#include        <cstdlib>
#include        <sstream>
#include        <stdexcept> // runtime_error()

#include	"Pneb.hpp"

#include "blas.h"

namespace pwdft {



/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
//Pneb::Pneb(Parallel *inparall, Lattice *inlattice, Control2& control, int ispin, int *ne)
//  : PGrid(inparall, inlattice, control),
//    d1db(inparall,control.mapping1d(),ispin,ne)

Pneb::Pneb(Parallel *inparall, Lattice *inlattice, Control2& control, int ispin, int *ne)
  : d1db(inparall,control.mapping1d(),ispin,ne),
    PGrid(inparall, inlattice, control)
{

    int np_i = d1db::parall->np_i();
    int np_j = d1db::parall->np_j();
    parallelized = (np_j>1);

    s22 = new (std::nothrow) double[7*ne[0]*ne[0]]();
    s21 = &s22[1*ne[0]*ne[0]];
    s12 = &s22[2*ne[0]*ne[0]];
    s11 = &s22[3*ne[0]*ne[0]];
    sa1 = &s22[4*ne[0]*ne[0]];
    sa0 = &s22[5*ne[0]*ne[0]];
    st1 = &s22[6*ne[0]*ne[0]];

    if (parallelized)
    {
       mcq[0] = mcq[1] = 0;
       ncq[0] = ncq[1] = 0;
       ncqmax = 0;
       for (auto ms=0; ms<ispin; ++ms)
       {
          ma[ms]  = new (std::nothrow) int[np_i];
          ma1[ms] = new (std::nothrow) int[np_i];
          ma2[ms] = new (std::nothrow) int[np_i];
          mc[ms]  = new (std::nothrow) int[np_i];
          na[ms]  = new (std::nothrow) int[np_j];
          nc[ms]  = new (std::nothrow) int[np_j];
       }
    }

    g_rnd_algorihm = control.initial_psi_random_algorithm();
}



/*************************************
 *                                   *
 *      Pneb::g_generate1_random     *
 *                                   *
 *************************************/
void Pneb::g_generate1_random(double *psi)
{
   int ms,n,indx,i,pj,qj,taskid_j;
   int filling[4],nfft[3];
   double zvalue[2];
   //double *tmp2 = new (std::nothrow) double[n2ft3d]();

   double tmp2[n2ft3d];

   nfft[0] = nx;
   nfft[1] = ny;
   nfft[2] = nz;

   taskid_j = d1db::parall->taskid_j();
   for (ms=0; ms<ispin; ++ms)
   for (n=0; n<ne[ms]; ++n)
   {
      util_getfilling(n,nfft,filling,zvalue);

      qj = msntoindex(ms,n);
      pj = msntop(ms,n);

      if (pj==taskid_j)
      {
         r_zero(tmp2);
         d3db::c_setpw(filling, zvalue, tmp2);
         d3db::c_addrandom(tmp2);

         PGrid::c_pack(1,tmp2);
         indx = 2*PGrid::npack(1)*qj;
         PGrid::cc_pack_copy(1,tmp2,psi+indx);
         PGrid::c_pack_noimagzero(1,psi+indx);
      }
   }

   //delete [] tmp2;
}

void Pneb::g_generate2_random(double *psi)
{
   double tmp2[n2ft3d];

   int taskid_j = d1db::parall->taskid_j();
   for (auto ms=0; ms<ispin; ++ms)
   for (auto n=0; n<ne[ms]; ++n)
   {
      int qj = msntoindex(ms,n);
      int pj = msntop(ms,n);
      if (pj==taskid_j)
      {
         d3db::r_setrandom(tmp2);
         d3db::r_zero_ends(tmp2);
         d3db::rc_fft3d(tmp2);

         PGrid::c_pack(1,tmp2);
         int indx = 2*PGrid::npack(1)*qj;
         PGrid::cc_pack_copy(1,tmp2,psi+indx);
      }
   }

}


/*************************************
 *                                   *
 *      Pneb::g_generate_random      *
 *                                   *
 *************************************/
void Pneb::g_generate_random(double *psi)
{
   int taskid = d1db::parall->taskid();
   util_random(taskid+91);

   if (g_rnd_algorihm==1)
      this->g_generate1_random(psi);
   else
      this->g_generate2_random(psi);
}



/*************************************
 *                                   *
 *           Pneb::g_read            *
 *                                   *
 *************************************/
void Pneb::g_read(const int iunit, double *psi)
{
   int ms,n,indx,i,pj,qj,taskid_j;
   double *tmp2 = new (std::nothrow) double[n2ft3d]();

   taskid_j = d1db::parall->taskid_j();

   for (ms=0; ms<ispin; ++ms)
   for (n=0; n<ne[ms]; ++n)
   {
      qj = msntoindex(ms,n);
      pj = msntop(ms,n);
      c_read(iunit,tmp2,pj);
      if (pj==taskid_j)
      {
         indx = 2*PGrid::npack(1)*qj;
         PGrid::c_pack(1,tmp2);
         PGrid::cc_pack_copy(1,tmp2,psi+indx);
      }
   }

   delete [] tmp2;
}

void Pneb::g_write(const int iunit, double *psi)
{
   int ms,n,indx,i,pj,qj,taskid_j;
   double *tmp2 = new (std::nothrow) double[n2ft3d]();

   taskid_j = d1db::parall->taskid_j();
   for (ms=0; ms<ispin; ++ms)
   for (n=0; n<ne[ms]; ++n)
   {
      qj = msntoindex(ms,n);
      pj = msntop(ms,n);
      if (pj==taskid_j)
      {
         indx = 2*PGrid::npack(1)*qj;
         PGrid::cc_pack_copy(1,psi+indx,tmp2);
         PGrid::c_unpack(1,tmp2);
      }
      c_write(iunit,tmp2,pj);
   }

   delete [] tmp2;
}


double Pneb::gg_traceall(double *psi1, double *psi2)
{
   int n,indx;
   double sum=0.0;

   indx = 0;
   for (n=0; n<(neq[0]+neq[1]); ++n)
   {
      sum += PGrid::cc_pack_idot(1,psi1+indx,psi2+indx);
      indx += 2*PGrid::npack(1);
   }
   if (ispin==1) sum *= 2.0;

   return d3db::parall->SumAll(0,sum);
}


void Pneb::gg_copy(double *psi1, double *psi2)
{
   int one=1;
   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
   std::memcpy(psi2,psi1,nsize*sizeof(double));
   //DCOPY_PWDFT(nsize, psi1, one, psi2, one);
}
void Pneb::gg_SMul(double alpha,double *psi1, double *psi2)
{
   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
   for (int i=0; i<nsize; ++i)
      psi2[i] = alpha*psi1[i];
}
//void Pneb::g_SMul1(double alpha,double *psi1)
//{
//   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
//   for (int i=0; i<nsize; ++i)
//      psi1[i] *= alpha;
//}

void Pneb::gg_daxpy(double alpha, double *psi1, double *psi2)
{
   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
   for (int i=0; i<nsize; ++i)
      psi2[i] += alpha*psi1[i];
}

void Pneb::g_Scale(double alpha,double *psi1)
{
   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
   for (int i=0; i<nsize; ++i)
      psi1[i] *= alpha;
}

void Pneb::gg_Sum2(double *psi1, double *psi2)
{
   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
   for (int i=0; i<nsize; ++i)
      psi2[i] += psi1[i];
}

/*************************************
 *                                   *
 *         Pneb::gg_Minus2           *
 *                                   *
 *************************************/
void Pneb::gg_Minus2(double *psi1, double *psi2)
{
   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
   for (int i=0; i<nsize; ++i)
      psi2[i] -= psi1[i];
}

/*************************************
 *                                   *
 *         Pneb::ggg_Minus           *
 *                                   *
 *************************************/
void Pneb::ggg_Minus(double *psi1, double *psi2, double *psi3)
{
   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
   for (int i=0; i<nsize; ++i)
      psi3[i] = psi1[i] - psi2[i];
}

/*************************************
 *                                   *
 *         Pneb::g_zero              *
 *                                   *
 *************************************/
void Pneb::g_zero(double *psi2)
{
   int one=1;
   int zero=0;
   int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
   double rzero=0.0;

   //dcopy_(&nsize,&rzero,&zero,psi2,&one);
   std::memset(psi2, 0, nsize*sizeof(double));
}

/*************************************
 *                                   *
 *         Pneb::gh_fftb             *
 *                                   *
 *************************************/
void Pneb::gh_fftb(double *psi, double *psi_r)
{
   int n,done;
   int indx1,indx1n,shift1;
   int indx2,indx2n,shift2;


   n = neq[0]+neq[1];
   shift1 = 2*PGrid::npack(1);
   shift2 = n2ft3d;
   indx1  = indx1n = 0;
   indx2  = indx2n = 0;
   done = 0;
   while (!done)
   {
      if (indx1<n)
      {
         cr_pfft3b_queuein(1,psi+indx1n);
         indx1n += shift1;
         ++indx1;
      }
      if (cr_pfft3b_queuefilled() || (indx1>=n))
      {
         cr_pfft3b_queueout(1,psi_r+indx2n);
         indx2n += shift2;
         ++indx2;
      }
      done = ((indx1>=n) && (indx2>=n));
   }
}

/*************************************
 *                                   *
 *         Pneb::hr_aSumSqr          *
 *                                   *
 *************************************/
void Pneb::hr_aSumSqr(const double alpha, double *psir, double *dn)
{
   int n,ms,k,indx0,indx1;
   int one=1;
   int zero=0;
   int nsize = n2ft3d*ispin;
   double rzero = 0.0;

   //dcopy_(&nsize,&rzero,&zero,dn,&one);
   std::memset(dn,0,nsize*sizeof(double));

   indx0 = 0;
   indx1 = 0;
   for (ms=0; ms<ispin; ++ms)
   {
      for (n=0; n<(neq[ms]); ++n)
      {
         for (k=0; k<n2ft3d; ++k)
            dn[indx0+k] += alpha*psir[indx1+k]*psir[indx1+k];
         indx1 += n2ft3d;
      }
      indx0 += n2ft3d;
   }
   d3db::parall->Vector_SumAll(2,ispin*n2ft3d,dn);

}

/*************************************
 *                                   *
 *         Pneb::hhr_aSumMul         *
 *                                   *
 *************************************/
/*
   This routine calculates the perturbation density.

   dn12(r) = Sum_i (psi0_i(r)*psi1_i(r) + psi1_i(r)*psi0_i(r) )

   where

   psi_i(r) = psi0_i(r) + lmbda*psi1_i(r)
   dn12(r) = n(r) = Sum_i (psi_i(r)*psi_i(r)) 
                  = Sum_i ( psi0_i(r)*psi0_i(r) 
                          + lmbda*(psi0_i(r)*psi1_i(r) + psi1_i(r)*psi0_i(r)) 
                          + lmbda^2*psi1_i(r)*psi1_i(r) )
*/
void Pneb::hhr_aSumMul(const double alpha, const double *psir0, const double *psir1, double *dn12)
{
   int n,ms,k,indx0,indx1;
   int one=1;
   int zero=0;
   int nsize = n2ft3d*ispin;
   double rzero = 0.0;

   std::memset(dn12,0,nsize*sizeof(double));

   indx0 = 0;
   indx1 = 0;
   for (ms=0; ms<ispin; ++ms)
   {
      for (n=0; n<(neq[ms]); ++n)
      {
         for (k=0; k<n2ft3d; ++k)
            dn12[indx0+k] += alpha*(psir0[indx1+k]*psir1[indx1+k] + psir1[indx1+k]*psir0[indx1+k]);
         indx1 += n2ft3d;
      }
      indx0 += n2ft3d;
   }
   d3db::parall->Vector_SumAll(2,ispin*n2ft3d,dn12);

}



/*************************************
 *                                   *
 *      Pneb::ggm_sym_Multiply       *
 *                                   *
 *************************************/
void Pneb::ggm_sym_Multiply(double *psi1, double *psi2, double *hml)
{

   nwpw_timing_function ftimer(15);
   int ms,j,k,n,shift0,shift1,mshift0,mshift1;

   int one = 1;
   int ng  = 2*PGrid::npack(1);
   int ng0 = 2*PGrid::nzero(1);

   double rzero = 0.0;
   double rtwo  = 2.0;
   double rone =  1.0;
   double rmone = -1.0;

   if (parallelized)
   {
        std::ostringstream msg;
        msg << "NWPW Error: ggm_sym_Multiply() parallelized is NOT supported\n"
	    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
   }
   else
   {
      shift0  = 0;
      mshift0 = 0;
      for (ms=0; ms<ispin; ++ms)
      {
         n       = ne[ms];
         gdevice_TN1_dgemm(ng,n,rtwo,psi1+shift0,psi2+shift0,rzero,hml+mshift0);

         if (ng0>0)
         {
            shift1  = shift0;
            mshift1 = mshift0;
            for (k=1; k<=n; ++k)
            {
               //DGEMM_PWDFT((char *) "T",(char *) "N",k,one,ng,
               //        rtwo,
               //        &psi1[shift0], ng,
               //        &psi2[shift1],ng,
               //        rzero,
               //        &hml[mshift1],k);
                DGEMM_PWDFT((char *) "T",(char *) "N",k,one,ng0,
                       rmone,
                       psi1+shift0,ng,
                       psi2+shift1,ng,
                       rone,
                       hml+mshift1,k);
                shift1  += ng;
                mshift1 += n;
             }
          }
          for (k=0; k<n; ++k)
          for (j=k+1; j<n; ++j)
             hml[mshift0 + j + k*n] = hml[mshift0 + k + j*n];

         shift0  += ng*ne[0];
         mshift0 += ne[0]*ne[0];
      }
      d3db::parall->Vector_SumAll(1,ne[0]*ne[0] + ne[1]*ne[1],hml);
   }
}
/*************************************
 *                                   *
 *      Pneb::ffm_sym_Multiply       *
 *                                   *
 *************************************/
void Pneb::ffm_sym_Multiply(const int mb, double *psi1, double *psi2, double *hml)
{
   nwpw_timing_function ftimer(15);
   int ms,ms1,ms2,ishift2,j,k,n,shift0,shift1,mshift0,mshift1,nn;

   int one = 1;
   int ng  = 2*PGrid::npack(1);
   int ng0 = 2*PGrid::nzero(1);

   double rzero = 0.0;
   double rtwo  = 2.0;
   double rone =  1.0;
   double rmone = -1.0;

   if (parallelized)
   {
        std::ostringstream msg;
        msg << "NWPW Error: ffm_sym_Multiply() parallelized is NOT supported\n"
	    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
   }
   else
   {
      if (mb==-1)
      {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0];
          nn=ne[0]*ne[0]+ne[1]*ne[1];
          shift0  = 0;
          mshift0 = 0;
      }
      else
      {   ms1=mb; ms2=mb+1; ishift2=0;
          nn = ne[mb]*ne[mb];
          shift0  = mb*ne[0]*ng;
          mshift0 = 0;
      }
      for (ms=ms1; ms<ms2; ++ms)
      {
         n       = ne[ms];
         gdevice_TN1_dgemm(ng,n,rtwo,&psi1[shift0],&psi2[shift0],rzero,&hml[mshift0]);

         if (ng0>0)
         {
            shift1  = shift0;
            mshift1 = mshift0;
            for (k=1; k<=n; ++k)
            {
               //DGEMM_PWDFT((char *) "T",(char *) "N",k,one,ng,
               //        rtwo,
               //        &psi1[shift0],ng,
               //        &psi2[shift1],ng,
               //        rzero,
               //        &hml[mshift1],k);
               DGEMM_PWDFT((char *) "T",(char *) "N",k,one,ng0,
                       rmone,
                       &psi1[shift0],ng,
                       &psi2[shift1],ng,
                       rone,
                       &hml[mshift1],k);

               shift1  += ng;
               mshift1 += n;
            }
         }
         for (k=0; k<n; ++k)
         for (j=k+1; j<n; ++j)
            hml[mshift0 + j + k*n] = hml[mshift0 + k + j*n];

         shift0  += ng*ne[0];
         mshift0 += ne[0]*ne[0];
      }
      d3db::parall->Vector_SumAll(1,nn,hml);
   }
}

void Pneb::ffm_Multiply(const int mb, double *psi1, double *psi2, double *hml)
{
   nwpw_timing_function ftimer(15);
   int ms,ms1,ms2,ishift2,j,k,n,shift0,shift1,mshift0,mshift1,nn;

   int one = 1;
   int ng  = 2*PGrid::npack(1);
   int ng0 = 2*PGrid::nzero(1);

   double rzero = 0.0;
   double rtwo  = 2.0;
   double rone =  1.0;
   double rmone = -1.0;

   if (parallelized)
   {
        std::ostringstream msg;
        msg << "NWPW Error: ffm_sym_Multiply() parallelized is NOT supported\n"
            << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
   }
   else
   {
      if (mb==-1)
      {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0];
          nn=ne[0]*ne[0]+ne[1]*ne[1];
          shift0  = 0;
          mshift0 = 0;
      }
      else
      {   ms1=mb; ms2=mb+1; ishift2=0;
          nn = ne[mb]*ne[mb];
          shift0  = mb*ne[0]*ng;
          mshift0 = 0;
      }
      for (ms=ms1; ms<ms2; ++ms)
      {  
         n       = ne[ms];
         gdevice_TN1_dgemm(ng,n,rtwo,&psi1[shift0],&psi2[shift0],rzero,&hml[mshift0]);

         //DGEMM_PWDFT((char *) "T",(char *) "N",n,n,ng,
         //           rtwo,
         //           &psi1[shift0],ng,
         //           &psi2[shift0],ng,
         //           rzero,
         //           &hml[mshift0],n);
         if (ng0>0)
            DGEMM_PWDFT((char *) "T",(char *) "N",n,n,ng0,
                       rmone,
                       &psi1[shift0],ng,
                       &psi2[shift0],ng,
                       rone,
                       &hml[mshift0],n);
         shift0  += ng*ne[0];
         mshift0 += ne[0]*ne[0];
      }
      d3db::parall->Vector_SumAll(1,nn,hml);
   }
}

void Pneb::ffm3_sym_Multiply(const int mb, double *psi1, double *psi2,
			     double* s11, double* s21, double* s22)
{
   nwpw_timing_function ftimer(15);
   int ms,ms1,ms2,ishift2,j,k,n,shift0,shift1,mshift0,mshift1,nn;
   int one = 1;
   int ng  = 2*PGrid::npack(1);
   int ng0 = 2*PGrid::nzero(1);

   double rzero = 0.0;
   double rtwo  = 2.0;
   double rone =  1.0;
   double rmone = -1.0;

   if (parallelized)
   {
        std::ostringstream msg;
        msg << "NWPW Error: ffm3_sym_Multiply() parallelized is NOT supported\n"
	    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
   }
   else
   {
      if (mb==-1)
      {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0];
          nn=ne[0]*ne[0]+ne[1]*ne[1];
          shift0  = 0;
          mshift0 = 0;
      }
      else
      {   ms1=mb; ms2=mb+1; ishift2=0;
          nn = ne[mb]*ne[mb];
          shift0  = mb*ne[0]*ng;
          mshift0 = 0;
      }
      for (ms=ms1; ms<ms2; ++ms)
      {
         n       = ne[ms];

         gdevice_TN3_dgemm(ng, n, rtwo, &psi1[shift0], &psi2[shift0], rzero, &s11[mshift0], &s21[mshift0], &s22[mshift0]);

         if (ng0>0)
         {
            shift1  = shift0;
            mshift1 = mshift0;
            for (k=1; k<=n; ++k)
            {
		DGEMM_PWDFT((char *) "T",(char *) "N",k,one,ng0,
			    rmone,
			    &psi1[shift0],ng,
			    &psi1[shift1],ng,
			    rone,
			    &s11[mshift1],k);
		DGEMM_PWDFT((char *) "T",(char *) "N",k,one,ng0,
			    rmone,
			    &psi1[shift0],ng,
			    &psi2[shift1],ng,
			    rone,
			    &s21[mshift1],k);
		DGEMM_PWDFT((char *) "T",(char *) "N",k,one,ng0,
			    rmone,
			    &psi2[shift0],ng,
			    &psi2[shift1],ng,
			    rone,
			    &s22[mshift1],k);
               shift1  += ng;
               mshift1 += n;
            }
         }

         for (k=0; k<n; ++k) {
	   for (j=k+1; j<n; ++j) {
	     s11[mshift0 + j + k*n] = s11[mshift0 + k + j*n];
	     s21[mshift0 + j + k*n] = s21[mshift0 + k + j*n];
	     s22[mshift0 + j + k*n] = s22[mshift0 + k + j*n];
	   }
	 }

         shift0  += ng*ne[0];
         mshift0 += ne[0]*ne[0];
      }
      d3db::parall->Vector_SumAll(1,nn,s11);
      d3db::parall->Vector_SumAll(1,nn,s21);
      d3db::parall->Vector_SumAll(1,nn,s22);
   }
}

void Pneb::fmf_Multiply(const int mb, double *psi1, double *hml, double alpha, double *psi2, double beta)
{
   nwpw_timing_function ftimer(16);
   int ms,ms1,ms2,n,shift1,mshift1,ishift2;
   int ng  = 2*PGrid::npack(1);

   if (parallelized)
   {
        std::ostringstream msg;
        msg << "NWPW Error: fmf_Multiply() parallelized is NOT supported\n"
	    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
   }
   else
   {
      if (mb==-1)
      {  ms1=0; ms2=ispin; ishift2=ne[0]*ne[0];
         shift1  = 0;
         mshift1 = 0;
       }
      else
      {   ms1=mb; ms2=mb+1; ishift2=0;
          shift1  = mb*ne[0]*ng;
          mshift1 = 0;
      }
      for (ms=ms1; ms<ms2; ++ms)
      {
         n = ne[ms];

	 // DGEMM_PWDFT((char *) "N",(char *) "N",ng,n,n,
         //	     alpha,
	 //	     &psi1[shift1],ng,
	 //	     &hml[mshift1],n,
	 //	     beta,
	 //	     &psi2[shift1],ng);
         gdevice_NN_dgemm(ng,n,alpha,&psi1[shift1],&hml[mshift1],beta,&psi2[shift1]);

         shift1  += ne[0]*ng;
         mshift1 += ishift2;
      }
   }
}



void Pneb::ggm_SVD(double *A, double *U, double *S, double *V)
{
   int n,indx;
   double *tmp2 = new (std::nothrow) double [neq[0] + neq[1]]();

   /* generate V and Sigma^2 */
   ggm_sym_Multiply(A,A,V);
   m_diagonalize(V,S);

   /* generate U*Sigma */
   fmf_Multiply(-1,A,V,1.0,U,0.0);

   /* normalize U*sigma */
   indx = 0;
   for (n=0; n<(neq[0]+neq[1]); ++n)
   {
      tmp2[n] = PGrid::cc_pack_idot(1,&U[indx],&U[indx]);
      indx += 2*PGrid::npack(1);
   }
   d3db::parall->Vector_SumAll(1,neq[0]+neq[1],tmp2);

   for (n=0; n<(neq[0]+neq[1]); ++n)
      tmp2[n] = 1.0/sqrt(tmp2[n]);

   indx = 0;
   for (n=0; n<(neq[0]+neq[1]); ++n)
   {
      PGrid::c_pack_SMul(1,tmp2[n],&U[indx]);
      indx += 2*PGrid::npack(1);
   }

   /* calculated sqrt(S^2) */
   for (n=0; n<(neq[0]+neq[1]); ++n)
   {
     if (S[n]<0.0) S[n] = fabs(S[n]);
     S[n] = sqrt(S[n]);
   }

   delete [] tmp2;
}


void Pneb::m_scal(double alpha, double *hml)
{
  int one = 1;
  int nsize = ne[0]*ne[0] + ne[1]*ne[1];

  DSCAL_PWDFT(nsize,alpha,hml,one);
}

double Pneb::m_trace(double *hml)
{
   int ms,i;
   int mshift = 0;
   double sum = 0.0;
   for (ms=0; ms<ispin; ++ms)
   {
      for (i=0; i<ne[ms]; ++i)
         sum += hml[i+i*ne[ms]+mshift];
      mshift += ne[0]*ne[0];
   }
   return sum;
}

void Pneb::m_diagonalize(double *hml, double *eig)
{
    nwpw_timing_start(17);
    int shift1,shift2;
    int n,ierr;

    if (parallelized)
    {
        std::ostringstream msg;
        msg << "NWPW Error: m_diagonalize() parallelized is NOT supported\n"
	    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
    }
    else
    {
        //double *xmp1 = new (std::nothrow) double[nn]();
        int nn  = ne[0]*ne[0]+14;
        double xmp1[nn];

        shift1 = 0;
        shift2 = 0;
        for (int ms=0; ms<ispin; ++ms)
        {
            n = ne[ms];

            //eigen_(&n,&n,&hml[shift2],&eig[shift1],xmp1,&ierr);

            d3db::parall->Barrier();

            //ierr=LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, &hml[shift2], n, &eig[shift1]);
            EIGEN_PWDFT(n,hml+shift2,eig+shift1,xmp1,nn,ierr);
            if (ierr != 0) throw std::runtime_error(std::string("NWPW Error: EIGEN_PWDFT failed!"));

            eigsrt(eig+shift1,hml+shift2,n);
            shift1 += ne[0];
            shift2 += ne[0]*ne[0];
        }

        //delete [] xmp1;
    }
    nwpw_timing_end(17);
}

void Pneb::m_scale_s22_s21_s11(const int mb, const double dte, double *s22, double *s21, double *s11)
{
   int j,k,ms,ms1,ms2,ishift2,indx0,indx,indxt;

   if (mb==-1)
   {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0];}
   else
   {   ms1=mb; ms2=mb+1; ishift2=0;}

   for (ms=ms1; ms<ms2; ++ms)
   {
      indx0 = ms*ishift2;
      for (k=0; k<ne[ms]; ++k)
      {
          s22[indx0] = (1.0-s22[indx0])*(0.5/dte);
          s21[indx0] = (1.0-s21[indx0])*(0.5);
          s11[indx0] *= -0.5*dte;

          indx  = indx0 + 1;
          indxt = indx0 + ne[ms];
          for (j=(k+1); j<ne[ms]; ++j)
          {
             s22[indx]  *= (-0.5/dte);
             s22[indxt] *= (-0.5/dte);
             s21[indx]  *= -0.5;
             s21[indxt] *= -0.5;
             s11[indx]  *= -0.5*dte;
             s11[indxt] *= -0.5*dte;

             indx  += 1;
             indxt += ne[ms];
          }
          indx0 += (ne[ms]+1);
      }
   }
}

void Pneb::mm_SCtimesVtrans(const int mb, const double t, double *S, double *Vt, 
                            double *A, double *B, double *SA, double *SB)
{
   nwpw_timing_function ftimer(19);
   int ms,n,ms1,ms2,ishift2,shift2,ishift1,shift1,nj;
   int j,k,indx1,indx2;
   if (mb==-1)
   {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0]; ishift1=ne[0]; nj = ne[0]+ne[1];}
   else
   {   ms1=mb; ms2=mb+1; ishift2=0; ishift1 = 0; nj = ne[mb];}
   for (j=0; j<nj; ++j)
   {
       SA[j] = cos(S[j]*t);
       SB[j] = sin(S[j]*t);
   }
   for (ms=ms1; ms<ms2; ++ms)
   {
      shift1 = ms*ishift1;
      shift2 = ms*ishift2;
      for (k=0; k<ne[ms]; ++k)
      {
         indx1 = shift1;
         indx2 = shift2 + k*ne[ms];
         for (j=0; j<ne[ms]; ++j)
         {
            A[indx2] = SA[indx1]*Vt[indx2];
            B[indx2] = SB[indx1]*Vt[indx2];
            ++indx1;
            ++indx2;
         }
      }

   }
}

void Pneb::mm_SCtimesVtrans2(const int mb, const double t, double *S, double *Vt,
                             double *A, double *B, double *SA, double *SB)
{
   nwpw_timing_function ftimer(19);
   int ms,n,ms1,ms2,ishift2,shift2,ishift1,shift1,nj;
   int j,k,indx1,indx2;
   if (mb==-1)
   {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0]; ishift1=ne[0]; nj = ne[0]+ne[1];}
   else
   {   ms1=mb; ms2=mb+1; ishift2=0; ishift1 = 0; nj = ne[mb];}
   for (j=0; j<nj; ++j)
   {
       SA[j] = S[j]*sin(S[j]*t);
       SB[j] = S[j]*cos(S[j]*t);
   }
   for (ms=ms1; ms<ms2; ++ms)
   {
      shift1 = ms*ishift1;
      shift2 = ms*ishift2;
      for (k=0; k<ne[ms]; ++k)
      {
         indx1 = shift1;
         indx2 = shift2 + k*ne[ms];
         for (j=0; j<ne[ms]; ++j)
         {
            A[indx2] = SA[indx1]*Vt[indx2];
            B[indx2] = SB[indx1]*Vt[indx2];
            ++indx1;
            ++indx2;
         }
      }

   }
}

void Pneb::mm_SCtimesVtrans3(const int mb, const double t, double *S, double *Vt,
                             double *A, double *B, double *SA, double *SB)
{
   nwpw_timing_function ftimer(19);
   int ms,n,ms1,ms2,ishift2,shift2,ishift1,shift1,nj;
   int j,k,indx1,indx2;
   if (mb==-1)
   {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0]; ishift1=ne[0]; nj = ne[0]+ne[1];}
   else
   {   ms1=mb; ms2=mb+1; ishift2=0; ishift1 = 0; nj = ne[mb];}
   for (j=0; j<nj; ++j)
   {
       SA[j] = sin(S[j]*t);
       SB[j] = 1.0-cos(S[j]*t);
   }
   for (ms=ms1; ms<ms2; ++ms)
   {
      shift1 = ms*ishift1;
      shift2 = ms*ishift2;
      for (k=0; k<ne[ms]; ++k)
      {
         indx1 = shift1;
         indx2 = shift2 + k*ne[ms];
         for (j=0; j<ne[ms]; ++j)
         {
            A[indx2] = SA[indx1]*Vt[indx2];
            B[indx2] = SB[indx1]*Vt[indx2];
            ++indx1;
            ++indx2;
         }
      }

   }
}


void Pneb::mmm_Multiply(const int mb, double *a, double *b, double alpha, double *c, double beta)
{
   nwpw_timing_function ftimer(18);
   int ms,n,ms1,ms2,ishift2,shift2;
   if (mb==-1)
   {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0];}
   else
   {   ms1=mb; ms2=mb+1; ishift2=0;}
   for (ms=ms1; ms<ms2; ++ms)
   {
      n = ne[ms];
      if (n>0)
      {
         shift2 = ms*ishift2;

         DGEMM_PWDFT((char *) "N",(char *) "N",n,n,n,
		     alpha,
		     a+shift2, n,
		     b+shift2, n,
		     beta,
		     c+shift2, n);
      }
   }
}


void Pneb::mmm_Multiply2(const int mb, double *a, double *b, double alpha, double *c, double beta)
{
   nwpw_timing_function ftimer(18);
   int ms,n,ms1,ms2,ishift2,shift2;
   if (mb==-1)
   {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0];}
   else
   {   ms1=mb; ms2=mb+1; ishift2=0;}
   for (ms=ms1; ms<ms2; ++ms)
   {
      n = ne[ms];
      if (n>0)
      {
         shift2 = ms*ishift2;

         DGEMM_PWDFT((char *) "T",(char *) "N",n,n,n,
                     alpha,
                     a+shift2, n,
                     b+shift2, n,
                     beta,
                     c+shift2, n);
      }
   }
}

void Pneb::mm_transpose(const int mb, double *a, double *b)
{
   int i,j,indx,indxt;
   int ms,n,ms1,ms2,ishift2,shift2;
   if (mb==-1)
   {   ms1=0; ms2=ispin; ishift2=ne[0]*ne[0];}
   else
   {   ms1=mb; ms2=mb+1; ishift2=0;}
   for (ms=ms1; ms<ms2; ++ms)
   {
      shift2 = ms*ishift2;
      for (j=0; j<ne[ms]; ++j)
      for (i=0; i<ne[ms]; ++i)
      {
          indx  = i + j*ne[ms] + shift2;
          indxt = j + i*ne[ms] + shift2;
          b[indx] = a[indxt];
      }
   }
}

#define ITERLMD         220
#define CONVGLMD        1e-15
#define CONVGLMD2       1e-12

// Lagrange multiplier (expensive method)
void Pneb::ggm_lambda(double dte, double *psi1, double *psi2, double *lmbda)
{
    nwpw_timing_function ftimer(3);

    int one=1;
    double rmone = -1.0;
    double adiff = 0.0;

    for (int ms=0; ms<ispin; ++ms) {

	int nn = m_size(ms);

	ffm3_sym_Multiply(ms, psi1, psi2, s11, s21, s22);
	m_scale_s22_s21_s11(ms, dte, s22, s21, s11);

        int jj;
	int ii   = 0;
	int done = 0;

	//DCOPY_PWDFT(nn, s21, one, s12, one);
	//DCOPY_PWDFT(nn, s22, one, sa0, one);
   std::memcpy(s12,s21,nn*sizeof(double));
   std::memcpy(sa0,s22,nn*sizeof(double));

	while ((!done) && ((ii++)<ITERLMD)) {
	    //DCOPY_PWDFT(nn, s22, one, sa1, one);
       std::memcpy(sa1,s22,nn*sizeof(double));

	    //mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
	    //mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
	    //mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
	    //mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
            gdevice_MM6_dgemm(ne[ms],s12,s12,s11,sa0,sa1,st1);

	    //DCOPY_PWDFT(nn, sa1, one, st1, one);
       std::memcpy(st1,sa1,nn*sizeof(double));
	    DAXPY_PWDFT(nn, rmone, sa0, one, st1, one);
	    jj = IDAMAX_PWDFT(nn, st1, one);
	    adiff = fabs(st1[IDAMAX_PWDFT(nn, st1, one) - 1]);

	    if (adiff < CONVGLMD)
		   done = 1;
	    else
          std::memcpy(sa0,sa1,nn*sizeof(double)); //DCOPY_PWDFT(nn, sa1, one, sa0, one);
	}
	//printf("ierr=10 check nn=%d jj=%d adiff=%le ii=%d done=%d\n",nn,jj,adiff,ii,done);

	if (adiff>CONVGLMD2) {
	    if (!done) printf("ierr=10 adiff=%le\n",adiff);
	}

	//DCOPY_PWDFT(nn, sa1, one, &lmbda[ms*ne[0]*ne[0]], one);
   std::memcpy(&lmbda[ms*ne[0]*ne[0]],sa1,nn*sizeof(double));
        //std::memcpy(lmbda+ms*ne[0]*ne[0],sa1,nn*sizeof(double));

    } // for loop - ms

    /* correction due to contraint */
    fmf_Multiply(-1, psi1, lmbda, dte, psi2, 1.0);
}

/********************************
 *                              *
 *        Pneb::g_ortho         *
 *                              *
 ********************************/
/*
   Performs a Gram-Schmidt orthogonalization on psi
*/
void Pneb::g_ortho(double *psi)
{
   int indxj,indxk,ishift;
   double w;
   if (parallelized)
   {
        std::ostringstream msg;
        msg << "NWPW Error: g_ortho() parallelized is NOT supported\n"
	    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
        throw(std::runtime_error(msg.str()));
   }
   else
   {
      for (auto ms=0; ms<ispin; ++ms)
      {
         ishift = ms*ne[0]*2*PGrid::npack(1);
         for (auto k=ne[ms]-1; k>=0;  --k)
         {
            indxk = 2*PGrid::npack(1)*k + ishift;
            w = PGrid::cc_pack_dot(1,&psi[indxk],&psi[indxk]);

            w = 1.0/sqrt(w);
            PGrid::c_pack_SMul(1,w,&psi[indxk]);

            for (auto j=k-1; j>=0; --j)
            {
               indxj = 2*PGrid::npack(1)*j  + ishift;
               w = -PGrid::cc_pack_dot(1,&psi[indxk],&psi[indxj]);
               PGrid::cc_pack_daxpy(1,w,&psi[indxk],&psi[indxj]);
            }
         }
      }
   }
}

}


