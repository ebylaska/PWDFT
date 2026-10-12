/* cKinetic.C -
   Author - Eric Bylaska
*/

/*


#include        <cmath>
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>
#include	<string>
*/

#include "cKinetic.hpp"

namespace pwdft {

/*********************************************
 *                                           *
 *     cKinetic_Operator::cKinetic_Operator  *
 *                                           *
 *********************************************/
cKinetic_Operator::cKinetic_Operator(Cneb *mygrid) 
{
   mycneb = mygrid;

   int nbrillq    = mycneb->nbrillq;
   int npack1_max = mycneb->npack1_max();
   tg = new double[nbrillq*npack1_max];
   //double *tmp = new double[mygrid->nfft3d];
   std::memset(tg,0.0,nbrillq*npack1_max);
 
   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {
      int npack1     = mycneb->npack(1+nbq);
      double *Gpackx = mycneb->Gpackxyz(1+nbq,0);
      double *Gpacky = mycneb->Gpackxyz(1+nbq,1);
      double *Gpackz = mycneb->Gpackxyz(1+nbq,2);

      double *kv = mycneb->pbrill_kvector(nbq);
      double *tmp_tg = tg + nbq*npack1_max;
 
      for (auto k=0; k<npack1; ++k) 
      {
         double gx = Gpackx[k] + kv[0];
         double gy = Gpacky[k] + kv[1];
         double gz = Gpackz[k] + kv[2];
         double gg = gx*gx + gy*gy + gz*gz;
         tmp_tg[k] = -0.5 * gg;
      }
      //mycneb->tt_pack_copy(1, tmp, tg + nbq*npack1);
   }
   //delete[] tmp;
}


/*******************************************
 *                                         *
 *         cKinetic_Operator::ke           *
 *                                         *
 *******************************************/
void cKinetic_Operator::ke(const double *psi, double *tpsi) 
{
   int nbqsize    = mycneb->nbrillq;
   int nsize      = mycneb->neq[0] + mycneb->neq[1];
   int npack1_max = mycneb->npack1_max();
   int shift1 = 2*npack1_max;
   int indx1n = 0;

   for (auto nbq=0; nbq<nbqsize; ++nbq)
   {
      int npack1  = mycneb->npack(1+nbq);
      double *tmp_tg = tg + nbq*npack1_max;

      for (auto n=0; n<nsize; ++n)
      {
         const double *tmp_psi  = psi  + indx1n;
         double       *tmp_tpsi = tpsi + indx1n;
         int k1 = 0;
         int k2 = 1;
         for (auto k=0; k<npack1; ++k) 
         {
           //tmp_tpsi[k1] += tmp_tg[k]*tmp_psi[k1];
           //tmp_tpsi[k2] += tmp_tg[k]*tmp_psi[k2];
           tmp_tpsi[k1] = tmp_tg[k]*tmp_psi[k1];
           tmp_tpsi[k2] = tmp_tg[k]*tmp_psi[k2];
           k1 += 2;
           k2 += 2;
         }
         indx1n += shift1;
      }
   }
}


/*******************************************
 *                                         *
 *        cKinetic_Operator::ke_ave        *
 *                                         *
 *******************************************/
double cKinetic_Operator::ke_ave(const double *psi) 
{
 
   int nbqsize    = mycneb->nbrillq;
   int nsize      = mycneb->neq[0] + mycneb->neq[1];
   int npack1_max = mycneb->npack1_max();
   int shift1 = 2*npack1_max;
   int indx1n = 0;
 
   double ave = 0.0;
   for (auto nbq=0; nbq<nbqsize; ++nbq)
   {
      int npack1     = mycneb->npack(1+nbq);
      double weight  = mycneb->pbrill_weight(nbq);
      double *tmp_tg = tg + nbq*npack1_max;

      for (auto n=0; n<nsize; ++n) 
      {
         const double *tmp_psi = psi + indx1n;
         int k1 = 0; 
         int k2 = 1;
         for (auto k=0; k<npack1; ++k) 
         {
            ave += weight*tmp_tg[k]*(tmp_psi[k1]*tmp_psi[k1] + tmp_psi[k2]*tmp_psi[k2]);
            k1 += 2;
            k2 += 2;
         }
         indx1n += shift1;
      }
   }

   ave = mycneb->c3db::parall->SumAll(0, ave);

   if (mycneb->ispin == 1)
      ave *= 2.0;
   ave = -ave;
 
   return ave;
}


double cKinetic_Operator::ke_ave(const double *psi, const double *occ)
{

   int nbqsize    = mycneb->nbrillq;
   int nsize      = mycneb->neq[0] + mycneb->neq[1];
   int npack1_max = mycneb->npack1_max();
   int shift1 = 2*npack1_max;
   int indx1n = 0;

   double ave = 0.0;
   for (auto nbq=0; nbq<nbqsize; ++nbq)
   {
      int npack1     = mycneb->npack(1+nbq);
      double weight  = mycneb->pbrill_weight(nbq);
      double *tmp_tg = tg + nbq*npack1_max;

      for (auto n=0; n<nsize; ++n)
      {
         const double *tmp_psi = psi + indx1n;
         int k1 = 0;
         int k2 = 1;
         for (auto k=0; k<npack1; ++k)
         {
            ave += weight*tmp_tg[k]*(tmp_psi[k1]*tmp_psi[k1] + tmp_psi[k2]*tmp_psi[k2])*occ[n+nbq*nsize];
            k1 += 2;
            k2 += 2;
         }
         indx1n += shift1;
      }
   }

   ave = mycneb->c3db::parall->SumAll(0, ave);

   if (mycneb->ispin == 1)
      ave *= 2.0;
   ave = -ave;

   return ave;
}


/*******************************************
 *                                         *
 *         cKinetic_Operator::ke_orb       *
 *                                         *
 *******************************************/
void cKinetic_Operator::ke_orb(const int nbq1, const double *psi, double *tpsi) 
{
   int npack1_max = mycneb->npack1_max();
   double *tmp_tg = tg + (nbq1-1)*npack1_max;

   int npack1  = mycneb->npack(nbq1);

   int k1 = 0;
   int k2 = 1;
   for (auto k=0; k<npack1; ++k)
   {
     //tmp_tpsi[k1] += tmp_tg[k]*tmp_psi[k1];
     //tmp_tpsi[k2] += tmp_tg[k]*tmp_psi[k2];
     tpsi[k1] = tmp_tg[k]*psi[k1];
     tpsi[k2] = tmp_tg[k]*psi[k2];
     k1 += 2;
     k2 += 2;
   }
}


/***********************************************
 *                                             *
 *       cKinetic_Operator::ke_precondition     *
 *                                             *
 ***********************************************/
// **** My preconditioner ****
void cKinetic_Operator::ke_precondition(const int nbq1, const double Ep, const int neall, double *psi, double *tpsi) 
{
   int npack1  = (mycneb->npack(nbq1));
   int npack1_max = mycneb->npack1_max();
   int npack2  = 2*npack1_max;
   double *tmp    = new double[npack2];
   double *tg1    = tg + (nbq1-1);
   int k1 = 0;
   int k2 = 1;
   for (auto n=0; n<neall; ++n)
   {
      double *worb = psi + n*npack2;
      mycneb->tcc_pack_Mul(nbq1,tg1,worb,tmp);
      double sum = mycneb->cc_pack_dot(nbq1,tmp,worb);
      for (auto k=0; k<npack1; ++k) 
      {
         double x = tg1[k];
         x =  x*(worb[k1]*worb[k1] + worb[k2]*worb[k2])/sum;
         double cm = 27.00+(18.00+(12.00+8.00*x)*x)*x;
         double dm = (cm + 16.00* x*x*x*x);
         cm = cm/dm;
         tpsi[k1] *= cm;
         tpsi[k2] *= cm;
         k1 += 2;
         k2 += 2;
      }
   }
   delete [] tmp;
}


/*******************************************
 *                                         *
 *         cKinetic_Operator::euv          *
 *                                         *
 *******************************************/
//
// Kinetic contribution to the lattice stress tensor.
//
// Computes
//
//     Aus(u,s) = Σ_k w_k Σ_n occ_n Σ_G (G+k)_u (G+k)_s |ψ_{n,G+k}|²
//
// then symmetrizes Aus (Aus is built symmetric by construction here,
// since the summand is symmetric in u,s), reduces across MPI ranks,
// doubles if spin-restricted, and finally contracts with the reciprocal
// metric
//
//     hm(i,j) = (1/2π) * lattice->unitg(i,j)
//
// to produce
//
//     euv(u,v) = -Σ_s Aus(u,s) * hm(s,v).
//
// The occupation array, if provided, is indexed as occ[n + nbq*nsize]
// with n in [0, nsize).

namespace {

inline int gxyz_component(int u, double gx, double gy, double gz)
{
   switch (u) {
      case 0: return 0;  // unused; kept for clarity
      default: return 0;
   }
}

} // namespace

void cKinetic_Operator::euv(const double *psi, double euv[3][3])
{
   int nbqsize    = mycneb->nbrillq;
   int nsize      = mycneb->neq[0] + mycneb->neq[1];
   int npack1_max = mycneb->npack1_max();
   int shift1     = 2*npack1_max;

   // Accumulate the symmetric second-moment tensor in the (G+k) basis.
   double Aus[3][3] = {{0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0}};

   int indx1n = 0;
   for (int nbq = 0; nbq < nbqsize; ++nbq)
   {
      const int     npack1 = mycneb->npack(1 + nbq);
      const double  weight = mycneb->pbrill_weight(nbq);
      const double *Gpackx = mycneb->Gpackxyz(1 + nbq, 0);
      const double *Gpacky = mycneb->Gpackxyz(1 + nbq, 1);
      const double *Gpackz = mycneb->Gpackxyz(1 + nbq, 2);
      const double *kv     = mycneb->pbrill_kvector(nbq);

      const double kx = kv[0];
      const double ky = kv[1];
      const double kz = kv[2];

      for (int u = 0; u < 3; ++u)
      {
         for (int s = u; s < 3; ++s)
         {
            double asum = 0.0;

            for (int n = 0; n < nsize; ++n)
            {
               const double *tmp_psi = psi + indx1n + n*shift1;

               double orb_sum = 0.0;
               for (int k = 0; k < npack1; ++k)
               {
                  const double gx = Gpackx[k] + kx;
                  const double gy = Gpacky[k] + ky;
                  const double gz = Gpackz[k] + kz;

                  const double g_u = (u == 0) ? gx : (u == 1) ? gy : gz;
                  const double g_s = (s == 0) ? gx : (s == 1) ? gy : gz;

                  const double re = tmp_psi[2*k    ];
                  const double im = tmp_psi[2*k + 1];

                  orb_sum += g_u * g_s * (re*re + im*im);
               }
               asum += orb_sum;
            }
            Aus[u][s] += weight * asum;
         }
      }
      indx1n += nsize * shift1;
   }

   // Symmetrize (the diagonal and upper triangle are already correct).
   for (int u = 0; u < 3; ++u)
      for (int s = u + 1; s < 3; ++s)
         Aus[s][u] = Aus[u][s];

   // Reduce across ranks.  SumAll takes and returns the summed value.
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         Aus[i][j] = mycneb->c3db::parall->SumAll(0, Aus[i][j]);

   if (mycneb->ispin == 1)
      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            Aus[i][j] *= 2.0;

   // Reciprocal metric * 1/(2π).
   const double pi   = 4.0 * std::atan(1.0);
   const double scal = 1.0 / (2.0 * pi);

   double hm[3][3];
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         hm[i][j] = scal * mycneb->lattice->unitg(i, j);

   // euv = -Aus * hm  (matrix product).
   for (int u = 0; u < 3; ++u)
      for (int v = 0; v < 3; ++v)
      {
         double sum = 0.0;
         for (int s = 0; s < 3; ++s)
            sum += Aus[u][s] * hm[s][v];
         euv[u][v] = -sum;
      }
}



void cKinetic_Operator::euv(const double *psi, const double *occ,
                            double euv[3][3])
{
   int nbqsize    = mycneb->nbrillq;
   int nsize      = mycneb->neq[0] + mycneb->neq[1];
   int npack1_max = mycneb->npack1_max();
   int shift1     = 2*npack1_max;

   double Aus[3][3] = {{0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0},
                       {0.0, 0.0, 0.0}};

   int indx1n = 0;
   for (int nbq = 0; nbq < nbqsize; ++nbq)
   {
      const int     npack1 = mycneb->npack(1 + nbq);
      const double  weight = mycneb->pbrill_weight(nbq);
      const double *Gpackx = mycneb->Gpackxyz(1 + nbq, 0);
      const double *Gpacky = mycneb->Gpackxyz(1 + nbq, 1);
      const double *Gpackz = mycneb->Gpackxyz(1 + nbq, 2);
      const double *kv     = mycneb->pbrill_kvector(nbq);

      const double kx = kv[0];
      const double ky = kv[1];
      const double kz = kv[2];

      const double *occ_k = occ + nbq*nsize;

      for (int u = 0; u < 3; ++u)
      {
         for (int s = u; s < 3; ++s)
         {
            double asum = 0.0;

            for (int n = 0; n < nsize; ++n)
            {
               const double occ_n = occ_k[n];
               if (occ_n == 0.0) continue;

               const double *tmp_psi = psi + indx1n + n*shift1;

               double orb_sum = 0.0;
               for (int k = 0; k < npack1; ++k)
               {
                  const double gx = Gpackx[k] + kx;
                  const double gy = Gpacky[k] + ky;
                  const double gz = Gpackz[k] + kz;

                  const double g_u = (u == 0) ? gx : (u == 1) ? gy : gz;
                  const double g_s = (s == 0) ? gx : (s == 1) ? gy : gz;

                  const double re = tmp_psi[2*k    ];
                  const double im = tmp_psi[2*k + 1];

                  orb_sum += g_u * g_s * (re*re + im*im);
               }
               asum += occ_n * orb_sum;
            }
            Aus[u][s] += weight * asum;
         }
      }
      indx1n += nsize * shift1;
   }

   for (int u = 0; u < 3; ++u)
      for (int s = u + 1; s < 3; ++s)
         Aus[s][u] = Aus[u][s];

   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         Aus[i][j] = mycneb->c3db::parall->SumAll(0, Aus[i][j]);

   if (mycneb->ispin == 1)
      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            Aus[i][j] *= 2.0;

   const double pi   = 4.0 * std::atan(1.0);
   const double scal = 1.0 / (2.0 * pi);

   double hm[3][3];
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         hm[i][j] = scal * mycneb->lattice->unitg(i, j);

   for (int u = 0; u < 3; ++u)
      for (int v = 0; v < 3; ++v)
      {
         double sum = 0.0;
         for (int s = 0; s < 3; ++s)
            sum += Aus[u][s] * hm[s][v];
         euv[u][v] = -sum;
      }
}

} // namespace pwdft
