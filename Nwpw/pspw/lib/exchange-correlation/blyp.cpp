/* blyp.cpp -
   Author - Eric Bylaska
*/

#include	<cmath>
#include	"pbe96.hpp"

/****************************************************
 *                                                  *
 *             Becke_smalln_correction              *
 *                                                  *
 ****************************************************/
/*
   This routine does small n (larger gradient) correction in which
 the function is slowly switched to PBE when (fdnx-xe)>0.
*/

#define thrd		(1.00/3.00)
#define two_thrd	1.25992104989487319066
#define frthrd		(4.00/3.00)

static void Becke_smalln_correction(double n, double n_thrd, double fac, double beta, double lda_c,
                                    double chi, double chi2, double chiSQ, double K, double F1, double F2,
                                    double pbex, double dfpbednx, double dfpbedagrx,
                                    double *xe, double *fdnx, double *fdagrx)
{

   /* local variables */
   double nf,nf_thrd;
   double x,s,dsdx,tmp1,tmp2,dchi,dF1,dF2;
   double dfdnxdagr,ddf0,dxdn,dxdagr,dsdn,dsdagr,xeold;

   double pi = 4.0*atan(1.0);

   nf = n/fac;
   nf_thrd = n_thrd/pow(fac,thrd);

   x = 0.850*(*fdnx-*xe)/(nf_thrd);
   s = tanh(x);
   dsdx = 0.00;
   if (x<100.00) dsdx = pow((1.00/cosh(x)),2);
   
   tmp1 = 6.00*beta*chi/chiSQ;
   tmp2 = (1.00+chi*K);
   dchi = -frthrd*chi/nf;
   dF1 = -frthrd*chi2*F2/nf;
   dF2 = ((-tmp1/(1.00+chi2)+K)/tmp2
         -2.00*F2*(tmp1+K))*dchi/tmp2;
   dfdnxdagr = -beta*(chi*dF2+F2*dchi);
   ddf0 = (thrd/nf)*(*fdnx)
       -frthrd*beta*nf_thrd
        *(dF1-2.00*chi*F2*dchi-chi2*dF2);

   dxdn = 0.850*ddf0/nf_thrd - frthrd*x/nf;
   dxdagr = 0.850*(dfdnxdagr - (*fdagrx)/nf)/nf_thrd;

   dsdn   = dsdx * dxdn;
   dsdagr = dsdx*dxdagr;

   xeold = *xe;
   *xe   = (1.00-s)*xeold   + s*pbex;
   *fdnx = (1.00-s)*(*fdnx)    + s*dfpbednx
         + dsdn  *nf*(-xeold + pbex);
   *fdagrx = (1.00-s)*(*fdagrx)  + s*dfpbedagrx 
           + dsdagr*nf*(-xeold + pbex);
}

/* tolerance constants */
#define ETA	1.0e-20
#define ETA2	2.0e-20
#define DNS_CUT	1.0e-18
#define tolrho	2.0e-11
#define minagr	1.0e-12

/* LYP parameters */
#define a	0.04918
#define b	0.132
#define c	0.2533
#define d	0.349
#define Cf	2.87123400018819108225

/* collated LYP parameters */
#define ho1	((19.00/36.00)*a*b)
#define ho2	(( 7.00/24.00)*a*b)
#define ho3	((59.00/72.00)*a*b)
#define thrd_d	(thrd*d)
#define abCf	(a*b*Cf)

/* Becke88 parameters */
#define beta	0.0042 
#define lda_c	0.93052573634910018540

/* other constants */
#define thrd	(1.00/3.00)
#define twthrd	(2.00/3.00)
#define frthrd	(4.00/3.00)
#define fvthrd	(5.00/3.00)
#define snthrd	(7.00/3.00)
#define etthrd	(8.00/3.00)


/**************************************************
 *                                                *
 *           gen_BLYP_BW_unrestricted             *
 *                                                *
 **************************************************/
void gen_BLYP_BW_unrestricted(const int n2ft3d,
                              double *dn_in, double *agr_in,
                              const double x_parameter, const double c_parameter,
                              double *xce, double *fn, double *fdn)
{
   /* local variables */
   double nup,ndn,n,agrup,agrdn,agr;
   double agr2,agrup2,agrdn2;
   double n2,nup2,ndn2,nup_thrd,ndn_thrd,sdup,sddn;
   double gamma,gamma_u,gamma_uu;
   double gamma_d,gamma_dd,gamma_du;
   double F,F_u,F_uu;
   double F_d,F_dd,F_du;
   double n13d;
   double enthrd;
   double G,G_u,G_uu;
   double G_d,G_dd,G_du;
   double fc,fclda,H0,Hu,Hd;
   double ce;
   double Hu_u,Hu_d,Hd_d,Hd_u;
   double fc_u,fc_d;
   double H0_u,H0_d;
   double fclda_u,fclda_d;
   double fc_agr, fc_agrup,fc_agrdn;

   double chiup,chidn;
   double chiup2,chidn2,chiupSQ,chidnSQ;
   double Kup,Kdn,F1up,F1dn,F2up,F2dn;
   double xeup,xedn,xe;
   double fdnxup,fdnxdn,fdagrxup,fdagrxdn;
   double Q,Q1,P,P1,n_m,n2_m,n3_m,n4_m;
   double n13d_m,n13d2_m,n13d3_m,d3;
   double n_thrd,nupndn;
   double n_mthrd,n_mfrthrd,n_mfvthrd;
   double nup_etthrd,nup_fvthrd;
   double ndn_etthrd,ndn_fvthrd;
   double xeup_pbe,xedn_pbe;
   double fdnxup_pbe,fdnxdn_pbe;
   double fdagrxup_pbe,fdagrxdn_pbe;
 
   double pi = 4.00*atan(1.00);
   //double Cf = 3.00*pi*pi;
   //Cf = (3.00*(pow(Cf,(2.00/3.00)))/10.00);

   //double lda_c = (3.00/2.00) * pow((3.00/(4.00*pi)),(1.00/3.00));

   for (int j=0; j<n2ft3d; ++j)
   {
      nup     = dn_in[j]        + 0.50*ETA2;
      ndn     = dn_in[j+n2ft3d] + 0.50*ETA2;
      nup_thrd = pow(nup,thrd);
      ndn_thrd = pow(ndn,thrd);

      agrup   = agr_in[j];
      agrdn   = agr_in[j+n2ft3d];
      agrup2  = agrup*agrup;
      agrdn2  = agrdn*agrdn;

      n       = nup + ndn;

      n2      = n*n;
      nup2    = nup*nup;
      ndn2    = ndn*ndn;

      /* exchange part */
      if ((dn_in[j]+dn_in[j+n2ft3d]) < DNS_CUT)
      {
         xe       = 0.00;
         fdnxup   = 0.00;
         fdnxdn   = 0.00;
         fdagrxup = 0.00;
         fdagrxdn = 0.00;
      }
      else
      {
         /*  UP */
         if (dn_in[j] < DNS_CUT)
         {
            xeup     = 0.00;
            fdnxup   = 0.00;
            fdagrxup = 0.00;
         }
         else
         {
            sdup  = 1.00/(nup_thrd*nup);
            chiup = agrup*sdup;
            chiup2 = chiup*chiup;
            chiupSQ = sqrt(1.00+chiup2);

            Kup = 6.00*beta*log(chiup+chiupSQ);
            F1up = chiup2/(1.00 + chiup*Kup);
            xeup = -nup_thrd*(lda_c + beta*F1up);
            F2up = (2.00 + chiup*Kup 
                 - 6.00*beta*chiup2/chiupSQ)
                 /pow((1.00+chiup*Kup),2.0);
            fdnxup = -nup_thrd*(4.00/3.00)
                  *(lda_c+beta*(F1up-chiup2*F2up));
            fdagrxup = -beta*chiup*F2up ;
         
            if ((fdnxup-xeup) > 0.00)
            {
               gen_PBE96_x_unrestricted(nup,agrup,&xeup_pbe,&fdnxup_pbe,&fdagrxup_pbe);
               Becke_smalln_correction(nup,nup_thrd,1.00,
                                       beta,lda_c,chiup,chiup2,
                                       chiupSQ,Kup,F1up,F2up,
                                       xeup_pbe,fdnxup_pbe,fdagrxup_pbe, 
                                       &xeup,&fdnxup,&fdagrxup);
            }
         } 
         /* END UP */

         /* DOWN */
         if (dn_in[j+n2ft3d] < DNS_CUT)
         {
            xedn     = 0.00;
            fdnxdn   = 0.00;
            fdagrxdn = 0.00;
         }
         else
         {
            sddn  = 1.00/(ndn_thrd*ndn);
            chidn = agrdn*sddn;
            chidn2 = chidn*chidn;
            chidnSQ = sqrt(1.00+chidn2);

            Kdn = 6.00*beta*log(chidn+chidnSQ);
            F1dn = chidn2/(1.00 + chidn*Kdn);
            xedn = -ndn_thrd*(lda_c + beta*F1dn);
            F2dn = (2.00 + chidn*Kdn
                 - 6.00*beta*chidn2/chidnSQ)
                  /pow((1.00+chidn*Kdn),2.00);
            fdnxdn = -ndn_thrd*(4.00/3.00)
                    *(lda_c+beta*(F1dn-chidn2*F2dn));
            fdagrxdn = -beta*chidn*F2dn;
            if ((fdnxdn-xedn) > 0.00)
            {
               gen_PBE96_x_unrestricted(ndn,agrdn,&xedn_pbe,&fdnxdn_pbe,&fdagrxdn_pbe);
               Becke_smalln_correction(ndn,ndn_thrd,1.00,
                                       beta,lda_c,chidn,chidn2,
                                       chidnSQ,Kdn,F1dn,F2dn,
                                       xedn_pbe,fdnxdn_pbe,fdagrxdn_pbe, 
                                       &xedn,&fdnxdn,&fdagrxdn);
            }
         
         }
         /* END DOWN */

         xe = (xeup*nup + xedn*ndn)/n;
      }
      /* end excange part */


      /* correlation part */
      if ((dn_in[j]+dn_in[j+n2ft3d]) < DNS_CUT)
      {
         ce = 0.0;
         fc_u = 0.0;
         fc_d = 0.0;
         fc_agr   = 0.0;
         fc_agrup = 0.0;
         fc_agrdn = 0.0;
      }
      else
      {
         agr     = agr_in[j+2*n2ft3d];
         agr2    = agr*agr;

         n_thrd  = pow(n,thrd);
         n_m     = 1.00/n;
         n2_m    = n_m*n_m;
         n3_m    = n2_m*n_m;
         n4_m    = n3_m*n_m;
         nupndn  = nup*ndn;
         n_mthrd = 1.00/n_thrd;
         n_mfrthrd = n_mthrd*n_m;
         n_mfvthrd = n_mfrthrd*n_mthrd;

         nup_etthrd = pow(nup,etthrd);
         nup_fvthrd = pow(nup,fvthrd);
         ndn_etthrd = pow(ndn,etthrd);
         ndn_fvthrd = pow(ndn,fvthrd);

         gamma   = (4.00*nupndn)*n2_m;
         gamma_u = 4.00*ndn*n2_m - 8.00*nupndn*n3_m;
         gamma_d = 4.00*nup*n2_m - 8.00*nupndn*n3_m;

         gamma_uu = -16.00*ndn*n3_m + 24.00*nupndn*n4_m;
         gamma_dd = -16.00*nup*n3_m + 24.00*nupndn*n4_m;

         gamma_du =  (6.00*gamma-4.00)*n2_m;

         d3         = d/3.00;
         n13d_m     = 1.00/(1.00 + d*n_mthrd);
         n13d2_m    = n13d_m*n13d_m;
         n13d3_m    = n13d2_m*n13d_m;
         F        = gamma*n13d_m;
         F_u      = gamma_u*n13d_m + d3*gamma*n_mfrthrd*n13d2_m;
         F_d      = gamma_d*n13d_m + d3*gamma*n_mfrthrd*n13d2_m;

         F_uu   = gamma_uu*n13d_m + d3*(gamma_u+gamma_u)*n_mfrthrd*n13d2_m
              - (4.00/9.00)*d*gamma*n_mfrthrd*n_m*n13d2_m
              + (2.00/9.00)*d*d*gamma*n_mfrthrd*n_mfrthrd*n13d3_m;

         F_dd   = gamma_dd*n13d_m + d3*(gamma_d+gamma_d)*n_mfrthrd*n13d2_m
              - (4.00/9.00)*d*gamma*n_mfrthrd*n_m*n13d2_m
              + (2.00/9.00)*d*d*gamma*n_mfrthrd*n_mfrthrd*n13d3_m;

         F_du   = gamma_du*n13d_m + d3*(gamma_u+gamma_d)*n_mfrthrd*n13d2_m
              - (4.00/9.00)*d*gamma*n_mfrthrd*n_m*n13d2_m
              + (2.00/9.00)*d*d*gamma*n_mfrthrd*n_mfrthrd*n13d3_m;

         enthrd   = exp(-c*n_mthrd);
         Q        = enthrd*n_mfvthrd;
         Q1 = (Q/3.00)*(c*n_mfrthrd - 5.0*n_m);
         G        = F*Q;

         P  = (1.00/3.00)*c*n_mfrthrd - (5.00/3.00)*n_m;
         P1 = ((-4.00/9.00)*c*n_mfrthrd + (5.0/3.00)*n_m)*n_m;

         G_u      = F_u*Q + G*P;
         G_d      = F_d*Q + G*P;

         G_uu     = F_uu*Q + F_u*Q1 + G_u*P + G*P1;
         G_dd     = F_dd*Q + F_d*Q1 + G_d*P + G*P1;
         G_du     = F_du*Q + F_d*Q1 + G_u*P + G*P1;

         fclda = -a*F*n - 2.00*a*b*G*Cf*pow(2.00,twthrd)
                  *(nup_etthrd + ndn_etthrd);

         fclda_u = -a*F_u*n - a*F 
                   - 2.00*a*b*G_u*Cf*pow(2.00,twthrd)
                   *(nup_etthrd + ndn_etthrd) 
                   - 2.00*a*b*G*Cf*pow(2.00,twthrd)
                   *(8.00/3.00)*nup_fvthrd;

         fclda_d = -a*F_d*n - a*F 
                   - 2.00*a*b*G_d*Cf*pow(2.00,twthrd)
                   *(nup_etthrd + ndn_etthrd) 
                   - 2.00*a*b*G*Cf*pow(2.00,twthrd)
                   *(8.00/3.00)*ndn_fvthrd;

         H0 = (a*b/2.00)*(G 
            + (1.00/3.00)*(nup*G_d + ndn*G_u) 
            + (1.00/4.00)*(nup*G_u + ndn*G_d));

         H0_u=(a*b/2.00)*(G_u
             + (1.00/3.00)*(G_d + nup*G_du + ndn*G_uu)
             + (1.00/4.00)*(G_u + nup*G_uu + ndn*G_du));

         H0_d=(a*b/2.00)*(G_d
             + (1.00/3.00)*(G_u + ndn*G_du + nup*G_dd)
             + (1.00/4.00)*(G_d + ndn*G_dd + nup*G_du));

         Hu = (a*b/18.00)*(G + (15.00/4.00)*nup*G_u 
                              -  (9.00/4.00)*ndn*G_d
                              -        (3.00)*nup*G_d 
                              +  (3.00/2.00)*ndn*G_u);

         Hu_u = (a*b/18.00)*(G_u + (15.00/4.00)*(G_u + nup*G_uu)
                                 -  (9.00/4.00)*(ndn*G_du)
                                 -        (3.00)*(G_d + nup*G_du) 
                                 +  (3.00/2.00)*(ndn*G_uu));

         Hu_d = (a*b/18.00)*(G_d + (15.00/4.00)*(nup*G_du)
                                 -  (9.00/4.00)*(G_d+ndn*G_dd)
                                 -        (3.00)*(nup*G_dd)
                                 +  (3.00/2.00)*(G_u + ndn*G_du));

         Hd = (a*b/18.00)*(G + (15.00/4.00)*ndn*G_d 
                             -  (9.00/4.00)*nup*G_u
                             -        (3.00)*ndn*G_u 
                             +  (3.00/2.00)*nup*G_d);
         Hd_d = (a*b/18.00)*(G_d + (15.00/4.00)*(G_d + ndn*G_dd)
                                 -  (9.00/4.00)*(nup*G_du)
                                 -        (3.00)*(G_u + ndn*G_du) 
                                 +  (3.00/2.00)*(nup*G_dd));
         Hd_u = (a*b/18.00)*(G_u + (15.00/4.00)*(ndn*G_du)
                                 -  (9.00/4.00)*(G_u+nup*G_uu)
                                 -        (3.00)*(ndn*G_uu)
                                 +  (3.00/2.00)*(G_d + nup*G_du));

         fc = fclda + H0*agr2 + Hu*agrup2 + Hd*agrdn2;

         /* calculate derivatives w.r.t up and down density */
         fc_u = fclda_u + H0_u*agr2 + Hu_u*agrup2 + Hd_u*agrdn2;
         fc_d = fclda_d + H0_d*agr2 + Hu_d*agrup2 + Hd_d*agrdn2;

         /* calculate derivatives w.r.t. up,down and total density gradients */
         fc_agr   = 2.00*H0*agr;
         fc_agrup = 2.00*Hu*agrup;
         fc_agrdn = 2.00*Hd*agrdn;

         /* correlation energy dentsity*/
         ce = fc/n;

      }
      /* end correlation part */


      /* return blyp exchange correlation values */
      xce[j]       = x_parameter*xe + c_parameter*ce;
      fn[j]        = x_parameter*fdnxup + c_parameter*fc_u;
      fn[j+n2ft3d] = x_parameter*fdnxdn + c_parameter*fc_d;

      fdn[j]          = x_parameter*fdagrxup + c_parameter*fc_agrup;
      fdn[j+n2ft3d]   = x_parameter*fdagrxdn + c_parameter*fc_agrdn;
      fdn[j+2*n2ft3d] =                        c_parameter*fc_agr;
   }
}

/**************************************************
 *                                                *
 *            gen_BLYP_BW_restricted              *
 *                                                *
 **************************************************/
/* blyp restricted  calc.    

      subroutine gen_BLYP_BW_restricted(n2ft3d,rho_in,agr_in,xce,xcp,fn,fdn)
      input:  n2ft3d                  grid
              rho_in                  density
              agr_in                  absolute gradient of density
               x_parameter:            scale parameter for exchange
              c_parameter:            scale parameter for correlation
      output: xce                     exchange correlation energy density
              fn                      d(n*exc)/dn
              fdn                     d(n*exc)/d(|grad n|)
*/

void gen_BLYP_BW_restricted(const int n2ft3d, 
                            double *rho_in, double *agr_in,
                            const double x_parameter, const double c_parameter,
                            double *xce, double *fn, double *fdn)
{
   /* local declarations */
   double Fc,Gc,C1,C2;
   double n, n_thrd,n_fv,n_fr,n_tw;
   double n_m,n_mthrd,n_mfrthrd,n_mfvthrd;
   double agr,agr2, chi, chi2,chiSQ,sd;
   double K;
   double F1, F2      ;
   double xe,fdnx,fdagrx;
   double xe_pbe,fdnx_pbe,fdagrx_pbe;
   double ce,fdnc, fdagrc,fdnc_lda;
 
   double fc_lda,Ho,Ho_n,Gc_n,Gc_nn,Fc_n,Fc_nn;
   double P,P_n;
      
   double pi       = 4.0*atan(1.0);
   //double two_thrd = pow(2.0,thrd);
    
   for (int i=0; i<n2ft3d; ++i)
   {
      n      = rho_in[i] + ETA;
      agr    = agr_in[i];
      n_thrd = pow(n,thrd);

      if (rho_in[i] < DNS_CUT)
      {
         xe     = 0.00;
         fdnx   = 0.00;
         fdagrx = 0.00;
      }

      /* calc. becke exchange energy density, fnx, fdnx */
      else
      {
         sd    = 1.00/(n_thrd*n); 
         chi   = two_thrd*agr*sd;
         chi2  = chi*chi;
         chiSQ = sqrt(1.00+chi2);

         K  = 6.00*beta*log(chi+chiSQ);
         F1 = chi2/(1.00+chi*K);
         xe = -n_thrd*(lda_c+beta*F1)/two_thrd;
         F2 = (2.00 + chi*K-(chi2)*6.00*beta
                  /chiSQ)
                /((1.00+chi*K)*(1.00+chi*K));
         fdnx = -(n_thrd/two_thrd)*(4.00/3.00)
                 *(lda_c+beta*(F1-chi2*F2));
         fdagrx = -beta*chi*F2 ;
         if ((fdnx-xe) > 0.00)
         {
            gen_PBE96_x_restricted(n,agr,&xe_pbe,&fdnx_pbe,&fdagrx_pbe);

            Becke_smalln_correction(n,n_thrd,2.00,beta,
                                    lda_c,chi,chi2,chiSQ,K,F1,F2,
                                    xe_pbe,fdnx_pbe,fdagrx_pbe,
                                    &xe,&fdnx,&fdagrx);
         }
      }

      /* final result for restricted LYP */
      agr2 = agr*agr;

      n_m     = 1.00/n;
      n_mthrd = 1.00/n_thrd;
      n_mfrthrd = n_mthrd*n_m;
      n_mfvthrd = n_mfrthrd*n_mthrd;

      n_fv = n_thrd*n_thrd*n_thrd*n_thrd*n_thrd;
      n_fr = n_thrd*n_thrd*n_thrd*n_thrd;
      n_tw = n_thrd*n_thrd;

      Fc = (1.00/(1.00+d*n_mthrd));
      Fc_n = thrd_d*n_mfrthrd*Fc*Fc;
      Fc_nn = thrd_d*(-4.00*thrd)*n_mfrthrd*n_m*Fc*Fc
            + thrd_d*n_mfrthrd*2.00*Fc*Fc_n;

      Gc = Fc*exp(-c*n_mthrd)*n_mfvthrd;

      P  = (thrd_d*Fc*n_mfrthrd - 5.00*thrd*n_m + c*thrd*n_mfrthrd);

      P_n = thrd_d*Fc_n*n_mfrthrd
          - thrd_d*Fc  *n_mfrthrd*(4.00*thrd*n_m)
          + 5.00*thrd*n_m*n_m
          - 4.00*thrd*thrd*c*n_mfrthrd*n_m;

      Gc_n  = Gc*P;
      Gc_nn = Gc*P*P + Gc*P_n;

       
      fc_lda   = -a*Fc*n 
               - abCf*Gc*(n_fr*n_fr);
      fdnc_lda = -a*Fc_n*n 
               - a*Fc
               - abCf*Gc_n*n_fr*n_fr
               - 8.00*thrd*abCf*Gc*n_fv;

      Ho   = ho1*Gc   + ho2*Gc_n*n;
      Ho_n = ho3*Gc_n + ho2*Gc_nn*n;

      ce = (fc_lda + Ho*agr2)*n_m;
      fdnc = fdnc_lda + Ho_n*agr2;

      fdagrc = 2.00*Ho*agr;

      xce[i] = x_parameter*xe   + c_parameter*ce;
      fn[i]  = x_parameter*fdnx  + c_parameter*fdnc;
      fdn[i] = x_parameter*fdagrx + c_parameter*fdagrc;
   }

}
