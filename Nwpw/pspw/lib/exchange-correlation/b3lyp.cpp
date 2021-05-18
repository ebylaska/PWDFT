/* b3lyp.cpp -
   Author - Eric Bylaska
*/

#include        <cmath>
#include        "pbe96.hpp"

/* other constants */
#define thrd    (1.00/3.00)
#define twthrd  (2.00/3.00)
#define frthrd  (4.00/3.00)
#define fvthrd  (5.00/3.00)
#define snthrd  (7.00/3.00)
#define etthrd  (8.00/3.00)
#define one6th	(1.00/6.00)
#define two_thrd        1.25992104989487319066


/****************************************************
 *                                                  *
 *             Becke_smalln_correction              *
 *                                                  *
 ****************************************************/
/*
   This routine does small n (larger gradient) correction in which
 the function is slowly switched to PBE when (fdnx-xe)>0.
*/

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


/*---- parameters given by vosko et al -----------------*/
#define ap   3.109070e-02
#define af   1.554530e-02
#define x0p -1.049800e-01
#define x0f -3.250000e-01
#define bp   3.727440e+00
#define bf   7.060420e+00
#define cpt  1.293520e+01
#define cft  1.805780e+01
/*------------------------------------------------------*/

/*     constants calculated from vosko's parameters     */
#define xp   -4.581653e-01
#define xf   -5.772521e-01
#define qp    6.151991e+00
#define qf    4.730927e+00
#define xx0p  1.255491e+01
#define xx0f  1.586879e+01
#define cp1   3.109070e-02
#define cf1   1.554530e-02
#define cp2   9.690228e-04
#define cf2   2.247860e-03
#define cp3   1.049800e-01
#define cf3   3.250000e-01
#define cp4   3.878329e-02
#define cf4   3.878329e-02
#define cp5   3.075995e+00
#define cf5   2.365463e+00
#define cp6   1.863720e+00
#define cf6   3.530210e+00
#define dp1   6.218140e-02
#define df1   3.109060e-02
#define dp2   1.938045e-03
#define df2   4.495720e-03
#define dp3   1.049800e-01
#define df3   3.250000e-01
#define dp4  -3.205972e-02
#define df4  -1.779316e-02
#define dp5  -1.192972e-01
#define df5  -1.241661e-01
#define dp6   1.863720e+00
#define df6   3.530210e+00
#define dp7   9.461748e+00
#define df7   5.595417e+00
//#define fc    1.923661e+00
//#define fd    2.564881e+00
#define crs   7.876233e-01

/* tolerance constants */
#define ETA     1.0e-20
#define ETA2    2.0e-20
#define DNS_CUT 1.0e-18
#define tolrho  2.0e-11
#define minagr  1.0e-12

/* LYP parameters */
#define a       0.04918
#define b       0.132
#define c       0.2533
#define d       0.349
#define Cf      2.87123400018819108225

/* collated LYP parameters */
#define ho1     ((19.00/36.00)*a*b)
#define ho2     (( 7.00/24.00)*a*b)
#define ho3     ((59.00/72.00)*a*b)
#define thrd_d  (thrd*d)
#define abCf    (a*b*Cf)

/* Becke88 parameters */
#define BETA    0.0042 
#define LDA_C   0.93052573634910018540



/**************************************************
 *                                                *
 *           gen_B3LYP_BW_unrestricted            *
 *                                                *
 **************************************************/

void gen_B3LYP_BW_unrestricted(const int n2ft3d,
                               double  *dn_in, double *agr_in,
                               const double x_parameter, const double c_parameter,
                               double *xce, double *fn, double *fdn)
{

   double nup,ndn,n,agrup,agrdn,agr;
   double agr2, agrup2,agrdn2;
   double n2,nup2,ndn2;
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
   double xeup_pbe,fdnxup_pbe,fdagrxup_pbe;
   double xedn_pbe,fdnxdn_pbe,fdagrxdn_pbe;
   double xeup,xedn,xe;
   double fdnxup,fdnxdn,fdagrxup,fdagrxdn;
   double  Q,Q1,P,P1,n_m,n2_m,n3_m,n4_m;
   double  n13d_m,n13d2_m,n13d3_m,d3;
   double  n_thrd,nupndn;
   double  n_mthrd,n_mfrthrd,n_mfvthrd;
   double   nup_etthrd,nup_fvthrd;
   double   ndn_etthrd,ndn_fvthrd;
   double x,xi,ff,dff,xxp,dxxp,xxf,dxxf;
   double ex_p,ux_p,ex_f,ux_f;
   double ec_p,uc_p,ec_f,uc_f;
   double xe_dirac,fnxup_dirac,fnxdn_dirac;
   double ce_vosko,fncup_vosko,fncdn_vosko;

   double pi = 4.0*atan(1.0);
   double AA = pow(2.00,thrd)-1.0;

   for (int j=0; j<n2ft3d; ++j)
   {
      nup     = dn_in[j]        + 0.50*ETA2;
      ndn     = dn_in[j+n2ft3d] + 0.50*ETA2;

      agrup   = agr_in[j];
      agrdn   = agr_in[j+n2ft3d];
      agrup2  = agrup*agrup;
      agrdn2  = agrdn*agrdn;

      n       = nup + ndn;
      agr     = agr_in[j+2*n2ft3d];
      agr2    = agr*agr;

      n2      = n*n;
      nup2    = nup*nup;
      ndn2    = ndn*ndn;

      /* LSDA terms */
      x = pow(crs/n,one6th);
      xi = (nup-ndn)/n;

      ff =  ( pow((1.00+xi),frthrd) + pow((1.00-xi),frthrd) - 2.00)/(2.00*AA);
      dff = frthrd*(  pow((1.00+xi),thrd) - pow((1.00-xi),thrd)) /(2.00*AA);

      /* dirac exchange */
      ex_p   = (xp/pow(x,2));
      ux_p   = frthrd*(xp/pow(x,2));
      ex_f   = (xf/pow(x,2));
      ux_f   = frthrd*(xf/pow(x,2));
      xe_dirac    = ex_p + ff*(ex_f-ex_p);
      fnxup_dirac = ux_p + ff*(ux_f-ux_p) + (+1.00-xi)*dff*(ex_f-ex_p);
      fnxdn_dirac = ux_p + ff*(ux_f-ux_p) + (-1.00-xi)*dff*(ex_f-ex_p);

      /* vosko correlation */
      xxp  = pow(x,2) + bp*x + cpt;
      dxxp = 2.00*x + bp;
      ec_p = cp1*log(pow(x,2)/xxp) + cp2*log( (x+cp3)*(x+cp3)/xxp) + cp4*atan(cp5/(x+cp6));
      uc_p = ec_p - one6th*x*(dp1/x + dp2/(x+cp3) + dp4*dxxp/xxp + dp5/( (x+dp6)*(x+dp6)+dp7));
      
      xxf  = pow(x,2) + bf*x + cft;
      dxxf = 2.00*x + bf;
      ec_f = cf1*log(pow(x,2)/xxf) + cf2*log( (x+cf3)*(x+cf3)/xxf) + cf4*atan(cf5/(x+cf6));
      uc_f = ec_f - one6th*x*(  df1/x + df2/(x+cf3) + df4*dxxf/xxf + df5/( (x+df6)*(x+df6)+df7));
      ce_vosko    = ec_p + ff*(ec_f - ec_p);
      fncup_vosko = uc_p + ff*(uc_f-uc_p) + (+1.00-xi)*dff*(ec_f-ec_p);
      fncdn_vosko = uc_p + ff*(uc_f-uc_p) + (-1.00-xi)*dff*(ec_f-ec_p);

      /***********************************/
      /**********exchange part************/
      /***********************************/
      if ((dn_in[j]+dn_in[j+n2ft3d])<DNS_CUT)
      {
         xe       = 0.00;
         fdnxup   = 0.00;
         fdnxdn   = 0.00;
         fdagrxup = 0.00;
         fdagrxdn = 0.00;
      }
      else
      {
         /**************UP*******************/
         if (dn_in[j]<DNS_CUT)
         {
             xeup     = 0.00;
             fdnxup   = 0.00;
             fdagrxup = 0.00;
         }
         else
         {
            chiup = pow(agrup/nup,(4.00/3.00));
            chiup2 = chiup*chiup;
            chiupSQ = sqrt(1.00+chiup2);

            Kup = 6.00*BETA*log(chiup+chiupSQ);
            F1up = chiup2/(1.00 + chiup*Kup);
            xeup = -pow(nup,thrd)*(LDA_C + BETA*F1up);
            F2up = (2.00 + chiup*Kup - 6.00*BETA*chiup2/chiupSQ)/pow((1.00+chiup*Kup),2.00);
            fdnxup = -pow(nup,(thrd))*(4.00/3.00)*(LDA_C+BETA*(F1up-chiup2*F2up));
            fdagrxup = -BETA*chiup*F2up;
            if ((fdnxup-xeup)> 0.00)
            {
               gen_PBE96_x_unrestricted(nup,agrup,&xeup_pbe,&fdnxup_pbe,&fdagrxup_pbe);

               Becke_smalln_correction(nup,pow(nup,thrd),1.0,
                                            BETA,LDA_C,chiup,chiup2,
                                            chiupSQ,Kup,F1up,F2up,
                                      xeup_pbe,fdnxup_pbe,fdagrxup_pbe, 
                                      &xeup,&fdnxup,&fdagrxup);
            }
         } 
         /*************END UP*****************/

         /*************DOWN******************/
         if (dn_in[j+n2ft3d]<DNS_CUT)
         {
            xedn     = 0.00;
            fdnxdn   = 0.00;
            fdagrxdn = 0.00;
         }
         else
         {
            chidn = agrdn/pow(ndn,(4.00/3.00));
            chidn2 = chidn*chidn;
            chidnSQ = sqrt(1.00+chidn2);

            Kdn = 6.00*BETA*log(chidn+chidnSQ);
            F1dn = chidn2/(1.00 + chidn*Kdn);
            xedn = -pow(ndn,thrd)*(LDA_C + BETA*F1dn);
            F2dn = (2.00 + chidn*Kdn - 6.00*BETA*chidn2/chidnSQ)/pow((1.00+chidn*Kdn),2.00);
            fdnxdn = -pow(ndn,(thrd))*(4.00/3.00)*(LDA_C+BETA*(F1dn-chidn2*F2dn));
            fdagrxdn = -BETA*chidn*F2dn;
            if ((fdnxdn-xedn)>0.00) 
            {
               gen_PBE96_x_unrestricted(ndn,agrdn,&xedn_pbe,&fdnxdn_pbe,&fdagrxdn_pbe);
               Becke_smalln_correction(ndn,pow(ndn,thrd),1.00,
                                            BETA,LDA_C,chidn,chidn2,
                                            chidnSQ,Kdn,F1dn,F2dn,
                                      xedn_pbe,fdnxdn_pbe,fdagrxdn_pbe, 
                                      &xedn,&fdnxdn,&fdagrxdn);
            }
         }
         /*************END DOWN******************/

          xe = (xeup*nup + xedn*ndn)/n;

      } 
      /***********************************/
      /*********end excange part**********/
      /***********************************/
      

      /***********************************/
      /*********correlation part**********/
      /***********************************/
      n_thrd   = pow(n,thrd);
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

      fclda = -a*F*n - 2.00*a*b*G*Cf*(pow(2.00,twthrd))*(nup_etthrd + ndn_etthrd);

      fclda_u = -a*F_u*n - a*F 
                - 2.00*a*b*G_u*Cf*(pow(2.00,twthrd))
                *(nup_etthrd + ndn_etthrd) 
                - 2.00*a*b*G*Cf*(pow(2.00,twthrd))
                *(8.00/3.00)*nup_fvthrd;

      fclda_d = -a*F_d*n - a*F 
                - 2.00*a*b*G_d*Cf*(pow(2.00,twthrd))
                *(nup_etthrd + ndn_etthrd) 
                - 2.00*a*b*G*Cf*(pow(2.00,twthrd))
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
                          -       (3.00)*nup*G_d 
                          +  (3.00/2.00)*ndn*G_u);

      Hu_u = (a*b/18.00)*(G_u + (15.00/4.00)*(G_u + nup*G_uu)
                              -  (9.00/4.00)*(ndn*G_du)
                              -       (3.00)*(G_d + nup*G_du) 
                              +  (3.00/2.00)*(ndn*G_uu));

      Hu_d = (a*b/18.00)*(G_d + (15.00/4.00)*(nup*G_du)
                              -  (9.00/4.00)*(G_d+ndn*G_dd)
                              -       (3.00)*(nup*G_dd)
                              +  (3.00/2.00)*(G_u + ndn*G_du));

      Hd = (a*b/18.00)*(G + (15.00/4.00)*ndn*G_d 
                          -  (9.00/4.00)*nup*G_u
                          -       (3.00)*ndn*G_u 
                          +  (3.00/2.00)*nup*G_d);

      Hd_d = (a*b/18.00)*(G_d + (15.00/4.00)*(G_d + ndn*G_dd)
                              -  (9.00/4.00)*(nup*G_du)
                              -       (3.00)*(G_u + ndn*G_du) 
                              +  (3.00/2.00)*(nup*G_dd));

      Hd_u = (a*b/18.00)*(G_u + (15.00/4.00)*(ndn*G_du)
                              -  (9.00/4.00)*(G_u+nup*G_uu)
                              -       (3.00)*(ndn*G_uu)
                              +  (3.00/2.00)*(G_d + nup*G_du));


      fc = fclda + H0*agr2 + Hu*agrup2 + Hd*agrdn2;

      /* calculate derivatives w.r.t up and down density */
      fc_u = fclda_u + H0_u*agr2 + Hu_u*agrup2 + Hd_u*agrdn2;
      fc_d = fclda_d + H0_d*agr2 + Hu_d*agrup2 + Hd_d*agrdn2;

      /* calculate derivatives w.r.t. up,down and total density gradients */
      fc_agr   = 2.00*H0*agr;
      fc_agrup = 2.00*Hu*agrup;
      fc_agrdn = 2.00*Hd*agrdn;

      ce = fc/n;

      /***********************************/
      /*******end correlation part********/
      /***********************************/


      /**** a*HF + (1-a)*(0.1*Dirac + 0.9*Becke88) + 0.81*LYP + 0.19*Vosko ****/
      /**** x_parameter = (1-a)                                            ****/

      /* return blyp exchange correlation values */
      xce[j]       = x_parameter*0.90*xe     + c_parameter*0.810*ce;
      fn[j]        = x_parameter*0.90*fdnxup + c_parameter*0.810*fc_u;
      fn[j+n2ft3d] = x_parameter*0.90*fdnxdn + c_parameter*0.810*fc_d;

      fdn[j]          = x_parameter*0.90*fdagrxup + c_parameter*0.81*fc_agrup;
      fdn[j+n2ft3d]   = x_parameter*0.90*fdagrxdn + c_parameter*0.81*fc_agrdn;
      fdn[j+2*n2ft3d] =                             c_parameter*0.81*fc_agr;

      xce[j]       = xce[j]        + x_parameter*0.100*xe_dirac    + c_parameter*0.190*ce_vosko;
      fn[j]        = fn[j]         + x_parameter*0.100*fnxup_dirac + c_parameter*0.190*fncup_vosko;
      fn[j+n2ft3d] = fn[j+n2ft3d]  + x_parameter*0.100*fnxdn_dirac + c_parameter*0.190*fncdn_vosko;
   }
}

/**************************************************
 *                                                *
 *            gen_B3LYP_BW_restricted             *
 *                                                *
 **************************************************/
/*
*    blyp restricted  calc.    
*      subroutine gen_B3LYP_BW_restricted(n2ft3d,rho_in,agr_in,xce,xcp,fn,fdn)
*      input:  n2ft3d                  grid
*              rho_in                  density
*              agr_in                  absolute gradient of density
*              x_parameter:            scale parameter for exchange
*              c_parameter:            scale parameter for correlation
*      output: xce                     exchange correlation energy density
*              fn                      d(n*exc)/dn
*              fdn                     d(n*exc)/d(|grad n|)
*/

void gen_B3LYP_BW_restricted(const int n2ft3d,double *rho_in,double *agr_in,
                             const double x_parameter, const double c_parameter,
                             double *xce, double *fn, double *fdn)
{
   /* local variables */
   double Fc,Gc,C1,C2;
   double n, n_thrd,n_fv,n_fr,n_tw;
   double n_m,n_mthrd,n_mfrthrd,n_mfvthrd;
   double agr,agr2, chi, chi2,chiSQ,sd;
   double K;
   double F1, F2;
   double xe_pbe,fdnx_pbe,fdagrx_pbe;
   double xe,fdnx,fdagrx,xe_dirac,fnx_dirac,x,xx;
   double ce,fdnc, fdagrc,fdnc_lda,ce_vosko,fnc_vosko;

   double fc_lda,Ho,Ho_n,Gc_n,Gc_nn,Fc_n,Fc_nn;
   double P,P_n;
      
   double pi = 4.0*atan(1.0);

   double p2 = -5.00*thrd;
   double p3 =  c*thrd;
   double p4 = -4.00*thrd*thrd_d;
   double p5 =  5.00*thrd;
   double p6 = -4.00*thrd*thrd*c;
      

   for (int i=0; i<n2ft3d; ++i)
   {
      n      = rho_in[i] + ETA;
      agr    = agr_in[i];
      n_thrd = pow(n,thrd);

      /* LDA terms */
      x = crs/pow(n,one6th);

      /* dirac exchange */
      xe_dirac   = (xp/pow(x,2));
      fnx_dirac = frthrd*(xp/pow(x,2));

      /* vosko correlation */
      xx=1.00/(x*(x+bp)+cpt);
      ce_vosko=cp1*log(xx*pow(x,2))+cp2*log(xx*pow((x+cp3),2)) + cp4*atan(cp5/(x+cp6));
      fnc_vosko = ce_vosko-one6th*x*(dp1/x+dp2/(x+dp3)+dp4*xx*(2.00*x+bp)+dp5/(pow((x+dp6),2)+dp7));

      /***********************************************************/
      /******calc. becke exchange energy density, fnx, fdnx*******/
      /***********************************************************/
      if (rho_in[i]<DNS_CUT)
      {
         xe     = 0.00;
         fdnx   = 0.00;
         fdagrx = 0.00;
      }
      else
      {
         sd       = 1.00/(n_thrd*n);
         chi      = two_thrd*agr*sd;
         chi2     = chi*chi;
         chiSQ    = sqrt(1.00+chi2);

         K        = 6.00*BETA*log(chi+chiSQ);
         F1       = chi2/(1.00+chi*K);
         xe       = -n_thrd*(LDA_C+BETA*F1)/two_thrd;
         F2       = (2.00 + chi*K-(chi2)*6.00*BETA/chiSQ)/((1.00+chi*K)*(1.00+chi*K));
         fdnx     = -(n_thrd/two_thrd)* (4.00/3.00) *(LDA_C+BETA*(F1-chi2*F2));
         fdagrx   = -BETA*chi*F2;
         if ((fdnx-xe)>0.00)
         {
            gen_PBE96_x_restricted(n,agr,&xe_pbe,&fdnx_pbe,&fdagrx_pbe);
            Becke_smalln_correction(n,pow(n,thrd),2.00,BETA,
                                         LDA_C,chi,chi2,chiSQ,K,F1,F2,
                                         xe_pbe,fdnx_pbe,fdagrx_pbe,
                                         &xe,&fdnx,&fdagrx);
         }
      }


      /*******final result for restricted LYP*********/
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
      Fc_nn = thrd_d*(-4.00*thrd)*n_mfrthrd*n_m*Fc*Fc + thrd_d*n_mfrthrd*2.00*Fc*Fc_n;

      Gc = Fc*exp(-c*n_mthrd)*n_mfvthrd;

      P  = (thrd_d*Fc*n_mfrthrd - 5.00*thrd*n_m + c*thrd*n_mfrthrd);

      P_n = thrd_d*Fc_n*n_mfrthrd 
          - thrd_d*Fc  *n_mfrthrd*(4.00*thrd*n_m)
          + 5.00*thrd*n_m*n_m
          - 4.00*thrd*thrd*c*n_mfrthrd*n_m;

      Gc_n  = Gc*P;
      Gc_nn = Gc*P*P + Gc*P_n;

       
      fc_lda   = -a*Fc*n - abCf*Gc*(n_fr*n_fr);
      fdnc_lda = -a*Fc_n*n 
               - a*Fc
               - abCf*Gc_n*n_fr*n_fr
               - 8.00*thrd*abCf*Gc*n_fv;

      Ho   = ho1*Gc   + ho2*Gc_n*n;
      Ho_n = ho3*Gc_n + ho2*Gc_nn*n;

      ce = (fc_lda + Ho*agr2)*n_m;
      fdnc = fdnc_lda + Ho_n*agr2;

      fdagrc = 2.00*Ho*agr;


      /*** a*HF + (1-a)*(0.1*Dirac + 0.9*Becke88) + 0.81*LYP + 0.19*Vosko ***/
      /*** x_parameter = (1-a)                                            ***/
 
      xce[i] = x_parameter*0.90*xe     + c_parameter*0.810*ce;
      fn[i]  = x_parameter*0.90*fdnx   + c_parameter*0.810*fdnc;
      fdn[i] = x_parameter*0.90*fdagrx + c_parameter*0.810*fdagrc;

      xce[i] = xce[i] + x_parameter*0.100*xe_dirac  + c_parameter*0.190*ce_vosko;
      fn[i]  = fn[i]  + x_parameter*0.100*fnx_dirac + c_parameter*0.190*fnc_vosko;
   }
}
