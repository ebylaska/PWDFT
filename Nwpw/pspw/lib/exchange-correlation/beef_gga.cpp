/* beef_gga.cpp
  Author - Eric Bylaska
*/

#include        <cmath>

namespace pwdft {



/* Density cutoff parameters */
#define DNS_CUT		(1.0e-20)
#define ETA		(1.0e-20)
#define ETA2		(1.0e-14)
#define alpha_zeta	(1.0e0-ETA2)
#define alpha_zeta2	(1.ee0-ETA2)

/* PBE96 GGA exchange constants */
#define MU		(0.2195149727645171e0)
#define KAPPA 		(0.8040000000000000e0)

/* PBE96 GGA correlation constants */
#define GAMMA	(0.031090690869655e0)
#define BETA	(0.066724550603149e0)
#define BOG	(BETA/GAMMA)

/* Perdew-Wang92 LDA correlation coefficients */
#define GAM	(0.519842099789746329e0)
#define iGAM	(1.0e0/GAM)
#define FZZ	(8.0e0/(9.0e0*GAM))
#define iFZZ	(0.125e0*9.0e0*GAM)

#define A_1	(0.0310907e0)
#define A1_1	(0.2137000e0)
#define B1_1	(7.5957000e0)
#define B2_1	(3.5876000e0)
#define B3_1	(1.6382000e0)
#define B4_1	(0.4929400e0)

#define A_2	(0.01554535e0)
#define A1_2	(0.20548000e0)
#define B1_2	(14.11890000e0)
#define B2_2	(6.19770000e0)
#define B3_2	(3.36620000e0)
#define B4_2	(0.62517000e0)

#define A_3	(0.0168869e0)
#define A1_3	(0.1112500e0)
#define B1_3	(10.3570000e0)
#define B2_3	(3.6231000e0)
#define B3_3	(0.8802600e0)
#define B4_3	(0.4967100e0)

/* Perdew-Wang92 LDA functional */
static void LSDT(double a, double a1, double b1, double b2, double b3, double b4, double srs, double *ec, double *ec_rs)
{
   double q0,q1,q1p,qd,ql;
   q0  = -2.00*a*(1.0e0+a1*srs*srs);
   q1  = 2.00*a*srs*(b1+srs*(b2+srs*(b3+srs*b4)));
   q1p = a*((b1/srs)+2.00*b2+srs*(3.00*b3+srs*4.00*b4));
   qd  =1.00/(q1*q1+q1);
   ql  = -log(qd*q1*q1);

   *ec = q0*ql;
   *ec_rs = -2.0e0*a*a1*ql-q0*q1p*qd;
}

/* other constants */
#define onethird        (1.00/3.00)
#define onethirdm       (-1.00/3.00)
#define twothird        (2.00/3.00)
#define fourthird       (4.00/3.00)
#define fivethird       (5.00/3.00)
#define onesixthm       (-1.00/6.00)
#define sevensixthm     (-7.00/6.00)
#define sevensixths     (7.00/6.00)


/* BEEF expansion parameters given by Wellendorff et al */
static double am[30] = { 1.516501714e0,   4.413532099e-1,-9.182135241e-2,
                        -2.352754331e-2,  3.418828455e-2, 2.411870076e-3,
                        -1.416381352e-2,  6.975895581e-4, 9.859205137e-3,
                        -6.737855051e-3, -1.573330824e-3, 5.036146253e-3,
                        -2.569472453e-3, -9.874953976e-4, 2.033722895e-3,
                        -8.018718848e-4, -6.688078723e-4, 1.030936331e-3,
                        -3.673838660e-4, -4.213635394e-4, 5.761607992e-4,
                        -8.346503735e-5, -4.458447585e-4, 4.601290092e-4,
                        -5.231775398e-6, -4.239570471e-4, 3.750190679e-4,
                         2.114938125e-5, -1.904911565e-4, 7.384362421e-5 };


/****************************************
 *					*
 *       gen_BEEF_BW_unrestricted	*
 *					*
 ****************************************/
/*
     This function returns the BEEF exchange-correlation
   energy density, xce, and its derivatives with respect
   to nup, ndn, |grad nup|, |grad ndn|, and |grad n|.
 
    Entry - n2ft3d     : number of grid points
            dn_in(*,2) : spin densites nup and ndn
            agr_in(*,3): |grad nup|, |grad ndn|, and |grad n|
            x_parameter: scale parameter for exchange
            c_parameter: scale parameter for correlation
            alphac: scale parameter for BEEF
 
    Exit - xce(*)  : PBE96 energy density
         - fn(*,2) : d(n*xce)/dnup, d(n*xce)/dndn
         - fdn(*,3): d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|
                     d(n*xce)/d|grad n|
*/
void gen_BEEF_BW_unrestricted(const int n2ft3d,
                              double *dn_in, double *agr_in,
                              const double x_parameter, const double c_parameter, double alphac,
                              double *xce, double *fn, double *fdn)
{

   /* local variables */
   double n,agr,nup,agrup,ndn,agrdn;
   double kf,ks,s,P0,n_onethird;
   double rs;        // Wigner radius
   double rss;       // rss  = sqrt(rs)
   double rs_n;      // rs_n = n*drs/dn
   double t,t2,t4,t6;
   double t_nup;     // t_nup = n*dt/dnup
   double t_ndn;     // t_ndn = n*dt/dndn
   double t_agr;     // t_agr = n*dt/dagr
   double zet,twoksg;
   double zet_nup;   // zet_nup = n*dzet/dnup
   double zet_ndn;   // zet_nup = n*dzet/dnup
   double zetp_1_3,zetm_1_3;
   double zetpm_1_3,zetmm_1_3;
   double phi,phi3,phi4;
   double phi_zet;
   double A,A2;
   double A_phi,A_ec_lda;
   double Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8;
   double PON,FZ,z4;
   double tau;
   double F;
   double Fs;                // dF/ds
   double Hpbe;
   double Hpbe_t;            // dHpbe/dt
   double Hpbe_phi;          // dHpbe/dphi
   double Hpbe_ec_lda;       // dHpbe/d(ec_lda)
   double Hpbe_nup,Hpbe_ndn; // n*dHpbe/dnup, n*dHpbe/dndn
   double Ipbe;
   double Ipbe_t,Ipbe_A;     // dIpbe/dt, dIpbe/dA

   double exup,exdn,ex,ex_lda;
   double ecu,ecp,eca,ec,ec_lda;
   double ecu_rs,ecp_rs,eca_rs;
   double ec_lda_rs,ec_lda_zet;  // d(ec_lda)/drs, d(ec_lda)/dzet
   double ec_lda_nup,ec_lda_ndn; // n*d(ec_lda)/dnup, n*d(ec_lda)/dndn
   double fnxup,fdnxup;          // d(n*ex)/dnup, d(n*ex)/dndn
   double fnxdn,fdnxdn;          // d(n*ex)/d|grad nup|, d(n*ex)/d|grad ndn|
   double fncup,fncdn;           // d(n*ec)/dnup, d(n*ec)/dndn

   double p0,p1,p2,dp,s2,oneovers2,Ft,dp2,sgn;
   double malphac=1.00-alphac;

   double pi = 4.00*atan(1.00);
   double rs_scale = pow((0.750/pi),onethird);
   double fdnx_const = -3.00/(8.00*pi);



   for (int i=0; i<n2ft3d; ++i)
   {
      nup     = dn_in[i]+ETA;
      agrup   = agr_in[i];
 
      ndn     = dn_in[i+n2ft3d]+ETA;
      agrdn   = agr_in[i+n2ft3d];
 
     /****************************************************************
      ***** calculate polarized Exchange energies and potentials *****
      ****************************************************************/

      /**** up ****/
      n     = 2.00*nup;
      agr   = 2.00*agrup;

      n_onethird = pow((3.00*n/pi),onethird);
      ex_lda     = -0.750*n_onethird;

      kf = pow((3.00*pi*pi*n),onethird);
      s  = agr/(2.00*kf*n);

      t  = 2.00*s*s/(4.00+s*s) - 1.00;
      if (std::abs(t)>1.00) t = (t>0.0) ? 1.0 : -1.0;
      F  = 0.00;
      Ft = 0.00;
         
      s2   = t*t - 1.00;
      if (std::abs(s2)<1.0e-12) 
      {
         if (t>0.00) 
         {
            for (int j=0; j<30; ++j)
            {
               F  = F  + am[j];
               Ft = Ft + am[j]*0.50*((double) (j*(j+1)));
            }
         }
         else
         {
            sgn = 1.00;
            for (int j=0; j<30; ++j)
            {
               F  = F  + sgn*am[j];
               Ft = Ft + sgn*am[j]*0.50*((double) (j*(j+1)));
               sgn = -sgn;
            }
         }
      }
      else
      {
         oneovers2 = 1.00/s2;
         p0 = 1.00;
         dp = 0.0;
         F  = F  + am[0]*p0;
         Ft = Ft + am[0]*dp;
         p1 = t;
         dp = 1.00;
         F  = F  + am[1]*t;
         Ft = Ft + am[1]*dp;
         for (int j=2; j<30; ++j)
         {
            p2    = (1.00/((double) (j+1))) * ((2*j+1)*t*p1 - ((double) j)*p0);
            dp2   = ((double) j)*oneovers2*(t*p2-p1);
            F  = F  + am[j]*p2;
            Ft = Ft + am[j]*dp2;
            p0 = p1;
            p1 = p2;
         }
      }
      Fs  = (16.00*s/pow((4.00+s*s),2)) * Ft;
        
      exup = ex_lda*F;
      fnxup = fourthird*(exup - ex_lda*Fs*s);
      fdnxup = fdnx_const*Fs;


     /**** down ****/
      n     = 2.00*ndn;
      agr   = 2.00*agrdn;

      n_onethird = pow((3.00*n/pi),onethird);
      ex_lda     = -0.750*n_onethird;

      kf = pow((3.00*pi*pi*n),onethird);
      s  = agr/(2.00*kf*n);

      t  = 2.00*s*s/(4.00+s*s) - 1.00;
      if (std::abs(t)>1.00) t = (t>0.0) ? 1.0 : -1.0;
      F  = 0.00;
      Ft = 0.00;
         
      s2   = t*t - 1.00;
      if (std::abs(s2)<1.0e-12)
      {
         if (t>0.00)
         {
            for (int j=0; j<30; ++j)
            {
               F  = F  + am[j];
               Ft = Ft + am[j]*0.50* ((double) (j*(j+1)));
            }
         }
         else
         {
            sgn = 1.0;
            for (int j=0; j<30; ++j)
            {
               F  = F  + sgn*am[j];
               Ft = Ft + sgn*am[j]*0.50* ((double) (j*(j+1)));
               sgn = -sgn;
            }
         }
      }
      else
      {
         oneovers2 = 1.00/s2;
         p0 = 1.00;
         dp = 0.0;
         F  = F  + am[0]*p0;
         Ft = Ft + am[0]*dp;
         p1 = t;
         dp = 1.00;
         F  = F  + am[1]*t;
         Ft = Ft + am[1]*dp;
         for (int j=2; j<30; ++j)
         {
            p2    = (1.00/((double)(j+1))) * ((2*j+1)*t*p1 - ((double) j)*p0);
            dp2   = ((double) j)*oneovers2*(t*p2-p1);
            F  = F  + am[j]*p2;
            Ft = Ft + am[j]*dp2;
            p0 = p1;
            p1 = p2;
         }
      }
      Fs  = (16.00*s/pow((4.00+s*s),2)) * Ft;

      exdn   = ex_lda*F;
      fnxdn  = fourthird*(exdn - ex_lda*Fs*s);
      fdnxdn = fdnx_const*Fs;

      n = nup+ndn;

      ex = (exup*nup+ exdn*ndn)/ n;

     /*******************************************************************
      ***** calculate polarized correlation energies and potentials *****
      *******************************************************************/
      agr   = agr_in[i+2*n2ft3d];

      zet = (nup-ndn)/n;
      zet_nup = -(zet - 1.00);
      zet_ndn = -(zet + 1.00);
      zetpm_1_3 = pow((1.00+zet*alpha_zeta),onethirdm);
      zetmm_1_3 = pow((1.00-zet*alpha_zeta),onethirdm);
      zetp_1_3  = pow((1.00+zet*alpha_zeta)*zetpm_1_3,2);;
      zetm_1_3  = pow((1.00-zet*alpha_zeta)*zetmm_1_3,2);

      zetp_1_3  = (1.00+zet*alpha_zeta)*zetpm_1_3; zetp_1_3 *= zetp_1_3;
      zetm_1_3  = (1.00-zet*alpha_zeta)*zetmm_1_3; zetm_1_3 *= zetm_1_3;

      phi = 0.50*(zetp_1_3*zetp_1_3 + zetm_1_3*zetm_1_3);
      phi_zet = alpha_zeta*( zetpm_1_3 - zetmm_1_3)/3.00;
      F  =( (1.00+zet*alpha_zeta)*zetp_1_3
          + (1.00-zet*alpha_zeta)*zetm_1_3 - 2.00)*iGAM;
      FZ = (zetp_1_3 - zetm_1_3)*(alpha_zeta*fourthird*iGAM);

      /* calculate Wigner radius */
      rs  = rs_scale/pow(n,onethird);
      rss = sqrt(rs);

      /* calculate n*drs/dn */
      rs_n = onethirdm*rs/n;
      rs_n = onethirdm*rs;

      /* calculate t */
      kf = pow((3.00*pi*pi*n),onethird);
      ks = sqrt(4.00*kf/pi);

      twoksg = 2.00*ks*phi;

      t  = agr/(twoksg*n);

      /* calculate n*dt/dnup, n*dt/dndn, n*dt/d|grad n| */
      t_nup = sevensixthm*t - (phi_zet)*(zet_nup)*t/phi;
      t_ndn = sevensixthm*t - (phi_zet)*(zet_ndn)*t/phi;
      t_agr  = 1.00/(twoksg);

     /**************************************************
      ***** compute LSDA correlation energy density ****
      **************************************************/
      LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,&ecu,&ecu_rs);
      LSDT(A_2,A1_2,B1_2,B2_2,B3_2,B4_2,rss,&ecp,&ecp_rs);
      LSDT(A_3,A1_3,B1_3,B2_3,B3_3,B4_3,rss,&eca,&eca_rs);

      z4 = zet*zet*zet*zet;

      ec_lda = ecu*(1.00-F*z4)
             + ecp*F*z4
             - eca*F*(1.00-z4)/FZZ;

      ec_lda_rs = ecu_rs*(1.00-F*z4)
                + ecp_rs*F*z4
                - eca_rs*F*(1.00-z4)/FZZ;

      ec_lda_zet = (4.00*(zet*zet*zet)*F + FZ*z4)*(ecp-ecu+eca*iFZZ)
                 - FZ*eca*iFZZ;
                  
      /********************************************
       **** calculate PBE96 correlation energy ****
       ********************************************/
      phi3 = phi*phi*phi;
      phi4 = phi3*phi;
      PON  = -ec_lda/(phi3*GAMMA);
      tau  = exp(PON);

      A = BOG/(tau-1.00+ETA);
      A2 = A*A;
      t2 = t*t;
      t4 = t2*t2;
      t6 = t4*t2;
      Q4 = 1.00 + A*t2;
      Q5 = 1.00 + 2.00*A*t2;
      Q6 = 2.00 + A*t2;
      Q7 = 1.00+A*t2+A2*t4;
      Q8 = Q7*Q7;

      Ipbe = 1.00 + BOG*t2*Q4/Q7;
      Hpbe = GAMMA*phi3*log(Ipbe);

      Ipbe_t =  BOG*(2.00*t)*Q5/Q8;
      Ipbe_A = -BOG*(A*t6)   *Q6/Q8;

      A_ec_lda  = tau/(BETA*phi3)*A2;
      A_phi     = -3.00*ec_lda*tau/(BETA*phi4)*A2;

      Hpbe_ec_lda = (GAMMA*phi3/Ipbe)*Ipbe_A*A_ec_lda;

      Hpbe_phi    = 3.00*Hpbe/phi
                  + (GAMMA*phi3/Ipbe)*Ipbe_A*A_phi;

      Hpbe_t      = (GAMMA*phi3/Ipbe)*Ipbe_t;

      ec_lda_nup = ec_lda_zet
                 - zet * ec_lda_zet
                 + rs_n * ec_lda_rs;
      ec_lda_ndn = -ec_lda_zet
                 - zet  * ec_lda_zet
                 + rs_n * ec_lda_rs;


      Hpbe_nup  = ec_lda_nup   * Hpbe_ec_lda
             + phi_zet*zet_nup * Hpbe_phi
             + t_nup           * Hpbe_t;

      Hpbe_ndn  = ec_lda_ndn   * Hpbe_ec_lda
             + phi_zet*zet_ndn * Hpbe_phi
             + t_ndn           * Hpbe_t;


      ec = alphac*ec_lda + malphac*Hpbe;

      fncup  = ec + (alphac*ec_lda_nup + malphac*Hpbe_nup);
      fncdn  = ec + (alphac*ec_lda_ndn + malphac*Hpbe_ndn);

      xce[i]       = x_parameter*ex     + c_parameter*ec;
      fn[i]        = x_parameter*fnxup  + c_parameter*fncup;
      fn[i+n2ft3d] = x_parameter*fnxdn  + c_parameter*fncdn;

      fdn[i]          = x_parameter*fdnxup;
      fdn[i+n2ft3d]   = x_parameter*fdnxdn;
      fdn[i+2*n2ft3d] = c_parameter*t_agr*Hpbe_t*malphac;
   }
}
      


/************************************
 *                                  *
 *        gen_BEEF_BW_restricted    *
 *                                  *
 ************************************/
/*
    This routine calculates the non-Langreth terms of the BEEF-vdw exchange-correlation 
    potential(xcp) and energy density(xce).
 
 
    Entry - n2ft3d     : number of grid points
            rho_in(*) :  density (nup+ndn)
            agr_in(*): |grad rho_in|
            x_parameter: scale parameter for exchange
            c_parameter: scale parameter for correlation
            alphac: scale parameter for BEEF

    Exit  - xce(n2ft3d) : PBE96 exchange correlation energy density
             fn(n2ft3d)  : d(n*xce)/dn
             fdn(n2ft3d) : d(n*xce/d|grad n|
*/
void gen_BEEF_BW_restricted(const int n2ft3d,
                             double *rho_in, double *agr_in,
                             const double x_parameter, const double c_parameter, const double alphac,
                             double *xce, double *fn, double *fdn)
{

   /* local variables */
   double n,agr;
   double kf,ks,s,P0,n_onethird;
   double rs,rss,t,t2,t4,t6;
   double Q0,Q1,Q2,Q3,Q4,Q5,Q8,Q9,B;
   double Ht;
   double B_ec,Hrs,H_B;
   double F,Fs;

   double ex_lda,ec_lda;
   double ec_lda_rs;
   double ex,ec,H;
   double fnx,fdnx,fnc,fdnc;

   double pi         = 4.00*atan(1.00);
   double rs_scale   = pow((0.750/pi),onethird);
   double fdnx_const = -3.00/(8.00*pi);
   double  malphac=1.00-alphac;

   double p0,p1,p2,dp,s2,oneovers2,Ft,dp2,sgn;

   for (int i=0; i<n2ft3d; ++i)
   {
      n   = rho_in[i]+ETA;
      agr = agr_in[i];

      /* calculate unpolarized Exchange energies and potentials */
      n_onethird = pow((3.00*n/pi),onethird);
      ex_lda     = -0.750*n_onethird;

      kf = pow((3.00*pi*pi*n),onethird);
      s  = agr/(2.00*kf*n);

      t  = 2.00*s*s/(4.00+s*s) - 1.00;
      if (std::abs(t)>1.00) t = (t>0.0) ? 1.0 : -1.0;
      F  = 0.0;
      Ft = 0.0;
      s2   = t*t - 1.0;
      if (std::abs(s2)<1.0e-12)
      {
         if (t>0.0) 
         {
            for (int j=0; j<30; ++j)
            {
               F  = F  + am[j];
               Ft = Ft + am[j]*0.5*((double) (j*(j+1)));
            }
         }
         else
         {
            sgn = 1.00;
            for (int j=0; j<30; ++j)
            {
               F  = F  + sgn*am[j];
               Ft = Ft + sgn*am[j]*0.50*((double) (j*(j+1)));
               sgn = -sgn;
            }
         }
      }
      else
      {
         oneovers2 = 1.00/s2;
         p0 = 1.00;
         dp = 0.0;
         F  = F  + am[0]*p0;
         Ft = Ft + am[0]*dp;
         p1 = t;
         dp = 1.00;
         F  = F  + am[1]*t;
         Ft = Ft + am[1]*dp;
         for (int j=2; j<30; ++j)
         {
            p2    = (1.00/((double) (j+1))) *((2*j+1)*t*p1 - ((double) j)*p0);
            dp2   = ((double) j)*oneovers2*(t*p2-p1);
            F  = F  + am[j]*p2;
            Ft = Ft + am[j]*dp2;
            p0 = p1;
            p1 = p2;
         }
      }
      Fs  = (16.00*s/pow((4.00+s*s),2)) * Ft;

      ex   = ex_lda*F;
      fnx  = fourthird*(ex - ex_lda*Fs*s);
      fdnx = fdnx_const*Fs;


      /* calculate unpolarized correlation energies and potentials */

      /* calculate rs and t */
      rs    = rs_scale/pow(n,onethird);
      rss   = sqrt(rs);

      kf = pow((3.00*pi*pi*n),onethird);
      ks = sqrt(4.00*kf/pi);
      t  = agr/(2.00*ks*n);



      /**** unpolarized LDA correlation energy ****
       **** ec_p = correlation energy          ****
       ****   ec_p_rs = dec_p/drs              ****
       ****   uc_p    = dec_p/dn               ****/
      LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,&ec_lda,&ec_lda_rs);

      /* PBE96 correlation energy  corrections */
      t2 = t*t;
      t4 = t2*t2;
      B = -ec_lda/GAMMA;
      B = BOG/(exp(B)-1.00+ETA);
      Q4 = 1.00 + B*t2;
      Q5 = 1.00 + B*t2 + B*B*t4;
      H = GAMMA*log(1.00 + BOG*Q4*t2/Q5);

      /* PBE96 correlation fdn and fdnc derivatives */
      t6   = t4*t2;

      B_ec = (B/BETA)*(BOG+B);

      Q8  = Q5*Q5+BOG*Q4*Q5*t2;
      Q9  = 1.00+2*B*t2;
      H_B  = -BETA*B*t6*(2.00+B*t2)/Q8 ;
      Hrs  = H_B*B_ec*ec_lda_rs;

      Ht  = 2.00*BETA*Q9/Q8*t;

      ec   = alphac*ec_lda + malphac*H;
      fnc = ec  - alphac*(onethird*rs*ec_lda_rs)
                - malphac*(onethird*rs*Hrs)
                - malphac*(sevensixths*t*Ht);
      fdnc = malphac*(0.50* Ht/ks);

      xce[i] = x_parameter*ex   + c_parameter*ec;
      fn[i]  = x_parameter*fnx  + c_parameter*fnc;
      fdn[i] = x_parameter*fdnx + c_parameter*fdnc;

   }

}

}
