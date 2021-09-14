/* hsepbe.cpp
  Author - Eric Bylaska
*/

#include	<cmath>

namespace pwdft {
using namespace pwdft;

/* Constants for HJS hole */
#define AA	7.57211e-1
#define BB	-1.06364e-1
#define CC	-1.18649e-1
#define DD	6.09650e-1
#define EE	-4.77963e-2

/*     Constants for fit of H(s) (PBE hole) */
/*     Taken from JCTC_5_754 (2009)         */
#define ha2	1.59941e-2
#define ha3	8.52995e-2
#define ha4	-1.60368e-1
#define ha5	1.52645e-1
#define ha6	-9.71263e-2
#define ha7	4.22061e-2

#define hb1	5.33319e0
#define hb2	-12.4780e0
#define hb3	11.0988e0
#define hb4	-5.11013e0
#define hb5	1.71468e0
#define	hb6	-6.10380e-1
#define hb7	3.07555e-1
#define hb8	-7.70547e-2
#define	hb9	3.34840e-2


/* Whole numbers used during evaluation */
#define zero	0.0
#define one	1.0
#define two	2.0
#define three	3.0
#define four	4.0
#define five	5.0
#define six	6.0
#define seven	7.0
#define eight	8.0
#define	nine	9.0
#define	ten	10.0
       
#define r11	11.0
#define r12	12.0
#define r14	14.0
#define	r15	15.0
#define r16	16.0
#define r18	18.0
#define r20	20.0
#define r24	24.0
#define r27	27.0
#define r30	30.0
#define r32	32.0

#define r35	35.0
#define r42	42.0
#define r48	48.0
#define r56	56.0
#define r64	64.0
#define r72	72.0

/* Fractions used during evaluation */
#define f12	0.50

#define tollz	1.0e-16

/**********************************************
 *                                            *
 *               FsrxHSE                      *
 *                                            *
 **********************************************/
/*
  HSE evaluates the Heyd et al. Screened Coulomb
  Exchange Functional
 
  Calculates the enhancement factor
*/
static void FsrxHSE(double s, double nu, double *Fxhse, double *d10Fxhse, double *d01Fxhse)
{
   double smax,strans,sconst;

   double H,hnum,hden;
   double d1H,d1hnum,d1hden;
   double s2,s3,s4,s5,s6,s7,s8,s9;
   double Fs,d1Fs;
   double zeta, lambda, eta, kf, chi, lambda2;
   double d1zeta,d1lambda,d1eta,d1nu,d1chi,d1lambda2;
   double EGs,d1EGs;
   double nu2,L2,L3,nu3,nu4,nu5,nu6;
   double Js,Ks,Ms,Ns;
   double d1Js,d1Ks,d1Ms,d1Ns;

   double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
   double tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15;
   double Fxhse1,Fxhse2,Fxhse3,Fxhse4,Fxhse5,Fxhse6;
   double d1Fxhse1,d1Fxhse2,d1Fxhse3,d1Fxhse4,d1Fxhse5;
   double d1Fxhse6,d1Fxhse7;

   double pi,pi2,srpi,s02;
   double f13,f32,f52,f72,f92;
   double faczeta;


   /* General constants */
   f13   = one/three;
   f32   = three/two;
   f52   = five/two;
   f72   = seven/two;
   f92   = nine/two;
   //pi    = acos(-one);
   pi    = 4.0*atan(one);
   pi2   = pi*pi;
   srpi = sqrt(pi);

   /* Calculate prelim variables */
   s2 = s*s;
   s02 = s2/four;
   s3 = s2*s;
   s4 = s3*s;
   s5 = s4*s;
   s6 = s5*s;
   s7 = s6*s;
   s8 = s7*s;
   s9 = s8*s;


   /* Calculate H(s) the model exhange hole */
   hnum = ha2*s2 + ha3*s3 + ha4*s4 + ha5*s5 + ha6*s6 + ha7*s7;
   hden = one + hb1*s + hb2*s2 + hb3*s3 + hb4*s4 + hb5*s5 +
          hb6*s6 + hb7*s7 + hb8*s8 + hb9*s9;
   H = hnum/hden;


   /* Calculate helper variables */
   zeta = s2*H;
   eta = AA + zeta;
   lambda = DD + zeta;

   chi = nu/sqrt(lambda+pow(nu,two));
   lambda2 = (one+chi)*(lambda+pow(nu,two));


   /* Calculate F(H(s)) for the model exhange hole */
   Fs = one-s2/(r27*CC*(one+s02))-zeta/(two*CC);

   /* Calculate EG(s) */
   EGs = -(two/five)*CC*Fs*lambda - (four/r15)*BB*pow(lambda,two) -
         (six/five)*AA*pow(lambda,three) - 
         (four/five)*srpi*pow(lambda,(seven/two)) -
         (r12/five)*(pow(lambda,(seven/two)))*(sqrt(zeta)-sqrt(eta));
 
   /* Calculate the denominators needed */
   nu2 = nu*nu;
   Js = (sqrt(zeta+nu2)+sqrt(eta+nu2))*(sqrt(zeta+nu2)+nu);
   Ks = (sqrt(zeta+nu2)+sqrt(eta+nu2))*(sqrt(eta+nu2)+nu);
   Ms = (sqrt(zeta+nu2)+sqrt(lambda+nu2))*(sqrt(lambda+nu2)+nu);
   Ns = (sqrt(eta+nu2)+sqrt(lambda+nu2))*(sqrt(lambda+nu2)+nu);


   /* The final value for the enhancement factor is */
   tmp1 = one + f12*chi;
   tmp2 = one + (nine/eight)*chi + (three/eight)*pow(chi,two);
   Fxhse1  = AA*(zeta/Js + eta/Ks);
   Fxhse2  = -(four/nine)*BB/lambda2;
   Fxhse3  = -(four/nine)*CC*Fs*tmp1/pow(lambda2,two);
   Fxhse4  = -(eight/nine)*EGs*tmp2/pow(lambda2,three);
   Fxhse5  = two*zeta*log(one -DD/Ms);
   Fxhse6  = -two*eta*log(one -(DD-AA)/Ns);

   *Fxhse = Fxhse1+Fxhse2+Fxhse3+Fxhse4+Fxhse5+Fxhse6;


   /* Calculate the first derivative of H with respect to the reduced density gradient, s.*/
   d1hnum = two*ha2*s + three*ha3*s2 + four*ha4*s3 +
             five*ha5*s4 + six*ha6*s5 + seven*ha7*s6;

   d1hden  = hb1 + two*hb2*s +three*hb3*s2 + four*hb4*s3 +
             five*hb5*s4 + six*hb6*s5 + seven*hb7*s6 +
             eight*hb8*s7 + nine*hb9*s8;

   d1H =   ((d1hnum -hnum*(d1hden/hden))/hden);

   /* calculate first derivative of variables needed with respect to s */
   d1zeta = two*s*H + s2*d1H;
   d1eta  = d1zeta;
   d1lambda = d1zeta;
   d1chi = -f12*nu*d1zeta/pow((lambda + nu2),f32);
   d1lambda2 = d1chi*(lambda + pow(nu,two)) + (one+chi)*d1lambda;

   /* calculate the first derivative of Fs with respect to s */
   d1Fs = -two*s/(r27*CC*pow((one+s02),two)) - d1zeta/(two*CC);

   /* Calculate the first derivate of EGs with respect to s */
   faczeta = 0.0;
   if(std::abs(zeta)>tollz) faczeta = faczeta + d1zeta/sqrt(zeta);
   faczeta = faczeta - d1eta/sqrt(eta);
   d1EGs = -(two/five)*CC*(d1Fs*lambda + Fs*d1lambda) -
           (eight/r15)*BB*lambda*d1lambda -
           (r18/five)*AA*lambda*lambda*d1lambda -
           (r14/five)*srpi*d1lambda*pow(lambda,f52) -
           (r42/five)*(pow(lambda,f52))*
           d1lambda*(sqrt(zeta)-sqrt(eta))-
           (six/five)*(pow(lambda,(seven/two)))*
           faczeta;


   /* Calculate the first derivate of denominators needed with respect to s */
   tmp1 = (sqrt(zeta+nu2)+nu)/(sqrt(eta+nu2));
   tmp2 = (sqrt(eta+nu2)+nu)/(sqrt(zeta+nu2));

   d1Js = f12*d1zeta*(two+tmp1+tmp2);
   d1Ks = d1Js;

   tmp3 = (sqrt(zeta+nu2)+nu)/(sqrt(lambda+nu2));
   tmp4 = (sqrt(lambda+nu2)+nu)/(sqrt(zeta+nu2));
   d1Ms = f12*d1zeta*(two +tmp3+tmp4);

   tmp5 = (sqrt(lambda+nu2)+nu)/(sqrt(eta+nu2));
   tmp6 = (sqrt(eta+nu2)+nu)/(sqrt(lambda+nu2));
   d1Ns = f12*d1zeta*(two + tmp5+tmp6);


   /* Calculate the derivative of the 08-Fxhse with respect to s */
   L2 = lambda2*lambda2;
   L3 = lambda2*lambda2*lambda2;
   d1Fxhse1  = AA*( (Js*d1zeta - zeta*d1Js)/(Js*Js) +
                    (Ks*d1zeta - eta*d1Ks)/(Ks*Ks) );

   d1Fxhse2  = (four/nine)*BB*d1lambda2/L2;

   tmp9 = d1lambda2/lambda2;
   tmp7 = d1Fs - two*Fs*tmp9;
   tmp8 = one + f12*chi;
   tmp10 =  f12*Fs*d1chi;

   d1Fxhse3 = -(four*CC/(nine*L2))*(tmp7*tmp8+tmp10);

   tmp7 = one + (nine/eight)*chi+(three/eight)*chi*chi;
   tmp8 = (nine/eight)*d1chi + (six/eight)*chi*d1chi;

   d1Fxhse4 = -(eight/(nine*L3))*((d1EGs-three*EGs*tmp9)*tmp7 + EGs*tmp8);

   d1Fxhse5  = two*d1zeta*log(one-DD/Ms) + two*zeta*DD*d1Ms/(Ms*Ms*(one-DD/Ms));

   d1Fxhse6  = -two*d1eta*log(one- (DD-AA)/Ns) - two*eta*(DD-AA)*d1Ns/(Ns*Ns*(one-(DD-AA)/Ns));
 
   *d10Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6;


   /* Calculate the derivative of 08-Fxhse with respect to nu */
   nu3 = nu2*nu;

   d1Fxhse1 = -((AA*(nu + sqrt(eta + nu2))*zeta)/
               (sqrt(eta + nu2)*sqrt(nu2 + zeta)*
               (nu + sqrt(nu2 + zeta))*
               (sqrt(eta + nu2) + sqrt(nu2 + zeta))));

   d1Fxhse2 = -((AA*eta*(nu/sqrt(eta + nu2) + nu/
               sqrt(nu2 + zeta)))/
               ((nu + sqrt(eta + nu2))*
               pow((sqrt(eta + nu2) + sqrt(nu2 + zeta)),two))) -
               (AA*eta*(one + nu/sqrt(eta + nu2)))/
               (pow((nu + sqrt(eta + nu2)),two)*
               (sqrt(eta + nu2) + sqrt(nu2 + zeta)));

   d1Fxhse3 = (four*BB)/(nine*pow((lambda + nu2),(f32)));

   d1Fxhse4 = (two*CC*Fs)/(three*pow((lambda + nu2),(f52)));

   d1Fxhse5 = (five*EGs*(pow(lambda,two) + four*nu3*
               (nu + sqrt(lambda + nu2)) +
               lambda*nu*(five*nu + three*sqrt(lambda + nu2))))/
      (three*pow((lambda + nu2),four)* pow((nu + sqrt(lambda + nu2)),three));

   d1Fxhse6 = (two*DD*zeta*(nu + sqrt(nu2 + zeta)))/
                 (sqrt(lambda + nu2)*sqrt(nu2 + zeta)*
                 (-DD + lambda + (nu + sqrt(lambda + nu2))*
                 (nu + sqrt(nu2 + zeta))));

   d1Fxhse7 = (two*(AA - DD)*eta*(nu + sqrt(eta + nu2)))/
              (sqrt(eta + nu2)*sqrt(lambda + nu2)*
              (AA - DD + lambda + nu2 + nu*sqrt(eta + nu2) +
              nu*sqrt(lambda + nu2) +
              sqrt(eta + nu2)*sqrt(lambda + nu2)));

   *d01Fxhse = d1Fxhse1+d1Fxhse2+d1Fxhse3+d1Fxhse4+d1Fxhse5+d1Fxhse6+d1Fxhse7;

}


/* Density cutoff parameters */
#define DNS_CUT	1.0e-20
#define ETA	1.0e-20
#define ETA2	1.0e-14
#define alpha_zeta	(1.0e0-ETA2)
#define alpha_zeta2	(1.0e0-ETA2)

/* HSEPBE GGA exchange constants */
#define MU	0.2195149727645171e0
#define KAPPA	0.8040000000000000e0
#define OMEGA	0.2070000000000000e0

/* PBE96 GGA correlation constants */
#define GAMMA	0.031090690869655e0
#define BETA	0.066724550603149e0
#define BOG	(BETA/GAMMA)



/* Perdew-Wang92 LDA correlation coefficients */
#define GAM	0.519842099789746329e0
#define iGAM	(1.0/GAM)
#define FZZ	(8.0/(9.0*GAM)) 
#define iFZZ	(0.125*9.0*GAM)

#define A_1	0.0310907e0
#define A1_1	0.2137000e0
#define B1_1	7.5957000e0
#define B2_1	3.5876000e0
#define B3_1	1.6382000e0
#define B4_1	0.4929400e0

#define A_2 	0.01554535e0
#define A1_2	0.20548000e0
#define B1_2	14.11890000e0
#define B2_2	6.19770000e0
#define B3_2	3.36620000e0
#define B4_2	0.62517000e0

#define A_3 	0.0168869e0
#define A1_3	0.1112500e0
#define B1_3	10.3570000e0
#define B2_3	3.6231000e0
#define B3_3	0.8802600e0
#define B4_3	0.4967100e0

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
#define onethird	(1.00/3.00)
#define onethirdm	(-1.00/3.00)
#define twothird	(2.00/3.00)
#define fourthird	(4.00/3.00)
#define fivethird	(5.00/3.00)
#define onesixthm	(-1.00/6.00)
#define sevensixthm	(-7.00/6.00)
#define sevensixths	(7.00/6.00)

/****************************************
 *					*
 *	 gen_HSE_BW_unrestricted	*
 *					*
 ****************************************/
/*
     This function returns the PBE96 exchange-correlation
   energy density, xce, and its derivatives with respect
   to nup, ndn, |grad nup|, |grad ndn|, and |grad n|.
 
    Entry - n2ft3d     : number of grid points
            dn_in(*,2) : spin densites nup and ndn
            agr_in(*,3): |grad nup|, |grad ndn|, and |grad n|
            x_parameter: scale parameter for exchange
            c_parameter: scale parameter for correlation
 
    Exit - xce(*)  : PBE96 energy density
         - fn(*,2) : d(n*xce)/dnup, d(n*xce)/dndn
         - fdn(*,3): d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|
                     d(n*xce)/d|grad n|
*/
void gen_HSE_BW_unrestricted(const int n2ft3d,
                             double *dn_in, double *agr_in,
                             const double x_parameter, const double c_parameter,
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

   double pi = 4.00*atan(1.00);
   double rs_scale = pow((0.750/pi),onethird);
   double fdnx_const = -3.00/(8.00*pi);

   double nu,Fxhse,Fxhse_s,Fxhse_nu;

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
      P0 = 1.00 + (MU/KAPPA)*s*s;

      F   = (1.00 + KAPPA - KAPPA/P0);
      Fs  = 2.00*MU/(P0*P0)*s;

      /* shortrange-HSEPBE96 */
      nu = OMEGA/kf;
      FsrxHSE(s,nu,&Fxhse,&Fxhse_s,&Fxhse_nu);

      exup = ex_lda*(F - 0.25*Fxhse);
      fnxup = fourthird*(exup - ex_lda*(Fs - 0.20*Fxhse_s)*s);
      fnxup = fnxup + 0.25*onethird*ex_lda*Fxhse_nu*nu;
      fdnxup = fdnx_const*(Fs - 0.25*Fxhse_s);

     /**** down ****/
      n     = 2.00*ndn;
      agr   = 2.00*agrdn;

      n_onethird = pow((3.00*n/pi),onethird);
      ex_lda     = -0.750*n_onethird;

      kf = pow((3.00*pi*pi*n),onethird);
      s  = agr/(2.00*kf*n);
      P0 = 1.00 + (MU/KAPPA)*s*s;

      F   = (1.00 + KAPPA - KAPPA/P0);
      Fs  = 2.00*MU/(P0*P0)*s;

      /* shortrange-HSEPBE96 */
      nu = OMEGA/kf;
      FsrxHSE(s,nu,&Fxhse,&Fxhse_s,&Fxhse_nu);

      exdn   = ex_lda*(F - 0.25*Fxhse);
      fnxdn  = fourthird*(exdn - ex_lda*(Fs-0.25*Fxhse_s)*s);
      fnxdn  = fnxdn + 0.25*onethird*ex_lda*Fxhse_nu*nu;
      fdnxdn = fdnx_const*(Fs - 0.25*Fxhse_s);


      n = nup+ndn;

      ex = (exup*nup + exdn*ndn)/ n;
                  
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



      ec = ec_lda + Hpbe;

      fncup  = ec + (ec_lda_nup + Hpbe_nup);
      fncdn  = ec + (ec_lda_ndn + Hpbe_ndn);

      xce[i]       = x_parameter*ex     + c_parameter*ec;
      fn[i]        = x_parameter*fnxup  + c_parameter*fncup;
      fn[i+n2ft3d] = x_parameter*fnxdn  + c_parameter*fncdn;

      fdn[i]          = x_parameter*fdnxup;
      fdn[i+n2ft3d]   = x_parameter*fdnxdn;
      fdn[i+2*n2ft3d] = c_parameter*t_agr*Hpbe_t;

   }
}


/****************************************
 *					*
 *	    gen_HSE_BW_restricted	*
 *					*
 ****************************************/
/*
   This routine calculates the PBE96 exchange-correlation 
   potential(xcp) and energy density(xce).

   Entry - n2ft3d     : number of grid points
           rho_in(*) :  density (nup+ndn)
           agr_in(*): |grad rho_in|
           x_parameter: scale parameter for exchange
           c_parameter: scale parameter for correlation

     Exit  - xce(n2ft3d) : PBE96 exchange correlation energy density
             fn(n2ft3d)  : d(n*xce)/dn
             fdn(n2ft3d) : d(n*xce/d|grad n| 
*/
void gen_HSE_BW_restricted(const int n2ft3d,
                           double *rho_in, double *agr_in,
                           const double x_parameter, const double c_parameter,
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

   double nu,Fxhse,Fxhse_s,Fxhse_nu;

   double pi         = 4.00*atan(1.00);
   double rs_scale   = pow((0.750/pi),onethird);
   double fdnx_const = -3.00/(8.00*pi);

      
   for (int i=0; i<n2ft3d; ++i)
   {
      n   = rho_in[i]+ETA;
      agr = agr_in[i];
        
      /* calculate unpolarized Exchange energies and potentials */
      n_onethird = pow((3.00*n/pi),onethird);
      ex_lda     = -0.750*n_onethird;

      kf = pow((3.00*pi*pi*n),onethird);
      s  = agr/(2.00*kf*n);
      P0 = 1.00 + (MU/KAPPA)*s*s;

      F   = (1.00 + KAPPA - KAPPA/P0);
      Fs  = 2.00*MU/(P0*P0)*s;

      /* shortrange-HSEPBE96 */
      nu = OMEGA/kf;
      FsrxHSE(s,nu,&Fxhse,&Fxhse_s,&Fxhse_nu);

      ex   = ex_lda*(F - 0.20*Fxhse);
      fnx  = fourthird*(ex - ex_lda*(Fs - 0.20*Fxhse_s)*s);
      fnx  = fnx + 0.250*onethird*ex_lda*Fxhse_nu*nu;
      fdnx = fdnx_const*(Fs - 0.250*Fxhse_s);


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

      ec   = ec_lda + H;
      fnc = ec  - (onethird*rs*ec_lda_rs)
                - (onethird*rs*Hrs)
                - (sevensixths*t*Ht);
      fdnc = 0.50* Ht/ks;

      xce[i] = x_parameter*ex   + c_parameter*ec;
      fn[i]  = x_parameter*fnx  + c_parameter*fnc;
      fdn[i] = x_parameter*fdnx + c_parameter*fdnc;
   }
}

      
}



