
#include        <cmath>
#include        <cstdlib>
#include	"nwpw_timing.hpp"
#include	"v_exc.hpp"

namespace pwdft {
using namespace pwdft;

/**********************************************
 *                                            *
 *                v_exc                       *              *
 *                                            *
 **********************************************/

/* computes the vosko lda xc energy density and potential
*/

/*---- parameters given by vosko et al -----------------*/
#define ap   3.109070e-02
#define af   1.554530e-02
#define x0p -1.049800e-01
#define x0f -3.250000e-01
#define bp   3.727440e+00
#define bf   7.060420e+00
#define cp   1.293520e+01
#define cf   1.805780e+01
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
#define fc    1.923661e+00
#define fd    2.564881e+00
#define crs   7.876233e-01

/* other constants */
#define one3rd (1.00/3.00)
#define for3rd (4.00/3.00)
#define one6th (1.00/6.00)
#define dncut  1.0e-20

void v_exc(const int ispin, const int n2ft3d, double *dn, 
           double *xcp, double *xce, double *x)
{
   nwpw_timing_start(4);

   /* local variables */
   int    k;
   double pi,xx,xx1,rho,zup,zdw,f,df;
   double *rhoup,*xcpup,*xceup;
   double *rhodn,*xcpdn,*xcedn;

   pi      = 4.00*atan(1.00);

   /* define arrays and such */
   rhoup = dn;   rhodn = &dn[ (ispin-1)*n2ft3d];
   xcpup = xcp;  xcpdn = &xcp[(ispin-1)*n2ft3d];
   xceup = xce;  xcedn = &xce[(ispin-1)*n2ft3d];
   

   /* square root of wigner radius */
   for (k=0; k<n2ft3d; ++k)
   {
        rho=rhoup[k]+rhodn[k]+dncut;
        x[k]=crs/pow(rho,one6th);
   }
  

   /* paramagnetic correlation energy & potential */
   for (k=0; k<n2ft3d; ++k)
   {
      xx=1.00/(x[k]*(x[k]+bp)+cp);
      xx1 = x[k] + cp3;
      xx1 *= xx1;
      xceup[k]= cp1*log(xx*x[k]*x[k])
            + cp2*log(xx*xx1)
            + cp4*atan(cp5/(x[k]+cp6));

      xx1 = x[k]+dp6;
      xx1 *= xx1;
      xcpup[k]=xceup[k]
            -one6th*x[k]*(
                 dp1/x[k]+dp2/(x[k]+dp3)
                +dp4*xx*(2.00*x[k]+bp)
                +dp5/(xx1+dp7) );
   }

   /* paramagnetic exchange energy & potential */
   for (k=0; k<n2ft3d; ++k)
   {
      xceup[k]=xceup[k]+(xp/(x[k]*x[k]));
      xcpup[k]=xcpup[k]+for3rd*(xp/(x[k]*x[k]));
   }

   /* do unrestricted part */
   if (ispin==2)
   {
      /* ferromagnetic correlation energy & potential */
      for (k=0; k<n2ft3d; ++k)
      {
        xx=1.00/(x[k]*(x[k]+bf)+cf);
        xx1 = x[k]+cf6;
        xx1 *= xx1;
        xcedn[k]=cf1*log(xx*x[k]*x[k])
                     +cf2*log(xx*xx1)
                     +cf4*atan(cf5/(x[k]+cf6));

        xx1 = x[k]+df6;
        xx1 *= xx1;
        xcpdn[k]=xcedn[k]
                     -one6th*x[k]*(
                        df1/x[k]+df2/(x[k]+df3)
                       +df4*xx*(2.00*x[k]+bf)
                       +df5/(xx1+df7) );
      }

      /* ferromagnetic exchange-energy & potential */
      for (k=0; k<n2ft3d; ++k)
      {
         xcedn[k]=xcedn[k]+(xf/(x[k]*x[k]));
         xcpdn[k]=xcpdn[k]+for3rd*(xf/(x[k]*x[k]));
      }

      /* spin polarized exchange-correlation potential */
      for (k=0; k<n2ft3d; ++k)
      {
         rho=rhoup[k]+rhodn[k]+dncut;
         zup=2.00*rhoup[k]/rho;
         zdw=2.00*rhodn[k]/rho;
         f=(zup*(pow(zup,one3rd))+zdw*(pow(zdw,one3rd))-2.00)*fc;
         xcpup[k]=(1.00-f)*xcpup[k]+f*xcpdn[k];
         df=(pow(zup,one3rd)-pow(zdw,one3rd))*(xceup[k]-xcedn[k])*fd;
         xcpdn[k]=xcpup[k]-zup*df;
         xcpup[k]=xcpup[k]+zdw*df;
         xceup[k]=xceup[k]+f*(xcedn[k]-xceup[k]);
      }
   }

   nwpw_timing_end(4);
}

}

