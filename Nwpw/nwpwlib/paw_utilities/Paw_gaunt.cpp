/* Paw_gaunt.cpp - 
   Author - Eric Bylaska
*/

#include        <iostream>
#include        <cstdio>
#include        <cstdlib>
#include        <cstring>
#include        <cmath>
#include        "blas.h"

using namespace std;


#include        "util.hpp"
#include        "Paw_gaunt.hpp"

namespace pwdft {
using namespace pwdft;

/*******************************************
 *                                         *
 *            gaunt_sub_cmplx              *
 *                                         *
 *******************************************/
/* If iscmplx=.true. this routine computes the Gaunt coefficient
*
*  Gaunt(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | dconjg(Y_lm(theta,phi)) * Y_l1m1(theta,phi) * dconjg(Y_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*      /2pi  / pi
*     |     |
*   = |     | Y_lm(theta,phi) * dconjg(Y_l1m1(theta,phi)) * Y_l2m2(theta,phi)  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*   = sqrt( (2*l+1)*(2*l2+1)/(4*pi*(2*l1+1)) ) * Clebsch(l l2 l1; m m2 m1) * Clebsch(l l2 l1; 0 0 0)
*/
static double gaunt_sub_cmplx(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3)
{
   double coeff;
   double pi = 4.0*atan(1.0);
   double twopi   = 2.00*pi;
   double fourpi  = 4.00*pi;
   double piover2 = pi/2.00;

   // **** Error Checking ****
   if ((l1<0) || (l2<0) || (l3<0))
      std::cout << "Invalid parameter in gen_gaunt_coeff, negative l" << std::endl;

   if ((l1<abs(m1)) || (l2<abs(m2)) || (l3<abs(m3)))
      std::cout  << "Invalid parameter in gen_gaunt_coeff, m > l" << std::endl;

   if (((-m1)+m2-m3) != 0) 
      return 0.0;
    else
       coeff = twopi;

   // **** Check the triangle rule ****
   if ((l3>(l1+l2)) || (l3<abs(l1-l2))) return 0.0;

   // **** Check if the integrand is odd function==>integral is zero ****
   if (((l1+l2+l3)%2) == 1)  return 0.0;

   // **** hANDLE THE EXEPTIONAL CASE ****
   if ((l1==0) && (l2==0) && (l3==0)) return (1.00/sqrt(fourpi));

   int order = l1+l2+l3;
   double x[order],w[order];
   double x1 = -1.0
   double x2 =  1.0

   // **** Generate weights and coordinates for Gauss-Legendre integration ****
   util_gauss_weights(x1,x2,x,w,order);
   double gg = 0.0;
   for (auto i=0; i<order; ++i)
      gg += w[i]*util_ytheta_lm(l1,m1,x[i])
                *util_ytheta_lm(l2,m2,x[i])
                *util_ytheta_lm(l3,m3,x[i]);
   return (coeff*gg);
}


/*******************************************
 *                                         *
 *             gaunt_sub_real              *
 *                                         *
 *******************************************/
/* or if iscmplx=.false. this routine computes the Taunt coefficient
*
*  Taunt(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | T_lm(theta,phi)) * T_l1m1(theta,phi) * T_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*/
static double gaunt_sub_real(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3)
{
   double pi = 4.0*atan(1.0);
   double twopi   = 2.00*pi;
   double fourpi  = 4.00*pi;
   double piover2 = pi/2.00;

   // **** Error Checking ****
   if ((l1<0) || (l2<0) || (l3<0))
      std::cout << "Invalid parameter in gen_gaunt_coeff, negative l" << std::endl;

   if ((l1<abs(m1)) || (l2<abs(m2)) || (l3<abs(m3)))
      std::cout  << "Invalid parameter in gen_gaunt_coeff, m > l" << std::endl;

   int mm[3] = {m1,m2,m3};
   int n=3;
   while (n>1)
   {
      int newn = 0
      for (auto i=1; i<n; ++i)
      {
         if (mm[i-1] < mm[i])
         {
            int itmp = mm[i-1];
            mm[i-1]  = mm[i];
            mm[i]    = itmp;
            newn     = i;
         }
      }
      n = newn;
   }
   int tp = (abs(mm[0])==(abs(mm[1])+abs(mm[2])) );
   int tm = (abs(mm[0])==abs((abs(mm[1])-abs(mm[2]))) );

   if ((mm[0]>0)       && (mm[1]>0)  && (mm[2]>0) && tp)
      coeff = piover2;
   else if ((mm[0]>0)  && (mm[1]>0)  && (mm[2]==0) && tp)
      coeff = pi;
   else if ((mm[0]>0)  && (mm[1]<0)  && (mm[2]<0) && tp)
      coeff = -piover2;
   else if ((mm[0]>0)  && (mm[1]<0)  && (mm[2]<0) && tm)
      coeff = piover2;
   else if ((mm[0]==0) && (mm[1]==0) && (mm[2]==0)) 
      coeff = twopi;
   else if ((mm[0]==0) && (mm[1]<0)  && (mm[2]<0) && (tm || tp))
      coeff = pi;
   else
       return 0.0;

   // **** Check the triangle rule **** 
   if ((l3>(l1+l2)) || (l3<abs(l1-l2))) return 0.0;

   // **** Check if the integrand is odd function==>integral is zero ****
   if (((l1+l2+l3)%2) == 1)  return 0.0;

   // **** hANDLE THE EXEPTIONAL CASE ****
   if ((l1==0) && (l2==0) && (l3==0)) return (1.00/sqrt(fourpi));

   int order = l1+l2+l3;
   double x[order],w[order];
   double x1 = -1.0
   double x2 =  1.0

   // **** Generate weights and coordinates for Gauss-Legendre integration ****
   util_gauss_weights(x1,x2,x,w,order);
   double gg = 0.0;   
   for (auto i=0; i<order; ++i)
      gg += w[i]*util_rtheta_lm(l1,m1,x[i])
                *util_rtheta_lm(l2,m2,x[i])
                *util_rtheta_lm(l3,m3,x[i]);
   return (coeff*gg);
}


/* Constructors */

/*******************************************
 *                                         *
 *         Paw_gaunt::Paw_gaunt            *
 *                                         *
 *******************************************/
Paw_gaunt::Paw_gaunt(const bool iscmplx, const int lmax0)
{
   gaunt_iscmplx = iscmplx
   gaunt_lmax    = lmax0
   int   sizel   = (gaunt_lmax+1);

   gaunt_sizel2 = sizel*sizel;

   gaunt_coeff = new double [gaunt_sizel2*gaunt_sizel2*gaunt_sizel2];

   for (auto l1=0; l1<=gaunt_lmax; ++l1)
   for (auto l2=0; l2<=gaunt_lmax; ++l2)
   for (auto l3=0; l3<=gaunt_lmax; ++l3)
   {
      for (auto m1=-l1; m1<=l1; ++m1)
      for (auto m2=-l2; m2<=l2; ++m2)
      for (auto m3=-l3; m3<=l3; ++m3)
      {
         int i = l1*l1 + (l1+m1);
         int j = l2*l2 + (l2+m2);
         int k = l3*l3 + (l3+m3);
         int indx = i + j*gaunt_sizel2 + k*gaunt_sizel2*gaunt_sizel2;
         
         if (gaunt_iscmplx)
            gaunt_coeff[indx] = gaunt_sub_cmplx(l1,m1,l2,m2,l3,m3);
         else
            gaunt_coeff[indx] = gaunt_sub_real(l1,m1,l2,m2,l3,m3);
      }
   }
}


/*******************************************
 *                                         *
 *            Paw_gaunt::gaunt             *
 *                                         *
 *******************************************/

double Paw_gaunt::gaunt(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3)
{

   double gg;

   if ((l1>gaunt_lmax) || (l2>gaunt_lmax) || (l3>gaunt_lmax))
   {
      if (gaunt_iscmplx)
         gg = gaunt_sub_cmplx(l1,m1,l2,m2,l3,m3);
      else
         gg = gaunt_sub_real(l1,m1,l2,m2,l3,m3);
   }
   else
   {
      int i = l1*l1 + (l1+m1);
      int j = l2*l2 + (l2+m2);
      int k = l3*l3 + (l3+m3);
      int indx = i + j*gaunt_sizel2 + k*gaunt_sizel2*nwpw_gaunt_sizel2;
      gg = dbl_mb(nwpw_gaunt_coeff(ic,1)+indx)
   }
   return gg;
}



/*******************************************
 *                                         *
 *            Paw_gaunt::gaunt2            *
 *                                         *
 *******************************************/
/*
* If iscmplx=.true. this routine computes the Gaunt coefficient
*
*  Gaunt2(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | dconjg(dY_lm(theta,phi)/dtheta) * dY_l1m1(theta,phi)/dtheta * dconjg(Y_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*      /2pi  / pi
*     |     |
*   = |     | dY_lm(theta,phi)/dtheta * dconjg(dY_l1m1(theta,phi)/dtheta) * Y_l2m2(theta,phi)  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*
* or if iscmplx=.false. this routine computes the Taunt coefficient
*
*  Taunt2(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | dT_lm(theta,phi))/dtheta * dT_l1m1(theta,phi)/dtheta * T_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*/
double Paw_gaunt::gaunt2(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3)
{

      logical iscmplx
      integer l1,l2,l3
      integer m1,m2,m3
      integer i
      integer order
      real*8  x1,x2
     
*     *** local variables ****
      double precision pi
      double precision piover2,twopi,fourpi
c      parameter(pi = 3.14159265358979323846264338327950288419d0)
c      parameter(piover2= pi*0.5d0)
c      parameter(twopi  = pi*2.0d0)
c      parameter(fourpi = pi*4.0d0)

c     work arrays for integration
      double precision x(100),w(100),coeff

      logical tp,tm
      integer mm(3),newn,n,itmp
      real*8   rtheta_lm,rtheta2_lm,ytheta_lm,ytheta2_lm
      external rtheta_lm,rtheta2_lm,ytheta_lm,ytheta2_lm

      pi = 4.0d0*datan(1.0d0)
      twopi   = 2.0d0*pi
      fourpi  = 4.0d0*pi
      piover2 = pi/2.0d0

      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) 
     > write(*,*)'Invalid parameter in nwpw_gaunt2, negative l'

      If (l1.lt.(abs(m1)+1).or.l2.lt.(abs(m2)+1).or.l3.lt.abs(m3)) then
        nwpw_gaunt2 = 0.0d0
        return
      end if

      !**** Do integration over angle phi ****
      if (iscmplx) then
         if ((-m1) + m2 - m3 .ne. 0) then
            nwpw_gaunt2 = 0.0d0
            return 
         else
           coeff = twopi
         endif
      else
         mm(1) = m1
         mm(2) = m2
         mm(3) = m3
         n = 3
         do while (n>1) 
            newn = 0
            do i=1,n-1
               if (mm(i).lt.mm(i+1)) then
                  itmp    = mm(i) 
                  mm(i)   = mm(i+1)
                  mm(i+1) = itmp
                  newn    = i
               end if
            end do
            n = newn
         end do
         tp = (abs(mm(1)).eq.(abs(mm(2))+abs(mm(3))) )
         tm = (abs(mm(1)).eq.abs((abs(mm(2))-abs(mm(3)))) )

         if ((mm(1).gt.0).and.(mm(2).gt.0).and.(mm(3).gt.0).and.tp) then
            coeff = piover2
         else if ((mm(1).gt.0).and.(mm(2).gt.0)
     >                        .and.(mm(3).eq.0).and.tp) then
            coeff = pi
         else if ((mm(1).gt.0).and.(mm(2).lt.0) 
     >                        .and.(mm(3).lt.0).and.tp) then
            coeff = -piover2
         else if ((mm(1).gt.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.tm) then
            coeff = piover2
         else if ((mm(1).eq.0).and.(mm(2).eq.0)
     >                        .and.(mm(3).eq.0)) then
            coeff = twopi
         else if ((mm(1).eq.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.(tm.or.tp)) then
            coeff = pi
         else
            nwpw_gaunt2 = 0.0d0
            return 
         end if
      endif

      !**** Check the triangle rule ****
      if (l3.gt.l1+l2 .or. l3.lt.abs(l1-l2)) then
         nwpw_gaunt2 = 0.0d0
         return 
      endif

      !**** Check if the integrand is odd function==>integral is zero ****
      if (mod(l1 + l2 + l3,2) .eq. 1) then
         nwpw_gaunt2 = 0.0d0
         return 
      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l1.eq.0 .and. l2.eq.0 .and. l3.eq.0) then
         nwpw_gaunt2 = 0.0d0
         return 
      endif
      x1 = -1.0
      x2 =  1.0
      order = l1 + l2 + l3

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL nwpw_gauss_weights(x1,x2,x,w,order)
      nwpw_gaunt2 = 0.0d0
      if (iscmplx) then
         do i = 1, order
            nwpw_gaunt2 = nwpw_gaunt2
     >                 + w(i)*ytheta2_lm(l1,m1,x(i)) 
     >                       *ytheta2_lm(l2,m2,x(i))
     >                       *ytheta_lm(l3,m3,x(i))
         end do
      else
         do i = 1, order
            nwpw_gaunt2 = nwpw_gaunt2
     >                 + w(i)*rtheta2_lm(l1,m1,x(i)) 
     >                       *rtheta2_lm(l2,m2+1,x(i))
     >                       *rtheta_lm(l3,m3,x(i))
         end do
      end if

      nwpw_gaunt2 = nwpw_gaunt2*coeff

}



/*******************************************
 *                                         *
 *            Paw_gaunt::gaunt3            *
 *                                         *
 *******************************************/
/*
* If iscmplx=.true. this routine computes the Gaunt coefficient
*
*  Gaunt3(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | dconjg(dY_lm(theta,phi)/dphi) * dY_l1m1(theta,phi)/dphi * dconjg(Y_l2m2(theta,phi)) * 1/(sin(theta))**2 *  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*      /2pi  / pi
*     |     |
*   = |     | dY_lm(theta,phi)/dphi * dconjg(dY_l1m1(theta,phi)/dphi) * Y_l2m2(theta,phi) 1/(sin(theta))**2  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*
* or if iscmplx=.false. this routine computes the Taunt coefficient
*
*  Taunt3(l,m,l1,m1,l2,m2) =
*
*      /2pi  / pi
*     |     |
*   = |     | (1/sin(theta))**2 * dT_lm(theta,phi))/dphi * dT_l1m1(theta,phi)/dphi * T_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*/

double Paw_gaunt::gaunt3(const int l1, const int m1, const int l2, const int m2, const int l3, const int m3)
{
      implicit none
      logical iscmplx
      integer l1,l2,l3
      integer m1,m2,m3
      integer i
      integer order
      real*8  x1,x2
     
*     *** local variables ****
      double precision pi
      double precision piover2,twopi,fourpi
c      parameter(pi = 3.14159265358979323846264338327950288419d0)
c      parameter(piover2= pi*0.5d0)
c      parameter(twopi  = pi*2.0d0)
c      parameter(fourpi = pi*4.0d0)

c     work arrays for integration
      double precision x(100),w(100),coeff

      logical tp,tm
      integer mm(3),newn,n,itmp
      real*8   rtheta_lm,rtheta_lm_div,ytheta_lm,theta_lm_div
      external rtheta_lm,rtheta_lm_div,ytheta_lm,theta_lm_div

      pi = 4.0d0*datan(1.0d0)
      twopi   = 2.0d0*pi
      fourpi  = 4.0d0*pi
      piover2 = pi/2.0d0

      !**** Error Checking ****
      if (l1.lt.0 .or. l2.lt.0 .or. l3.lt.0) 
     > write(*,*)'Invalid parameter in nwpw_gaunt3, negative l'
      If (l1.lt.abs(m1) .or. l3.lt.abs(m3) .or. l2.lt.abs(m2))
     1 write(*,*) 'Invalid parameter in nwpw_gaunt3, m > l'

      !**** Do integration over angle phi ****
      if (iscmplx) then
         if ((-m1) + m2 - m3 .ne. 0) then
            nwpw_gaunt3 = 0.0d0
            return 
         else
           coeff = twopi
         endif
      else
         m1 = -m1
         m2 = -m2
         mm(1) = m1
         mm(2) = m2
         mm(3) = m3
         n = 3
         do while (n>1) 
            newn = 0
            do i=1,n-1
               if (mm(i).lt.mm(i+1)) then
                  itmp    = mm(i) 
                  mm(i)   = mm(i+1)
                  mm(i+1) = itmp
                  newn    = i
               end if
            end do
            n = newn
         end do
         tp = (abs(mm(1)).eq.(abs(mm(2))+abs(mm(3))) )
         tm = (abs(mm(1)).eq.abs((abs(mm(2))-abs(mm(3)))) )

         if ((mm(1).gt.0).and.(mm(2).gt.0).and.(mm(3).gt.0).and.tp) then
            coeff = piover2
         else if ((mm(1).gt.0).and.(mm(2).gt.0)
     >                        .and.(mm(3).eq.0).and.tp) then
            coeff = pi
         else if ((mm(1).gt.0).and.(mm(2).lt.0) 
     >                        .and.(mm(3).lt.0).and.tp) then
            coeff = -piover2
         else if ((mm(1).gt.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.tm) then
            coeff = piover2
         else if ((mm(1).eq.0).and.(mm(2).eq.0)
     >                        .and.(mm(3).eq.0)) then
            coeff = twopi
         else if ((mm(1).eq.0).and.(mm(2).lt.0)
     >                        .and.(mm(3).lt.0).and.(tm.or.tp)) then
            coeff = pi
         else
            nwpw_gaunt3 = 0.0d0
            return 
         end if
      endif

c      !**** Check the triangle rule ****
c      if (l3.gt.l1+l2 .or. l3.lt.abs(l1-l2)) then
c         nwpw_gaunt3 = 0.0d0
c         return 
c      endif

c      !**** Check if the integrand is odd function==>integral is zero ****
c      if (mod(l1 + l2 + l3,2) .eq. 1) then
c         nwpw_gaunt3 = 0.0d0
c         return 
c      endif

      !**** hANDLE THE EXEPTIONAL CASE ****
      if (l1.eq.0 .and. l2.eq.0 .and. l3.eq.0) then
         nwpw_gaunt3 = 0.0d0
         return 
      endif
      x1 = -1.0
      x2 =  1.0
      order = l1 + l2 + l3

      !**** Generate weights and coordinates for Gauss-Legendre integration ****
      CALL nwpw_gauss_weights(x1,x2,x,w,order)
      nwpw_gaunt3 = 0.0d0
      if (iscmplx) then
         do i = 1, order
            nwpw_gaunt3 = nwpw_gaunt3
     >                 - w(i)*theta_lm_div(l1,m1,x(i)) * (m1)
     >                       *theta_lm_div(l2,m2,x(i)) * (m2)
     >                       *ytheta_lm(l3,m3,x(i))
         end do
      else
         do i = 1, order
            nwpw_gaunt3 = nwpw_gaunt3
     >                 + w(i)*rtheta_lm_div(l1,-m1,x(i)) * (m1)
     >                       *rtheta_lm_div(l2,-m2,x(i)) * (m2)
     >                       *rtheta_lm(l3,m3,x(i))
         end do
      end if

      nwpw_gaunt3 = nwpw_gaunt3*coeff

}



}
