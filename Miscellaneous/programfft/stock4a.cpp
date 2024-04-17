
#include 	<cmath>
#include 	<iostream>
#include 	<complex>
#include 	<stdexcept>  // Required for std::runtime_error

typedef std::complex<double> complex_t;

void fft0_fac2(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    //std::cout << "n=" << n << " s=" << s << " eo=" << eo  << std::endl;
    const int m = n/2;
    const double theta0 = 2*M_PI/n;
    //const complex_t u12 = complex_t(cos(2*M_PI/2.0), isgn*sin(2*M_PI/2.0));

    if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
    else {
        for (int p = 0; p < m; p++) {
            const complex_t wp = complex_t(cos(p*theta0), isgn*sin(p*theta0));
            for (int q = 0; q < s; q++) {
                const complex_t a = x[q + s*(p + 0)];
                const complex_t b = x[q + s*(p + m)];
                y[q + s*(2*p + 0)] =  a + b;
                y[q + s*(2*p + 1)] = (a - b) * wp;
                //y[q + s*(2*p + 1)] = (a + b*u12) * wp;
            }
        }
        //fft0_fac2(n/2, 2*s, !eo, y, x);
    }
}

void fft0_fac3(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
   const int m = n/3;
   const double theta0 = 2*M_PI/n;
   const complex_t u13 = complex_t(cos(2*M_PI/3.0), isgn*sin(2*M_PI/3.0));
   const complex_t u23 = complex_t(cos(4*M_PI/3.0), isgn*sin(4*M_PI/3.0));
   std::cout << "FAC3 U=" << u13 << " " << u23 << std::endl;
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              y[q + s*(3*p + 0)] =  a + b + c;
              y[q + s*(3*p + 1)] = (a + b*u13 + c*u23) * wp;
              y[q + s*(3*p + 2)] = (a + b*u23 + c*u13) * wp2;
          }
      }
      //fft0_fac3(n/2, 2*s, !eo, y, x);
   }
}

         
void fft0_fac3_twiddle(int isgn, int n, int s, bool eo, complex_t* twiddle, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{             
   const int m = n/3;
   const double theta0 = 2*M_PI/n; 
   //const complex_t u13 = complex_t(cos(2*M_PI/3.0), isgn*sin(2*M_PI/3.0));
   //const complex_t u23 = complex_t(cos(4*M_PI/3.0), isgn*sin(4*M_PI/3.0));

   complex_t u13 = twiddle[n+1];
   complex_t u23 = twiddle[n+2];
   if (isgn<0)
   {
      u13 = conj(u13);
      u23 = conj(u23);
   }
    

   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          //const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          //const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));

          complex_t wp  = twiddle[3*p+1];
          complex_t wp2 = twiddle[3*p+2];
          if (isgn<0)
          {
             wp  = conj(wp);
             wp2 = conj(wp2);
          }

          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              y[q + s*(3*p + 0)] =  a + b + c;
              y[q + s*(3*p + 1)] = (a + b*u13 + c*u23) * wp;
              y[q + s*(3*p + 2)] = (a + b*u23 + c*u13) * wp2;
          }
      }
      //fft0_fac3(n/2, 2*s, !eo, y, x);
   }
}

void fft0_fac4(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{  
   const int m = n/4;
   const double theta0 = 2*M_PI/n;
   const complex_t u14 = complex_t(cos(2*M_PI/4.0), isgn*sin(2*M_PI/4.0)); //i
   const complex_t u24 = complex_t(cos(4*M_PI/4.0), isgn*sin(4*M_PI/4.0)); //-1
   const complex_t u34 = complex_t(cos(6*M_PI/4.0), isgn*sin(6*M_PI/4.0)); //-i
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else { 
      for (int p = 0; p < m; p++) {
          const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          for (int q = 0; q < s; q++) { 
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              y[q + s*(4*p + 0)] =  a + b + c + d;
              y[q + s*(4*p + 1)] = (a + b*u14 + c*u24 + d*u34) * wp;
              y[q + s*(4*p + 2)] = (a + b*u24 + c     + d*u24) * wp2;
              y[q + s*(4*p + 3)] = (a + b*u34 + c*u24 + d*u14) * wp3;
          }
      }
      //fft0_fac4(n/4, 4*s, !eo, y, x);
   }
}


void fft0_fac4_twiddle(int isgn, int n, int s, bool eo, complex_t* twiddle, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{ 
   const int m = n/4;
   const double theta0 = 2*M_PI/n;
   //const complex_t u14 = complex_t(cos(2*M_PI/4.0), isgn*sin(2*M_PI/4.0)); //i
   //const complex_t u24 = complex_t(cos(4*M_PI/4.0), isgn*sin(4*M_PI/4.0)); //-1
   //const complex_t u34 = complex_t(cos(6*M_PI/4.0), isgn*sin(6*M_PI/4.0)); //-i
   complex_t u14 = twiddle[n+1];
   complex_t u24 = twiddle[n+2];
   complex_t u34 = twiddle[n+3];
   if (isgn<0)
   {
      u14 = conj(u14);
      u24 = conj(u24);
      u34 = conj(u34);
   }
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          //const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          //const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          //const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          complex_t wp  = twiddle[4*p+1];
          complex_t wp2 = twiddle[4*p+2];
          complex_t wp3 = twiddle[4*p+3];
          if (isgn<0)
          {
             wp  = conj(wp);
             wp2 = conj(wp2);
             wp3 = conj(wp3);
          }
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              y[q + s*(4*p + 0)] =  a + b + c + d;
              y[q + s*(4*p + 1)] = (a + b*u14 + c*u24 + d*u34) * wp;
              y[q + s*(4*p + 2)] = (a + b*u24 + c     + d*u24) * wp2;
              y[q + s*(4*p + 3)] = (a + b*u34 + c*u24 + d*u14) * wp3;
          }
      }
      //fft0_fac4(n/4, 4*s, !eo, y, x);
   }
}




void fft0_fac5(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
   const int m = n/5;
   const double theta0 = 2*M_PI/n;
   const complex_t u15 = complex_t(cos(2*M_PI/5.0), isgn*sin(2*M_PI/5.0));
   const complex_t u25 = complex_t(cos(4*M_PI/5.0), isgn*sin(4*M_PI/5.0));
   const complex_t u35 = complex_t(cos(6*M_PI/5.0), isgn*sin(6*M_PI/5.0));
   const complex_t u45 = complex_t(cos(8*M_PI/5.0), isgn*sin(8*M_PI/5.0));

   std::cout << "FAC5 U=" << u15 << " " << u25 << " " << u35 << " " << u45 <<  std::endl;
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          const complex_t wp4 = complex_t(cos(4*p*theta0), isgn*sin(4*p*theta0));
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              const complex_t e = x[q + s*(p + 4*m)];
              y[q + s*(5*p + 0)] =  a + b + c + d + e;
              y[q + s*(5*p + 1)] = (a + b*u15 + c*u25 + d*u35 + e*u45) * wp;
              y[q + s*(5*p + 2)] = (a + b*u25 + c*u45 + d*u15 + e*u35) * wp2;
              y[q + s*(5*p + 3)] = (a + b*u35 + c*u15 + d*u45 + e*u25) * wp3;
              y[q + s*(5*p + 4)] = (a + b*u45 + c*u35 + d*u25 + e*u15) * wp4;
          }
      }
      //fft0_fac5(n/5, 5*s, !eo, y, x);
   }
}


void fft0_fac5_twiddle(int isgn, int n, int s, bool eo, complex_t* twiddle, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
   const int m = n/5;
   const double theta0 = 2*M_PI/n;
   //const complex_t u15 = complex_t(cos(2*M_PI/5.0), isgn*sin(2*M_PI/5.0));
   //const complex_t u25 = complex_t(cos(4*M_PI/5.0), isgn*sin(4*M_PI/5.0));
   //const complex_t u35 = complex_t(cos(6*M_PI/5.0), isgn*sin(6*M_PI/5.0));
   //const complex_t u45 = complex_t(cos(8*M_PI/5.0), isgn*sin(8*M_PI/5.0));
   complex_t u15 = twiddle[n+1];
   complex_t u25 = twiddle[n+2];
   complex_t u35 = twiddle[n+3];
   complex_t u45 = twiddle[n+4];
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          //const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          //const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          //const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          //const complex_t wp4 = complex_t(cos(4*p*theta0), isgn*sin(4*p*theta0));
          complex_t wp  = twiddle[5*p+1];
          complex_t wp2 = twiddle[5*p+2];
          complex_t wp3 = twiddle[5*p+3];
          complex_t wp4 = twiddle[5*p+4];
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              const complex_t e = x[q + s*(p + 4*m)];
              y[q + s*(5*p + 0)] =  a + b + c + d + e;
              y[q + s*(5*p + 1)] = (a + b*u15 + c*u25 + d*u35 + e*u45) * wp;
              y[q + s*(5*p + 2)] = (a + b*u25 + c*u45 + d*u15 + e*u35) * wp2;
              y[q + s*(5*p + 3)] = (a + b*u35 + c*u15 + d*u45 + e*u25) * wp3;
              y[q + s*(5*p + 4)] = (a + b*u45 + c*u35 + d*u25 + e*u15) * wp4;
          }
      }
      //fft0_fac5(n/5, 5*s, !eo, y, x);
   }
}



void fft0_fac6(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
   const int m = n/6;
   const double theta0 = 2*M_PI/n;
   const complex_t u16 = complex_t(cos( 2*M_PI/6.0), isgn*sin( 2*M_PI/6.0)); //  0.5000 + 0.8660i
   const complex_t u26 = complex_t(cos( 4*M_PI/6.0), isgn*sin( 4*M_PI/6.0)); // -0.5000 + 0.8660i
   const complex_t u36 = complex_t(cos( 6*M_PI/6.0), isgn*sin( 6*M_PI/6.0)); // -1.0
   const complex_t u46 = complex_t(cos( 8*M_PI/6.0), isgn*sin( 8*M_PI/6.0)); // -0.5000 - 0.8660i
   const complex_t u56 = complex_t(cos(10*M_PI/6.0), isgn*sin(10*M_PI/6.0)); //  0.5000 - 0.8660i
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          const complex_t wp4 = complex_t(cos(4*p*theta0), isgn*sin(4*p*theta0));
          const complex_t wp5 = complex_t(cos(5*p*theta0), isgn*sin(5*p*theta0));
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              const complex_t e = x[q + s*(p + 4*m)];
              const complex_t f = x[q + s*(p + 5*m)];
              y[q + s*(6*p + 0)] =  a + b + c + d + e + f;;
              y[q + s*(6*p + 1)] = (a + b*u16 + c*u26 + d*u36 + e*u46 + f*u56) * wp;
              y[q + s*(6*p + 2)] = (a + b*u26 + c*u46 + d     + e*u26 + f*u46) * wp2;
              y[q + s*(6*p + 3)] = (a + b*u36 + c     + d*u36 + e     + f*u36) * wp3;
              y[q + s*(6*p + 4)] = (a + b*u46 + c*u26 + d     + e*u46 + f*u26) * wp4;
              y[q + s*(6*p + 5)] = (a + b*u56 + c*u46 + d*u36 + e*u26 + f*u16) * wp5;
          }
      }
      //fft0_fac6(n/6, 6*s, !eo, y, x);
   }
}

void fft0_fac6_twiddle(int isgn, int n, int s, bool eo, complex_t* twiddle,  complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
   const int m = n/6;
   const double theta0 = 2*M_PI/n;
   //const complex_t u16 = complex_t(cos( 2*M_PI/6.0), isgn*sin( 2*M_PI/6.0)); //  0.5000 + 0.8660i
   //const complex_t u26 = complex_t(cos( 4*M_PI/6.0), isgn*sin( 4*M_PI/6.0)); // -0.5000 + 0.8660i
   //const complex_t u36 = complex_t(cos( 6*M_PI/6.0), isgn*sin( 6*M_PI/6.0)); // -1.0
   //const complex_t u46 = complex_t(cos( 8*M_PI/6.0), isgn*sin( 8*M_PI/6.0)); // -0.5000 - 0.8660i
   //const complex_t u56 = complex_t(cos(10*M_PI/6.0), isgn*sin(10*M_PI/6.0)); //  0.5000 - 0.8660i
   complex_t u16 = twiddle[n+1];
   complex_t u26 = twiddle[n+2];
   complex_t u36 = twiddle[n+3];
   complex_t u46 = twiddle[n+4];
   complex_t u56 = twiddle[n+5];
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          //const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          //const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          //const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          //const complex_t wp4 = complex_t(cos(4*p*theta0), isgn*sin(4*p*theta0));
          //const complex_t wp5 = complex_t(cos(5*p*theta0), isgn*sin(5*p*theta0));
          complex_t wp  = twiddle[6*p+1];
          complex_t wp2 = twiddle[6*p+2];
          complex_t wp3 = twiddle[6*p+3];
          complex_t wp4 = twiddle[6*p+4];
          complex_t wp5 = twiddle[6*p+5];
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              const complex_t e = x[q + s*(p + 4*m)];
              const complex_t f = x[q + s*(p + 5*m)];
              y[q + s*(6*p + 0)] =  a + b + c + d + e + f;;
              y[q + s*(6*p + 1)] = (a + b*u16 + c*u26 + d*u36 + e*u46 + f*u56) * wp;
              y[q + s*(6*p + 2)] = (a + b*u26 + c*u46 + d     + e*u26 + f*u46) * wp2;
              y[q + s*(6*p + 3)] = (a + b*u36 + c     + d*u36 + e     + f*u36) * wp3;
              y[q + s*(6*p + 4)] = (a + b*u46 + c*u26 + d     + e*u46 + f*u26) * wp4;
              y[q + s*(6*p + 5)] = (a + b*u56 + c*u46 + d*u36 + e*u26 + f*u16) * wp5;
          }
      }
      //fft0_fac6(n/6, 6*s, !eo, y, x);
   }
}



void fft0_fac7(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
   const int m = n/7;
   const double theta0 = 2*M_PI/n;
   const complex_t u17 = complex_t(cos( 2*M_PI/7.0), isgn*sin( 2*M_PI/7.0));
   const complex_t u27 = complex_t(cos( 4*M_PI/7.0), isgn*sin( 4*M_PI/7.0));
   const complex_t u37 = complex_t(cos( 6*M_PI/7.0), isgn*sin( 6*M_PI/7.0));
   const complex_t u47 = complex_t(cos( 8*M_PI/7.0), isgn*sin( 8*M_PI/7.0));
   const complex_t u57 = complex_t(cos(10*M_PI/7.0), isgn*sin(10*M_PI/7.0));
   const complex_t u67 = complex_t(cos(12*M_PI/7.0), isgn*sin(12*M_PI/7.0));
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          const complex_t wp4 = complex_t(cos(4*p*theta0), isgn*sin(4*p*theta0));
          const complex_t wp5 = complex_t(cos(5*p*theta0), isgn*sin(5*p*theta0));
          const complex_t wp6 = complex_t(cos(6*p*theta0), isgn*sin(6*p*theta0));
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              const complex_t e = x[q + s*(p + 4*m)];
              const complex_t f = x[q + s*(p + 5*m)];
              const complex_t g = x[q + s*(p + 6*m)];
              y[q + s*(7*p + 0)] =  a + b + c + d + e + f + g;
              y[q + s*(7*p + 1)] = (a + b*u17 + c*u27 + d*u37 + e*u47 + f*u57 + g*u67) * wp;
              y[q + s*(7*p + 2)] = (a + b*u27 + c*u47 + d*u67 + e*u17 + f*u37 + g*u57) * wp2;
              y[q + s*(7*p + 3)] = (a + b*u37 + c*u67 + d*u27 + e*u57 + f*u17 + g*u47) * wp3;
              y[q + s*(7*p + 4)] = (a + b*u47 + c*u17 + d*u57 + e*u27 + f*u67 + g*u37) * wp4;
              y[q + s*(7*p + 5)] = (a + b*u57 + c*u37 + d*u17 + e*u67 + f*u47 + g*u27) * wp5;
              y[q + s*(7*p + 6)] = (a + b*u67 + c*u57 + d*u47 + e*u37 + f*u27 + g*u17) * wp6;
          }
      }
      //fft0_fac7(n/7, 7*s, !eo, y, x);
   }
}

void fft0_fac7_twiddle(int isgn, int n, int s, bool eo, complex_t* twiddle, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride  
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{           
   const int m = n/7;
   const double theta0 = 2*M_PI/n;
   //const complex_t u17 = complex_t(cos( 2*M_PI/7.0), isgn*sin( 2*M_PI/7.0));
   //const complex_t u27 = complex_t(cos( 4*M_PI/7.0), isgn*sin( 4*M_PI/7.0));
   //const complex_t u37 = complex_t(cos( 6*M_PI/7.0), isgn*sin( 6*M_PI/7.0));
   //const complex_t u47 = complex_t(cos( 8*M_PI/7.0), isgn*sin( 8*M_PI/7.0));
   //const complex_t u57 = complex_t(cos(10*M_PI/7.0), isgn*sin(10*M_PI/7.0));
   //const complex_t u67 = complex_t(cos(12*M_PI/7.0), isgn*sin(12*M_PI/7.0));
   complex_t u17 = twiddle[n+1];
   complex_t u27 = twiddle[n+2];
   complex_t u37 = twiddle[n+3];
   complex_t u47 = twiddle[n+4];
   complex_t u57 = twiddle[n+5];
   complex_t u67 = twiddle[n+6];
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          //const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          //const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          //const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          //const complex_t wp4 = complex_t(cos(4*p*theta0), isgn*sin(4*p*theta0));
          //const complex_t wp5 = complex_t(cos(5*p*theta0), isgn*sin(5*p*theta0));
          //const complex_t wp6 = complex_t(cos(6*p*theta0), isgn*sin(6*p*theta0));
          complex_t wp  = twiddle[7*p+1];
          complex_t wp2 = twiddle[7*p+2];
          complex_t wp3 = twiddle[7*p+3];
          complex_t wp4 = twiddle[7*p+4];
          complex_t wp5 = twiddle[7*p+5];
          complex_t wp6 = twiddle[7*p+6];
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              const complex_t e = x[q + s*(p + 4*m)];
              const complex_t f = x[q + s*(p + 5*m)];
              const complex_t g = x[q + s*(p + 6*m)];
              y[q + s*(7*p + 0)] =  a + b + c + d + e + f + g;
              y[q + s*(7*p + 1)] = (a + b*u17 + c*u27 + d*u37 + e*u47 + f*u57 + g*u67) * wp;
              y[q + s*(7*p + 2)] = (a + b*u27 + c*u47 + d*u67 + e*u17 + f*u37 + g*u57) * wp2;
              y[q + s*(7*p + 3)] = (a + b*u37 + c*u67 + d*u27 + e*u57 + f*u17 + g*u47) * wp3;
              y[q + s*(7*p + 4)] = (a + b*u47 + c*u17 + d*u57 + e*u27 + f*u67 + g*u37) * wp4;
              y[q + s*(7*p + 5)] = (a + b*u57 + c*u37 + d*u17 + e*u67 + f*u47 + g*u27) * wp5;
              y[q + s*(7*p + 6)] = (a + b*u67 + c*u57 + d*u47 + e*u37 + f*u27 + g*u17) * wp6;
          }
      }
      //fft0_fac7(n/7, 7*s, !eo, y, x);
   }
}



void fft0_fac8(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{  
   const int m = n/8;
   const double theta0 = 2*M_PI/n;
   const complex_t u18 = complex_t(cos( 2*M_PI/8.0), isgn*sin( 2*M_PI/8.0)); //  sqrt2/2 + i*sqrt2/2
   const complex_t u28 = complex_t(cos( 4*M_PI/8.0), isgn*sin( 4*M_PI/8.0)); //  i
   const complex_t u38 = complex_t(cos( 6*M_PI/8.0), isgn*sin( 6*M_PI/8.0)); // -sqrt2/2 + i*sqrt2/2
   const complex_t u48 = complex_t(cos( 8*M_PI/8.0), isgn*sin( 8*M_PI/8.0)); // -1
   const complex_t u58 = complex_t(cos(10*M_PI/8.0), isgn*sin(10*M_PI/8.0)); // -sqrt2/2 - i*sqrt2/2
   const complex_t u68 = complex_t(cos(12*M_PI/8.0), isgn*sin(12*M_PI/8.0)); // -i
   const complex_t u78 = complex_t(cos(14*M_PI/8.0), isgn*sin(14*M_PI/8.0)); // sqrt2/2 - i*sqrt2/2
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else { 
      for (int p = 0; p < m; p++) {
          const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          const complex_t wp4 = complex_t(cos(4*p*theta0), isgn*sin(4*p*theta0));
          const complex_t wp5 = complex_t(cos(5*p*theta0), isgn*sin(5*p*theta0));
          const complex_t wp6 = complex_t(cos(6*p*theta0), isgn*sin(6*p*theta0));
          const complex_t wp7 = complex_t(cos(7*p*theta0), isgn*sin(7*p*theta0));
          for (int q = 0; q < s; q++) { 
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              const complex_t e = x[q + s*(p + 4*m)];
              const complex_t f = x[q + s*(p + 5*m)];
              const complex_t g = x[q + s*(p + 6*m)]; 
              const complex_t h = x[q + s*(p + 7*m)]; 
              y[q + s*(8*p + 0)] =  a + b + c + d + e + f + g + h;
              y[q + s*(8*p + 1)] = (a + b*u18 + c*u28 + d*u38 + e*u48 + f*u58 + g*u68 + h*u78) * wp;
              y[q + s*(8*p + 2)] = (a + b*u28 + c*u48 + d*u68 + e     + f*u28 + g*u48 + h*u68) * wp2;
              y[q + s*(8*p + 3)] = (a + b*u38 + c*u68 + d*u18 + e*u48 + f*u78 + g*u28 + h*u58) * wp3;
              y[q + s*(8*p + 4)] = (a + b*u48 + c     + d*u48 + e     + f*u48 + g     + h*u48) * wp4;
              y[q + s*(8*p + 5)] = (a + b*u58 + c*u28 + d*u78 + e*u48 + f*u18 + g*u68 + h*u38) * wp5;
              y[q + s*(8*p + 6)] = (a + b*u68 + c*u48 + d*u28 + e     + f*u68 + g*u48 + h*u28) * wp6;
              y[q + s*(8*p + 7)] = (a + b*u78 + c*u68 + d*u58 + e*u48 + f*u38 + g*u28 + h*u18) * wp7;
          }
      }
      //fft0_fac8(n/8, 8*s, !eo, y, x);
   }
}

void fft0_fac8_twiddle(int isgn, int n, int s, bool eo, complex_t* twiddle, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
   const int m = n/8;
   const double theta0 = 2*M_PI/n;
   //const complex_t u18 = complex_t(cos( 2*M_PI/8.0), isgn*sin( 2*M_PI/8.0)); //  sqrt2/2 + i*sqrt2/2
   //const complex_t u28 = complex_t(cos( 4*M_PI/8.0), isgn*sin( 4*M_PI/8.0)); //  i
   //const complex_t u38 = complex_t(cos( 6*M_PI/8.0), isgn*sin( 6*M_PI/8.0)); // -sqrt2/2 + i*sqrt2/2
   //const complex_t u48 = complex_t(cos( 8*M_PI/8.0), isgn*sin( 8*M_PI/8.0)); // -1
   //const complex_t u58 = complex_t(cos(10*M_PI/8.0), isgn*sin(10*M_PI/8.0)); // -sqrt2/2 - i*sqrt2/2
   //const complex_t u68 = complex_t(cos(12*M_PI/8.0), isgn*sin(12*M_PI/8.0)); // -i
   //const complex_t u78 = complex_t(cos(14*M_PI/8.0), isgn*sin(14*M_PI/8.0)); // sqrt2/2 - i*sqrt2/2

   complex_t u18 = twiddle[n+1];
   complex_t u28 = twiddle[n+2];
   complex_t u38 = twiddle[n+3];
   complex_t u48 = twiddle[n+4];
   complex_t u58 = twiddle[n+5];
   complex_t u68 = twiddle[n+6];
   complex_t u78 = twiddle[n+7];
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          //const complex_t wp  = complex_t(cos(p*theta0),   isgn*sin(p*theta0));
          //const complex_t wp2 = complex_t(cos(2*p*theta0), isgn*sin(2*p*theta0));
          //const complex_t wp3 = complex_t(cos(3*p*theta0), isgn*sin(3*p*theta0));
          //const complex_t wp4 = complex_t(cos(4*p*theta0), isgn*sin(4*p*theta0));
          //const complex_t wp5 = complex_t(cos(5*p*theta0), isgn*sin(5*p*theta0));
          //const complex_t wp6 = complex_t(cos(6*p*theta0), isgn*sin(6*p*theta0));
          //const complex_t wp7 = complex_t(cos(7*p*theta0), isgn*sin(7*p*theta0));
          complex_t wp  = twiddle[8*p+1];
          complex_t wp2 = twiddle[8*p+2];
          complex_t wp3 = twiddle[8*p+3];
          complex_t wp4 = twiddle[8*p+4];
          complex_t wp5 = twiddle[8*p+5];
          complex_t wp6 = twiddle[8*p+6];
          complex_t wp7 = twiddle[8*p+7];
          for (int q = 0; q < s; q++) {
              const complex_t a = x[q + s*(p + 0)];
              const complex_t b = x[q + s*(p + 1*m)];
              const complex_t c = x[q + s*(p + 2*m)];
              const complex_t d = x[q + s*(p + 3*m)];
              const complex_t e = x[q + s*(p + 4*m)];
              const complex_t f = x[q + s*(p + 5*m)];
              const complex_t g = x[q + s*(p + 6*m)];
              const complex_t h = x[q + s*(p + 7*m)];
              y[q + s*(8*p + 0)] =  a + b + c + d + e + f + g + h;
              y[q + s*(8*p + 1)] = (a + b*u18 + c*u28 + d*u38 + e*u48 + f*u58 + g*u68 + h*u78) * wp;
              y[q + s*(8*p + 2)] = (a + b*u28 + c*u48 + d*u68 + e     + f*u28 + g*u48 + h*u68) * wp2;
              y[q + s*(8*p + 3)] = (a + b*u38 + c*u68 + d*u18 + e*u48 + f*u78 + g*u28 + h*u58) * wp3;
              y[q + s*(8*p + 4)] = (a + b*u48 + c     + d*u48 + e     + f*u48 + g     + h*u48) * wp4;
              y[q + s*(8*p + 5)] = (a + b*u58 + c*u28 + d*u78 + e*u48 + f*u18 + g*u68 + h*u38) * wp5;
              y[q + s*(8*p + 6)] = (a + b*u68 + c*u48 + d*u28 + e     + f*u68 + g*u48 + h*u28) * wp6;
              y[q + s*(8*p + 7)] = (a + b*u78 + c*u68 + d*u58 + e*u48 + f*u38 + g*u28 + h*u18) * wp7;
          }
      }
      //fft0_fac8(n/8, 8*s, !eo, y, x);
   }
}


void fft0_fac9(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
/* matlab code:
printf("a0 "); uu = 0; s = 8;
for i=1:8
uu = mod(uu+s,9);
printf(" + a%d*u%d" , i,uu); printf("_9")
endfor
printf("\n")
*/
{  
   const int m = n/9;
   const double theta0 = 2*M_PI/n;
   const complex_t u1_9 = complex_t(cos( 2*M_PI/9.0), isgn*sin( 2*M_PI/9.0));
   const complex_t u2_9 = complex_t(cos( 4*M_PI/9.0), isgn*sin( 4*M_PI/9.0));
   const complex_t u3_9 = complex_t(cos( 6*M_PI/9.0), isgn*sin( 6*M_PI/9.0));
   const complex_t u4_9 = complex_t(cos( 8*M_PI/9.0), isgn*sin( 8*M_PI/9.0));
   const complex_t u5_9 = complex_t(cos(10*M_PI/9.0), isgn*sin(10*M_PI/9.0));
   const complex_t u6_9 = complex_t(cos(12*M_PI/9.0), isgn*sin(12*M_PI/9.0));
   const complex_t u7_9 = complex_t(cos(14*M_PI/9.0), isgn*sin(14*M_PI/9.0));
   const complex_t u8_9 = complex_t(cos(16*M_PI/9.0), isgn*sin(16*M_PI/9.0));

   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          const complex_t wp1 = complex_t(cos(   p*theta0), isgn*sin(   p*theta0));
          const complex_t wp2 = complex_t(cos( 2*p*theta0), isgn*sin( 2*p*theta0));
          const complex_t wp3 = complex_t(cos( 3*p*theta0), isgn*sin( 3*p*theta0));
          const complex_t wp4 = complex_t(cos( 4*p*theta0), isgn*sin( 4*p*theta0));
          const complex_t wp5 = complex_t(cos( 5*p*theta0), isgn*sin( 5*p*theta0));
          const complex_t wp6 = complex_t(cos( 6*p*theta0), isgn*sin( 6*p*theta0));
          const complex_t wp7 = complex_t(cos( 7*p*theta0), isgn*sin( 7*p*theta0));
          const complex_t wp8 = complex_t(cos( 8*p*theta0), isgn*sin( 8*p*theta0));
       for (int q = 0; q < s; q++) {
              const complex_t a0 = x[q + s*(p +  0)];
              const complex_t a1 = x[q + s*(p +  1*m)];
              const complex_t a2 = x[q + s*(p +  2*m)];
              const complex_t a3 = x[q + s*(p +  3*m)];
              const complex_t a4 = x[q + s*(p +  4*m)];
              const complex_t a5 = x[q + s*(p +  5*m)];
              const complex_t a6 = x[q + s*(p +  6*m)];
              const complex_t a7 = x[q + s*(p +  7*m)];
              const complex_t a8 = x[q + s*(p +  8*m)];

              y[q + s*(9*p + 0)] = a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8;
              y[q + s*(9*p + 1)]  = (a0 + a1*u1_9 + a2*u2_9 + a3*u3_9 + a4*u4_9 + a5*u5_9 + a6*u6_9 + a7*u7_9 + a8*u8_9) * wp1;
              y[q + s*(9*p + 2)]  = (a0 + a1*u2_9 + a2*u4_9 + a3*u6_9 + a4*u8_9 + a5*u1_9 + a6*u3_9 + a7*u5_9 + a8*u7_9) * wp2;
              y[q + s*(9*p + 3)]  = (a0 + a1*u3_9 + a2*u6_9 + a3      + a4*u3_9 + a5*u6_9 + a6      + a7*u3_9 + a8*u6_9) * wp3;
              y[q + s*(9*p + 4)]  = (a0 + a1*u4_9 + a2*u8_9 + a3*u3_9 + a4*u7_9 + a5*u2_9 + a6*u6_9 + a7*u1_9 + a8*u5_9) * wp4;
              y[q + s*(9*p + 5)]  = (a0 + a1*u5_9 + a2*u1_9 + a3*u6_9 + a4*u2_9 + a5*u7_9 + a6*u3_9 + a7*u8_9 + a8*u4_9) * wp5;
              y[q + s*(9*p + 6)]  = (a0 + a1*u6_9 + a2*u3_9 + a3      + a4*u6_9 + a5*u3_9 + a6      + a7*u6_9 + a8*u3_9) * wp6;
              y[q + s*(9*p + 7)]  = (a0 + a1*u7_9 + a2*u5_9 + a3*u3_9 + a4*u1_9 + a5*u8_9 + a6*u6_9 + a7*u4_9 + a8*u2_9) * wp7;
              y[q + s*(9*p + 8)]  = (a0 + a1*u8_9 + a2*u7_9 + a3*u6_9 + a4*u5_9 + a5*u4_9 + a6*u3_9 + a7*u2_9 + a8*u1_9) * wp8;
         }
      }
      //fft0_fac9(n/9, 9*s, !eo, y, x);
   }
}

          
void fft0_fac11(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft 
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
/* matlab code:
printf("a0 ");
 uu = 0; s = 10;
for i=1:10
uu = mod(uu+s,11);
printf(" + a%d*u%d" , i,uu);
printf("_11")
endfor
printf("\n")
*/
{
   const int m = n/11;
   const double theta0 = 2*M_PI/n;
   const complex_t u1_11 = complex_t(cos( 2*M_PI/11.0), isgn*sin( 2*M_PI/11.0));
   const complex_t u2_11 = complex_t(cos( 4*M_PI/11.0), isgn*sin( 4*M_PI/11.0));
   const complex_t u3_11 = complex_t(cos( 6*M_PI/11.0), isgn*sin( 6*M_PI/11.0));
   const complex_t u4_11 = complex_t(cos( 8*M_PI/11.0), isgn*sin( 8*M_PI/11.0));
   const complex_t u5_11 = complex_t(cos(10*M_PI/11.0), isgn*sin(10*M_PI/11.0));
   const complex_t u6_11 = complex_t(cos(12*M_PI/11.0), isgn*sin(12*M_PI/11.0));
   const complex_t u7_11 = complex_t(cos(14*M_PI/11.0), isgn*sin(14*M_PI/11.0));
   const complex_t u8_11 = complex_t(cos(16*M_PI/11.0), isgn*sin(16*M_PI/11.0));
   const complex_t u9_11 = complex_t(cos(18*M_PI/11.0), isgn*sin(18*M_PI/11.0));
   const complex_t u10_11= complex_t(cos(20*M_PI/11.0), isgn*sin(20*M_PI/11.0));

   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          const complex_t wp1 = complex_t(cos(   p*theta0), isgn*sin(   p*theta0));
          const complex_t wp2 = complex_t(cos( 2*p*theta0), isgn*sin( 2*p*theta0));
          const complex_t wp3 = complex_t(cos( 3*p*theta0), isgn*sin( 3*p*theta0));
          const complex_t wp4 = complex_t(cos( 4*p*theta0), isgn*sin( 4*p*theta0));
          const complex_t wp5 = complex_t(cos( 5*p*theta0), isgn*sin( 5*p*theta0));
          const complex_t wp6 = complex_t(cos( 6*p*theta0), isgn*sin( 6*p*theta0));
          const complex_t wp7 = complex_t(cos( 7*p*theta0), isgn*sin( 7*p*theta0));
          const complex_t wp8 = complex_t(cos( 8*p*theta0), isgn*sin( 8*p*theta0));
          const complex_t wp9 = complex_t(cos( 9*p*theta0), isgn*sin( 9*p*theta0));
          const complex_t wp10= complex_t(cos(10*p*theta0), isgn*sin(10*p*theta0));
       for (int q = 0; q < s; q++) {
              const complex_t a0 = x[q + s*(p +  0)];
              const complex_t a1 = x[q + s*(p +  1*m)];
              const complex_t a2 = x[q + s*(p +  2*m)];
              const complex_t a3 = x[q + s*(p +  3*m)];
              const complex_t a4 = x[q + s*(p +  4*m)];
              const complex_t a5 = x[q + s*(p +  5*m)];
              const complex_t a6 = x[q + s*(p +  6*m)];
              const complex_t a7 = x[q + s*(p +  7*m)];
              const complex_t a8 = x[q + s*(p +  8*m)];
              const complex_t a9 = x[q + s*(p +  9*m)];
              const complex_t a10 = x[q +s*(p + 10*m)];

              y[q + s*(11*p + 0)] = a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10;
              y[q + s*(11*p + 1)]  = (a0  + a1*u1_11 + a2*u2_11 + a3*u3_11 + a4*u4_11 + a5*u5_11 + a6*u6_11 + a7*u7_11 + a8*u8_11 + a9*u9_11 + a10*u10_11) * wp1;
              y[q + s*(11*p + 2)]  = (a0  + a1*u2_11 + a2*u4_11 + a3*u6_11 + a4*u8_11 + a5*u10_11 + a6*u1_11 + a7*u3_11 + a8*u5_11 + a9*u7_11 + a10*u9_11) * wp2;
              y[q + s*(11*p + 3)]  = (a0  + a1*u3_11 + a2*u6_11 + a3*u9_11 + a4*u1_11 + a5*u4_11 + a6*u7_11 + a7*u10_11 + a8*u2_11 + a9*u5_11 + a10*u8_11) * wp3;
              y[q + s*(11*p + 4)]  = (a0  + a1*u4_11 + a2*u8_11 + a3*u1_11 + a4*u5_11 + a5*u9_11 + a6*u2_11 + a7*u6_11 + a8*u10_11 + a9*u3_11 + a10*u7_11) * wp4;
              y[q + s*(11*p + 5)]  = (a0  + a1*u5_11 + a2*u10_11 + a3*u4_11 + a4*u9_11 + a5*u3_11 + a6*u8_11 + a7*u2_11 + a8*u7_11 + a9*u1_11 + a10*u6_11) * wp5;
              y[q + s*(11*p + 6)]  = (a0  + a1*u6_11 + a2*u1_11 + a3*u7_11 + a4*u2_11 + a5*u8_11 + a6*u3_11 + a7*u9_11 + a8*u4_11 + a9*u10_11 + a10*u5_11) * wp6;
              y[q + s*(11*p + 7)]  = (a0  + a1*u7_11 + a2*u3_11 + a3*u10_11 + a4*u6_11 + a5*u2_11 + a6*u9_11 + a7*u5_11 + a8*u1_11 + a9*u8_11 + a10*u4_11) * wp7;
              y[q + s*(11*p + 8)]  = (a0  + a1*u8_11 + a2*u5_11 + a3*u2_11 + a4*u10_11 + a5*u7_11 + a6*u4_11 + a7*u1_11 + a8*u9_11 + a9*u6_11 + a10*u3_11) * wp8;
              y[q + s*(11*p + 9)]  = (a0  + a1*u9_11 + a2*u7_11 + a3*u5_11 + a4*u3_11 + a5*u1_11 + a6*u10_11 + a7*u8_11 + a8*u6_11 + a9*u4_11 + a10*u2_11) * wp9;
              y[q + s*(11*p + 10)] = (a0  + a1*u10_11 + a2*u9_11 + a3*u8_11 + a4*u7_11 + a5*u6_11 + a6*u5_11 + a7*u4_11 + a8*u3_11 + a9*u2_11 + a10*u1_11) * wp10;
         }
      }
      //fft0_fac11(n/11, 11*s, !eo, y, x);
   }
}



void fft0_fac17(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
/*/ 
printf("a00 ");
 uu = 0; s = 16;
for i=1:16
uu = mod(uu+s,17);
printf(" + a%02d*u%02d" , i,uu);
printf("_17")
endfor
printf("\n")
*/
{  
   const int m = n/17;
   const double theta0 = 2*M_PI/n;
   const complex_t u01_17 = complex_t(cos( 2*M_PI/17.0), isgn*sin( 2*M_PI/17.0));
   const complex_t u02_17 = complex_t(cos( 4*M_PI/17.0), isgn*sin( 4*M_PI/17.0));
   const complex_t u03_17 = complex_t(cos( 6*M_PI/17.0), isgn*sin( 6*M_PI/17.0));
   const complex_t u04_17 = complex_t(cos( 8*M_PI/17.0), isgn*sin( 8*M_PI/17.0));
   const complex_t u05_17 = complex_t(cos(10*M_PI/17.0), isgn*sin(10*M_PI/17.0));
   const complex_t u06_17 = complex_t(cos(12*M_PI/17.0), isgn*sin(12*M_PI/17.0));
   const complex_t u07_17 = complex_t(cos(14*M_PI/17.0), isgn*sin(14*M_PI/17.0));
   const complex_t u08_17 = complex_t(cos(16*M_PI/17.0), isgn*sin(16*M_PI/17.0));
   const complex_t u09_17 = complex_t(cos(18*M_PI/17.0), isgn*sin(18*M_PI/17.0));
   const complex_t u10_17 = complex_t(cos(20*M_PI/17.0), isgn*sin(20*M_PI/17.0));
   const complex_t u11_17 = complex_t(cos(22*M_PI/17.0), isgn*sin(22*M_PI/17.0));
   const complex_t u12_17 = complex_t(cos(24*M_PI/17.0), isgn*sin(24*M_PI/17.0));
   const complex_t u13_17 = complex_t(cos(26*M_PI/17.0), isgn*sin(26*M_PI/17.0));
   const complex_t u14_17 = complex_t(cos(28*M_PI/17.0), isgn*sin(28*M_PI/17.0));
   const complex_t u15_17 = complex_t(cos(30*M_PI/17.0), isgn*sin(30*M_PI/17.0));
   const complex_t u16_17 = complex_t(cos(32*M_PI/17.0), isgn*sin(32*M_PI/17.0));

   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else {
      for (int p = 0; p < m; p++) {
          const complex_t wp01 = complex_t(cos(   p*theta0), isgn*sin(   p*theta0));
          const complex_t wp02 = complex_t(cos( 2*p*theta0), isgn*sin( 2*p*theta0));
          const complex_t wp03 = complex_t(cos( 3*p*theta0), isgn*sin( 3*p*theta0));
          const complex_t wp04 = complex_t(cos( 4*p*theta0), isgn*sin( 4*p*theta0));
          const complex_t wp05 = complex_t(cos( 5*p*theta0), isgn*sin( 5*p*theta0));
          const complex_t wp06 = complex_t(cos( 6*p*theta0), isgn*sin( 6*p*theta0));
          const complex_t wp07 = complex_t(cos( 7*p*theta0), isgn*sin( 7*p*theta0));
          const complex_t wp08 = complex_t(cos( 8*p*theta0), isgn*sin( 8*p*theta0));
          const complex_t wp09 = complex_t(cos( 9*p*theta0), isgn*sin( 9*p*theta0));
          const complex_t wp10 = complex_t(cos(10*p*theta0), isgn*sin(10*p*theta0));
          const complex_t wp11 = complex_t(cos(11*p*theta0), isgn*sin(11*p*theta0));
          const complex_t wp12 = complex_t(cos(12*p*theta0), isgn*sin(12*p*theta0));
          const complex_t wp13 = complex_t(cos(13*p*theta0), isgn*sin(13*p*theta0));
          const complex_t wp14 = complex_t(cos(14*p*theta0), isgn*sin(14*p*theta0));
          const complex_t wp15 = complex_t(cos(15*p*theta0), isgn*sin(15*p*theta0));
          const complex_t wp16 = complex_t(cos(16*p*theta0), isgn*sin(16*p*theta0));

          for (int q = 0; q < s; q++) {
              const complex_t a00= x[q + s*(p +  0)];
              const complex_t a01= x[q + s*(p +  1*m)];
              const complex_t a02= x[q + s*(p +  2*m)];
              const complex_t a03= x[q + s*(p +  3*m)];
              const complex_t a04= x[q + s*(p +  4*m)];
              const complex_t a05= x[q + s*(p +  5*m)];
              const complex_t a06= x[q + s*(p +  6*m)];
              const complex_t a07= x[q + s*(p +  7*m)];
              const complex_t a08= x[q + s*(p +  8*m)];
              const complex_t a09= x[q + s*(p +  9*m)];
              const complex_t a10= x[q + s*(p + 10*m)];
              const complex_t a11= x[q + s*(p + 11*m)];
              const complex_t a12= x[q + s*(p + 12*m)];
              const complex_t a13= x[q + s*(p + 13*m)];
              const complex_t a14= x[q + s*(p + 14*m)];
              const complex_t a15= x[q + s*(p + 15*m)];
              const complex_t a16= x[q + s*(p + 16*m)];

              y[q + s*(17*p + 0)] = a00 + a01 + a02 + a03 + a04 + a05 + a06 + a07 + a08 + a09 + a10 + a11 + a12 + a13 + a14 + a15 + a16;

              y[q + s*(17*p + 1)]  = (a00  + a01*u01_17 + a02*u02_17 + a03*u03_17 + a04*u04_17 + a05*u05_17 + a06*u06_17 
                                           + a07*u07_17 + a08*u08_17 + a09*u09_17 + a10*u10_17 + a11*u11_17 + a12*u12_17 
                                           + a13*u13_17 + a14*u14_17 + a15*u15_17 + a16*u16_17)*wp01;

              y[q + s*(17*p + 2)]  = (a00  + a01*u02_17 + a02*u04_17 + a03*u06_17 + a04*u08_17 + a05*u10_17 + a06*u12_17
                                           + a07*u14_17 + a08*u16_17 + a09*u01_17 + a10*u03_17 + a11*u05_17 + a12*u07_17
                                           + a13*u09_17 + a14*u11_17 + a15*u13_17 + a16*u15_17)*wp02;

              y[q + s*(17*p + 3)]  = (a00  + a01*u03_17 + a02*u06_17 + a03*u09_17 + a04*u12_17 + a05*u15_17 + a06*u01_17
                                           + a07*u04_17 + a08*u07_17 + a09*u10_17 + a10*u13_17 + a11*u16_17 + a12*u02_17 
                                           + a13*u05_17 + a14*u08_17 + a15*u11_17 + a16*u14_17)*wp03;

              y[q + s*(17*p + 4)]  = (a00  + a01*u04_17 + a02*u08_17 + a03*u12_17 + a04*u16_17 + a05*u03_17 + a06*u07_17 
                                           + a07*u11_17 + a08*u15_17 + a09*u02_17 + a10*u06_17 + a11*u10_17 + a12*u14_17 
                                           + a13*u01_17 + a14*u05_17 + a15*u09_17 + a16*u13_17)*wp04;

              y[q + s*(17*p + 5)]  = (a00  + a01*u05_17 + a02*u10_17 + a03*u15_17 + a04*u03_17 + a05*u08_17 + a06*u13_17 
                                           + a07*u01_17 + a08*u06_17 + a09*u11_17 + a10*u16_17 + a11*u04_17 + a12*u09_17 
                                           + a13*u14_17 + a14*u02_17 + a15*u07_17 + a16*u12_17)*wp05;

              y[q + s*(17*p + 6)]  = (a00  + a01*u06_17 + a02*u12_17 + a03*u01_17 + a04*u07_17 + a05*u13_17 + a06*u02_17 
                                           + a07*u08_17 + a08*u14_17 + a09*u03_17 + a10*u09_17 + a11*u15_17 + a12*u04_17 
                                           + a13*u10_17 + a14*u16_17 + a15*u05_17 + a16*u11_17)*wp06;

              y[q + s*(17*p + 7)]  = (a00  + a01*u07_17 + a02*u14_17 + a03*u04_17 + a04*u11_17 + a05*u01_17 + a06*u08_17 
                                           + a07*u15_17 + a08*u05_17 + a09*u12_17 + a10*u02_17 + a11*u09_17 + a12*u16_17 
                                           + a13*u06_17 + a14*u13_17 + a15*u03_17 + a16*u10_17)*wp07;

              y[q + s*(17*p + 8)]  = (a00  + a01*u08_17 + a02*u16_17 + a03*u07_17 + a04*u15_17 + a05*u06_17 + a06*u14_17 
                                           + a07*u05_17 + a08*u13_17 + a09*u04_17 + a10*u12_17 + a11*u03_17 + a12*u11_17 
                                           + a13*u02_17 + a14*u10_17 + a15*u01_17 + a16*u09_17)*wp08;

              y[q + s*(17*p + 9)]  = (a00  + a01*u09_17 + a02*u01_17 + a03*u10_17 + a04*u02_17 + a05*u11_17 + a06*u03_17 
                                           + a07*u12_17 + a08*u04_17 + a09*u13_17 + a10*u05_17 + a11*u14_17 + a12*u06_17 
                                           + a13*u15_17 + a14*u07_17 + a15*u16_17 + a16*u08_17)*wp09;

              y[q + s*(17*p + 10)] = (a00  + a01*u10_17 + a02*u03_17 + a03*u13_17 + a04*u06_17 + a05*u16_17 + a06*u09_17 
                                           + a07*u02_17 + a08*u12_17 + a09*u05_17 + a10*u15_17 + a11*u08_17 + a12*u01_17 
                                           + a13*u11_17 + a14*u04_17 + a15*u14_17 + a16*u07_17)*wp10;

              y[q + s*(17*p + 11)] = (a00  + a01*u11_17 + a02*u05_17 + a03*u16_17 + a04*u10_17 + a05*u04_17 + a06*u15_17 
                                           + a07*u09_17 + a08*u03_17 + a09*u14_17 + a10*u08_17 + a11*u02_17 + a12*u13_17 
                                           + a13*u07_17 + a14*u01_17 + a15*u12_17 + a16*u06_17)*wp11;

              y[q + s*(17*p + 12)] = (a00  + a01*u12_17 + a02*u07_17 + a03*u02_17 + a04*u14_17 + a05*u09_17 + a06*u04_17 
                                           + a07*u16_17 + a08*u11_17 + a09*u06_17 + a10*u01_17 + a11*u13_17 + a12*u08_17 
                                           + a13*u03_17 + a14*u15_17 + a15*u10_17 + a16*u05_17)*wp12;

              y[q + s*(17*p + 13)] = (a00  + a01*u13_17 + a02*u09_17 + a03*u05_17 + a04*u01_17 + a05*u14_17 + a06*u10_17 
                                           + a07*u06_17 + a08*u02_17 + a09*u15_17 + a10*u11_17 + a11*u07_17 + a12*u03_17 
                                           + a13*u16_17 + a14*u12_17 + a15*u08_17 + a16*u04_17)*wp13;

              y[q + s*(17*p + 14)] = (a00  + a01*u14_17 + a02*u11_17 + a03*u08_17 + a04*u05_17 + a05*u02_17 + a06*u16_17 
                                           + a07*u13_17 + a08*u10_17 + a09*u07_17 + a10*u04_17 + a11*u01_17 + a12*u15_17  
                                           + a13*u12_17 + a14*u09_17 + a15*u06_17 + a16*u03_17)*wp14;

              y[q + s*(17*p + 15)] = (a00  + a01*u15_17 + a02*u13_17 + a03*u11_17 + a04*u09_17 + a05*u07_17 + a06*u05_17 
                                           + a07*u03_17 + a08*u01_17 + a09*u16_17 + a10*u14_17 + a11*u12_17 + a12*u10_17 
                                           + a13*u08_17 + a14*u06_17 + a15*u04_17 + a16*u02_17)*wp15;

              y[q + s*(17*p + 16)] = (a00  + a01*u16_17 + a02*u15_17 + a03*u14_17 + a04*u13_17 + a05*u12_17 + a06*u11_17 
                                           + a07*u10_17 + a08*u09_17 + a09*u08_17 + a10*u07_17 + a11*u06_17 + a12*u05_17 
                                           + a13*u04_17 + a14*u03_17 + a15*u02_17 + a16*u01_17)*wp16;
          }
      }
      //fft0_fac17(n/17, 17*s, !eo, y, x);
   }
}

/*****************************************
 *                                       *
 *            fft_radix                  *
 *                                       *
 *****************************************/
/**
 * @brief Performs a single radix-based FFT transformation step.
 *
 * This function computes one stage of a Fast Fourier Transform using a specified radix. It supports
 * flexible radix sizes, making it suitable for various FFT schemes (e.g., mixed radix FFT).
 * The function efficiently handles complex input and output data arrays and uses precomputed
 * twiddle factors to perform the FFT calculations.
 *
 * The transformation can be applied forward (FFT) or inverse (IFFT), depending on the sign
 * of the `isgn` parameter. This function is designed to be called multiple times to complete
 * an entire FFT computation, once for each stage of the transform.
 *
 * @param isgn The sign of the exponent in the FFT formula. Use -1 for forward FFT and +1 for inverse FFT.
 * @param n The total number of complex samples in the input array `x`. Should be divisible by `radix`.
 * @param s The stride between successive samples in the input data, allowing processing of subarrays.
 * @param eo A boolean flag indicating the output buffer:
 *           true (`eo = 1`) means `y` is the output buffer,
 *           false (`eo = 0`) means `x` is the output buffer and `y` is used as a temporary buffer.
 * @param radix The radix to be used for this stage of the FFT, determining the size of the FFT "butterflies."
 * @param twiddle A pointer to an array containing the precomputed twiddle factors necessary for
 *                the FFT calculations. The array should contain all necessary factors for the computations
 *                and post-multiplications specific to the given `radix`.
 * @param x A pointer to the input array of complex samples. This array contains `n` complex numbers.
 * @param y A pointer to the output or workspace array of complex samples. Should be at least as large as `x`.
 *
 * @note The function assumes that the size of `x` and `y` is sufficient to hold all data required for the
 *       computations, including any necessary stride or alignment considerations.
 */

void fft_radix(const int n, const int s, bool eo, const int radix, complex_t* twiddle, complex_t* x, complex_t *y)
{
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else 
   {
      std::cout << "FFT U= " ;
      for (int r1=0; r1<radix; ++ r1)
         std::cout << twiddle[n+r1] << " ";
      std::cout << std::endl;
      const int m = n/radix;

      for (int p=0; p<m; ++p) 
         for (int q=0; q<s; ++q) 
            for (int r1=0; r1<radix; ++r1)
            {
               y[q + s*(radix*p+r1)] = 0.0;
               for (int r2=0; r2<radix; ++r2)
                  y[q + s*(radix*p+r1)] += x[q + s*(p+r2*m)] * twiddle[n + (r1*r2)%radix];
               y[q + s*(radix*p+r1)] *= twiddle[radix*p+r1];
            }
   }
}



void fft_twiddle(int n, complex_t* twiddle, complex_t* x) // Fourier transform
{
    complex_t* y = new complex_t[n];
    int eo = 0;
    int s  = 1;
    int nn = n;
    int nsize = 0;
    while (s<=n) {
       std::cout << "nn=" << nn << " s=" << s << " eo=" << eo  << std::endl;
          
       if ((nn%17)==0)
       {  
          std::cout << " - radix-17 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 17, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 17, twiddle+nsize, x, y);
          nsize += (nn+17);
          nn /= 17; s *= 17; eo = !eo; 
          
          std::cout << " - final nn=" << nn << " s=" << s << std::endl;
       }

       else if ((nn%11)==0)
       {
          std::cout << " - radix-11 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 11, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 11, twiddle+nsize, x, y);
          nsize += (nn+11);
          nn /= 11; s *= 11; eo = !eo;
       }
       else if ((nn%9)==0)
       {
          std::cout << " - radix-9 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 9, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 9, twiddle+nsize, x, y);
          nsize += (nn+9);
          nn /= 9; s *= 9; eo = !eo;
       }
       else if ((nn%8)==0)
       {
          std::cout << " - radix-8 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 8, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 8, twiddle+nsize, x, y);
          nsize += (nn+8);
          nn /= 8; s *= 8; eo = !eo;
       }
       else if ((nn%7)==0)
       {
          std::cout << " - radix-7 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 7, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 7, twiddle+nsize, x, y);
          nsize += (nn+7);
          nn /= 7; s *= 7; eo = !eo;
       }
       else if ((nn%6)==0)
       {
          std::cout << " - radix-6 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 6, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 6, twiddle+nsize, x, y);
          nsize += (nn+6);
          nn /= 6; s *= 6; eo = !eo;
       }
       else if ((nn%5)==0)
       {
          std::cout << " - radix-5 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 5, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 5, twiddle+nsize, x, y);
          nsize += (nn+5);
          nn /= 5; s *= 5; eo = !eo;
       }
       else if ((nn%4)==0)
       {
          std::cout << " - radix-4 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 4, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 4, twiddle+nsize, x, y);
          nsize += (nn+4);
          nn /= 4; s *= 4; eo = !eo;
       }
       else if ((nn%3)==0)
       {
          std::cout << " - radix-3 " << std::endl;
          if (eo)
             fft_radix(nn, s, eo, 3, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 3, twiddle+nsize, x, y);
          nsize += (nn+3);
          nn /= 3; s *= 3; eo = !eo;
       }
       else 
       {
          std::cout << " - radix-2 " << std::endl;
          if (nn==1) nsize = 0;
          if (eo)
             fft_radix(nn, s, eo, 2, twiddle+nsize, y, x);
          else
             fft_radix(nn, s, eo, 2, twiddle+nsize, x, y);
          nsize += (nn+2);
          nn /= 2; s *= 2; eo = !eo;
       }
    }

    delete[] y;
    //for (int k = 0; k < n; k++) x[k] /= n;

}



void fft(int isgn, int n, complex_t* x) // Fourier transform
// isgn: -1 forward fft, 1 inverse fft
// n : sequence length
// x : input/output sequence
{
    complex_t* y = new complex_t[n];
    int eo = 0;
    int s  = 1;
    int nn = n;
    while (s<=n) {
       std::cout << "nn=" << nn << " s=" << s << " eo=" << eo  << std::endl;

       if ((nn%17)==0)
       {
          std::cout << " - radix-17 " << std::endl;
          if (eo)
             fft0_fac17(isgn,nn, s, eo, y, x);
          else
             fft0_fac17(isgn,nn, s, eo, x, y);
          nn /= 17; s *= 17; eo = !eo;

          std::cout << " - final nn=" << nn << " s=" << s << std::endl;
       }

       else if ((nn%11)==0)
       {
          std::cout << " - radix-11 " << std::endl;
          if (eo)
             fft0_fac11(isgn,nn, s, eo, y, x);
          else
             fft0_fac11(isgn,nn, s, eo, x, y);
          nn /= 11; s *= 11; eo = !eo;
       }

       else if ((nn%9)==0)
       {
          std::cout << " - radix-9 " << std::endl;
          if (eo)
             fft0_fac9(isgn,nn, s, eo, y, x);
          else
             fft0_fac9(isgn,nn, s, eo, x, y);
          nn /= 9; s *= 9; eo = !eo;
       }


       else if ((nn%8)==0)
       {
          std::cout << " - radix-8 " << std::endl;
          if (eo)
             fft0_fac8(isgn,nn, s, eo, y, x);
          else
             fft0_fac8(isgn,nn, s, eo, x, y);
          nn /= 8; s *= 8; eo = !eo;
       }
       else if ((nn%7)==0)
       {
          std::cout << " - radix-7 " << std::endl;
          if (eo)
             fft0_fac7(isgn,nn, s, eo, y, x);
          else
             fft0_fac7(isgn,nn, s, eo, x, y);
          nn /= 7; s *= 7; eo = !eo;
       }
       else if ((nn%6)==0)
       {
          std::cout << " - radix-6 " << std::endl;
          if (eo)
             fft0_fac6(isgn,nn, s, eo, y, x);
          else
             fft0_fac6(isgn,nn, s, eo, x, y);
          nn /= 6; s *= 6; eo = !eo;
       }
       else if ((nn%5)==0)
       {
          std::cout << " - radix-5 " << std::endl;
          if (eo)
             fft0_fac5(isgn,nn, s, eo, y, x);
          else
             fft0_fac5(isgn,nn, s, eo, x, y);
          nn /= 5; s *= 5; eo = !eo;
       }
       else if ((nn%4)==0)
       {
          std::cout << " - radix-4 "<< std::endl;
          if (eo)
             fft0_fac4(isgn,nn, s, eo, y, x);
          else
             fft0_fac4(isgn,nn, s, eo, x, y);
          nn /= 4; s *= 4; eo = !eo;
       }
       else if ((nn%3)==0)
       {
          std::cout << " - radix-3 " << std::endl;
          if (eo)
             fft0_fac3(isgn,nn, s, eo, y, x);
          else
             fft0_fac3(isgn,nn, s, eo, x, y);
          nn /= 3; s *= 3; eo = !eo;
       }
       else
       {
          std::cout << " - radix-2 " << std::endl;
          //fft0_fac2(n, 1, 0, x, y);
          if (eo)
             fft0_fac2(isgn,nn, s, eo, y, x);
          else
             fft0_fac2(isgn,nn, s, eo, x, y);
          nn /= 2; s *= 2; eo = !eo;
       }
    }

    delete[] y;
    //for (int k = 0; k < n; k++) x[k] /= n;
}


/*****************************************
 *                                       *
 *            set_sub_fft_twiddle        *
 *                                       *
 *****************************************/
/**
 * @brief Computes twiddle factors for FFT based on the given FFT size and radix.
 *
 * This function calculates two sets of twiddle factors essential for FFT computations:
 * 1. Main twiddle factors: These are used in the FFT computation loops and are
 *    calculated for each combination of group indices and radix elements.
 * 2. Additional radix-specific twiddle factors: These may be used for specific FFT
 *    optimizations or configurations and are calculated separately.
 *
 * The twiddle factors are complex exponential values, which are precomputed to
 * optimize the FFT computation time. The computation involves the cosine and sine
 * functions to create complex exponential terms based on the angles derived from
 * the FFT size 'n' and the specified 'radix'.
 *
 * @param n The total size of the FFT. The number of elements in the main FFT computation.
 *          It is expected that 'n' is divisible by 'radix'.
 * @param radix The base size for the FFT computation groups. This determines how many
 *              groups the FFT data will be divided into and how many additional twiddle
 *              factors will be calculated.
 * @param twiddle A pointer to an array where the computed twiddle factors will be stored.
 *                This array should be of size 'n + radix' to hold all computed factors.
 *                The first 'n' entries will hold the main twiddle factors for the FFT,
 *                followed by 'radix' additional factors for radix-specific computations.
 */
void set_sub_fft_twiddle(const int isgn, const int n, const int radix, complex_t *twiddle)
{
   const int m = n/radix;
   const double theta0 = M_PI/((double) n);
   const double theta_radix = 2*M_PI/((double) radix);

   // first n = (n/radix)*radix = m*radix entries 
   // Calculate the main twiddle factors for the FFT
   for (int p=0; p<m; ++p)                       
      for (int r=0; r<radix; ++r)
         twiddle[r+p*radix] = complex_t(cos(r*p*theta0), isgn*sin(r*p*theta0));

   //next radix entries
   // Calculate additional radix-specific twiddle factors
   for (int r=0; r<radix; ++r)
      twiddle[r+n] = complex_t(cos(r*theta_radix), isgn*sin(r*theta_radix));

}


/*****************************************
 *                                       *
 *            size_fft_twiddle           *
 *                                       *
 *****************************************/
/**
 * @brief Calculates the required storage size for FFT twiddle factors based on factorization.
 *
 * This function iterates through potential radix values and calculates the cumulative
 * storage size required for FFT twiddle factors, considering each level of the factorization
 * process. The function throws an exception if the FFT size cannot be factorized by the given radices.
 *
 * @param n The total size of the FFT.
 * @return int The total storage size required for all twiddle factors.
 * @throws std::runtime_error If the FFT size cannot be completely factorized by the defined radices.
 */
int size_fft_twiddle(const int n)
{
   int nsize = 0;
   int s  = 1;
   int nn = n;
   while (s<=n) {
      if ((nn%17)==0)      { nsize += (nn+17); nn /= 17; s *= 17; }
      else if ((nn%11)==0) { nsize += (nn+11); nn /= 11; s *= 11; }
      else if ((nn%9)==0)  { nsize += (nn+9);  nn /= 9;  s *= 9; }
      else if ((nn%8)==0)  { nsize += (nn+8);  nn /= 8;  s *= 8; }
      else if ((nn%7)==0)  { nsize += (nn+7);  nn /= 7;  s *= 7; }
      else if ((nn%6)==0)  { nsize += (nn+6);  nn /= 6;  s *= 6; }
      else if ((nn%5)==0)  { nsize += (nn+5);  nn /= 5;  s *= 5; }
      else if ((nn%4)==0)  { nsize += (nn+4);  nn /= 4;  s *= 4; }
      else if ((nn%3)==0)  { nsize += (nn+3);  nn /= 3;  s *= 3; }
      else if ((nn%2)==0)  { nsize += (nn+2);  nn /= 2;  s *= 2; }
      else                 {                   nn /= 2;  s *= 2; }
   }
   return nsize;
}


/*****************************************
 *                                       *
 *            set_fft_twiddle            *
 *                                       *
 *****************************************/
void set_fft_twiddle(const int isgn, const int n, complex_t *twiddle)
{
   int nsize = 0;
   int s  = 1;
   int nn = n;
   while (s<=n) {
      if ((nn%17)==0) {
         //const int m = n/17;
         set_sub_fft_twiddle(isgn,nn,17,twiddle+nsize);
         //nsize += 17*m;
         nsize += nn+17;
         nn /= 17; s *= 17;
      }
      else if ((nn%11)==0) { set_sub_fft_twiddle(isgn,nn,11,twiddle+nsize); nsize += (nn+11); nn /= 11; s *= 11; }
      else if ((nn%9)==0)  { set_sub_fft_twiddle(isgn,nn,9, twiddle+nsize); nsize += (nn+9);  nn /= 9;  s *= 9; }
      else if ((nn%8)==0)  { set_sub_fft_twiddle(isgn,nn,8, twiddle+nsize); nsize += (nn+8);  nn /= 8;  s *= 8; }
      else if ((nn%7)==0)  { set_sub_fft_twiddle(isgn,nn,7, twiddle+nsize); nsize += (nn+7);  nn /= 7;  s *= 7; }
      else if ((nn%6)==0)  { set_sub_fft_twiddle(isgn,nn,6, twiddle+nsize); nsize += (nn+6);  nn /= 6;  s *= 6; }
      else if ((nn%5)==0)  { set_sub_fft_twiddle(isgn,nn,5, twiddle+nsize); nsize += (nn+5);  nn /= 5;  s *= 5; }
      else if ((nn%4)==0)  { set_sub_fft_twiddle(isgn,nn,4, twiddle+nsize); nsize += (nn+4);  nn /= 4;  s *= 4; }
      else if ((nn%3)==0)  { set_sub_fft_twiddle(isgn,nn,3, twiddle+nsize); nsize += (nn+3);  nn /= 3;  s *= 3; }
      else if ((nn%2)==0)  { set_sub_fft_twiddle(isgn,nn,2, twiddle+nsize); nsize += (nn+2);  nn /= 2;  s *= 2; }
      else                 {                                                                  nn /= 2;  s *= 2; }
   }
}






/*void ifft(int n, complex_t* x) // Inverse Fourier transform
// n : sequence length
// x : input/output sequence
{
    for (int p = 0; p < n; p++) x[p] = conj(x[p]);
    complex_t* y = new complex_t[n];
    fft0_fac2(n, 1, 0, x, y);
    delete[] y;
    for (int k = 0; k < n; k++) x[k] = conj(x[k]);
}
*/

//#define		N 136
//#define	N 11


//#define	N 8
//#define	N 11
//#define	N 9

//#define		N 81
//#define		N 17
//#define		N 729
//#define		N 121
//#define		N 49
//#define		N 343

//#define		N 2
//#define		N 17

//#define 	Nbig	94248
//#define		N 94248

#define		Nbig 15
#define		N 15

//#define		N 4913
//#define		N 289

int main()
{
   complex_t x[Nbig];
   complex_t x2[Nbig];

   for (int i=0; i<Nbig; ++i)
   {
      if (i>=N)
         x[i] = 0.0;
      else
         x[i] = (i+1)*1.0;
      std::cout << x[i] << " ";
   }
   std::cout << std::endl;

   for (int i=0; i<Nbig; ++i)
      x2[i] = x[i];

  int nsize = size_fft_twiddle(Nbig);
  std::cout << "nsize=" << nsize << std::endl;

  complex_t forward_twiddle[nsize];
  set_fft_twiddle(-1,Nbig,forward_twiddle);

   // forward fft
   fft(-1,Nbig,x);
   // backward fft
   //fft(1,N,x);

   std::cout << "x=" ;
   for (int i=0; i<Nbig; ++i)
      std::cout << x[i] << " ";
   std::cout << std::endl;
   std::cout << std::endl;

   fft_twiddle(Nbig,forward_twiddle,x2);
   std::cout << "x2=" ;
   for (int i=0; i<Nbig; ++i)
      std::cout << x2[i] << " ";
   std::cout << std::endl;

}
