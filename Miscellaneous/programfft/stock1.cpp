#include 	<complex>
#include 	<cmath>
#include	<iostream>

typedef std::complex<double> complex_t;

void fft0_fac2(int isgn, int n, int s, bool eo, complex_t* x, complex_t* y)
// isgn: -1 forward fft, 1 inverse fft
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    std::cout << "n=" << n << " s=" << s << " eo=" << eo  << std::endl;
    const int m = n/2;
    const double theta0 = 2*M_PI/n;

    if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
    else {
        for (int p = 0; p < m; p++) {
            const complex_t wp = complex_t(cos(p*theta0), isgn*sin(p*theta0));
            for (int q = 0; q < s; q++) {
                const complex_t a = x[q + s*(p + 0)];
                const complex_t b = x[q + s*(p + m)];
                y[q + s*(2*p + 0)] =  a + b;
                y[q + s*(2*p + 1)] = (a - b) * wp;
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

       if ((nn%7)==0)
       {
          if (eo)
             fft0_fac7(isgn,nn, s, eo, y, x);
          else
             fft0_fac7(isgn,nn, s, eo, x, y);
          nn /= 7; s *= 7; eo = !eo;
       }
       else if ((nn%5)==0)
       {
          if (eo)
             fft0_fac5(isgn,nn, s, eo, y, x);
          else
             fft0_fac5(isgn,nn, s, eo, x, y);
          nn /= 5; s *= 5; eo = !eo;
       }
       else if ((nn%3)==0)
       {
          if (eo)
             fft0_fac3(isgn,nn, s, eo, y, x);
          else
             fft0_fac3(isgn,nn, s, eo, x, y);
          nn /= 3; s *= 3; eo = !eo;
       }
       else
       {
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

#define	N 21
int main()
{
   complex_t x[N];

   for (int i=0; i<N; ++i)
   {
      x[i] = (i+1)*1.0;
      std::cout << x[i] << " ";
   }
   std::cout << std::endl;

   // forward fft
   //fft(-1,N,x);

   // backward fft
   fft(1,N,x);

   for (int i=0; i<N; ++i)
      std::cout << x[i] << " ";
   std::cout << std::endl;

}
