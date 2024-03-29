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
    //std::cout << "n=" << n << " s=" << s << " eo=" << eo  << std::endl;
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

       if ((nn%8)==0)
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

#define	N 36
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
   fft(-1,N,x);

   std::cout << std::endl << std::endl;
   // backward fft
   //fft(1,N,x);

   for (int i=0; i<N; ++i)
      std::cout << x[i] << " ";
   std::cout << std::endl;

}
