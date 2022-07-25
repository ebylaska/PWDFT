#include 	<complex>
#include 	<cmath>
#include	<iostream>

typedef std::complex<double> complex_t;

void fft0_fac2(int n, int s, bool eo, complex_t* x, complex_t* y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    std::cout << "n=" << n << " s=" << s << " eo=" << eo << std::endl;
    const int m = n/2;
    const double theta0 = 2*M_PI/n;

    if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
    else {
        for (int p = 0; p < m; p++) {
            const complex_t wp = complex_t(cos(p*theta0), -sin(p*theta0));
            for (int q = 0; q < s; q++) {
                const complex_t a = x[q + s*(p + 0)];
                const complex_t b = x[q + s*(p + m)];
                y[q + s*(2*p + 0)] =  a + b;
                y[q + s*(2*p + 1)] = (a - b) * wp;
            }
        }
        fft0_fac2(n/2, 2*s, !eo, y, x);
    }
}

void fft(int n, complex_t* x) // Fourier transform
// n : sequence length
// x : input/output sequence
{
    complex_t* y = new complex_t[n];
    fft0_fac2(n, 1, 0, x, y);
    delete[] y;
    //for (int k = 0; k < n; k++) x[k] /= n;
}

void ifft(int n, complex_t* x) // Inverse Fourier transform
// n : sequence length
// x : input/output sequence
{
    for (int p = 0; p < n; p++) x[p] = conj(x[p]);
    complex_t* y = new complex_t[n];
    fft0_fac2(n, 1, 0, x, y);
    delete[] y;
    for (int k = 0; k < n; k++) x[k] = conj(x[k]);
}

#define	N 64
int main()
{
   complex_t x[N];

   double sum0=0.0;
   for (int i=0; i<N; ++i)
   {
      x[i] = (i+1)*1.0;
      sum0 +=(i+1)*1.0;
   }
   std::cout << "sum0=" << sum0 << std::endl;

   fft(N,x);

   for (int i=0; i<N; ++i)
      std::cout << x[i] << " ";
   std::cout << std::endl;

}
