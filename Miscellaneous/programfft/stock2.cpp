#include <complex>
#include <cmath>

typedef std::complex<double> complex_t;

void fft1(int n, int s, complex_t* x, complex_t* y);

void fft0(int n, int s, complex_t* x, complex_t* y)
// n : sequence length
// s : stride
// x : input/output sequence
// y : work area
{
    const int m = n/2;
    const double theta0 = 2*M_PI/n;

    if (n == 1) {}
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
        fft1(n/2, 2*s, y, x);
    }
}

void fft1(int n, int s, complex_t* x, complex_t* y)
// n : sequence length
// s : stride
// x : input sequence
// y : output sequence
{
    const int m = n/2;
    const double theta0 = 2*M_PI/n;

    if (n == 1) { for (int q = 0; q < s; q++) y[q] = x[q]; }
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
        fft0(n/2, 2*s, y, x);
    }
}

void fft(int n, complex_t* x) // Fourier transform
// n : sequence length
// x : input/output sequence
{
    complex_t* y = new complex_t[n];
    fft0(n, 1, x, y);
    delete[] y;
    for (int k = 0; k < n; k++) x[k] /= n;
}

void ifft(int n, complex_t* x) // Inverse Fourier transform
// n : sequence length
// x : input/output sequence
{
    for (int p = 0; p < n; p++) x[p] = conj(x[p]);
    complex_t* y = new complex_t[n];
    fft0(n, 1, x, y);
    delete[] y;
    for (int k = 0; k < n; k++) x[k] = conj(x[k]);
}
