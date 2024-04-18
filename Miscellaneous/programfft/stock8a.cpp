/**
 * @brief Performs a mixed-radix Fast Fourier Transform (FFT) using a radix-based decomposition.
 *
 * The `fft_twiddle` function orchestrates the execution of FFT computations by dynamically 
 * selecting the appropriate radix for each level of decomposition. It leverages a buffer for 
 * intermediate results and toggles between input and output buffers to apply radix-based FFT
 * transformations in a ping-pong fashion. This function requires precomputed twiddle factors 
 * and handles the FFT computation entirely in-place on the input data array.
 *
 * @param n The total number of elements in the input array `x`. The value should be decomposable 
 *          by the radices supported by the function (from the set {2, 3, 4, ..., 17}).
 * @param twiddle A pointer to an array of precomputed twiddle factors, which must be appropriately 
 *                sized and initialized before calling this function. These factors are crucial 
 *                for the radix-specific transformations within the FFT computation.
 * @param x Pointer to the complex number array representing the input data. This array is modified 
 *          in-place to contain the output Fourier coefficients upon completion of the FFT.
 *
 * @note The function assumes that memory management (allocation and deallocation) for `x` and `twiddle`
 *       is handled externally. It internally allocates and deallocates the secondary buffer `y`, 
 *       used for intermediate FFT results.
 */

#define		Nbig 256
#define		N 256

#include 	<cmath>
#include 	<iostream>
#include 	<complex>
#include 	<stdexcept>  // Required for std::runtime_error
#include	<initializer_list> 

typedef std::complex<double> complex_t;

// Define a constants for the radix values
constexpr int radix_values[] = {32, 17, 16, 11, 9, 8, 7, 6, 5, 4, 3, 2};



/*****************************************
 *                                       *
 *            fft_radix                  *
 *                                       *
 *****************************************/
/**
 * @brief Performs a Radix-based Fast Fourier Transform (FFT) on complex input data.
 *
 * This function computes the FFT using a radix-based approach, where the input data
 * is decomposed into smaller chunks using the specified radix. The transformation
 * is applied recursively until base cases are reached. The function also handles
 * the twiddle factor multiplications specific to FFT computations.
 *
 * @param n The total number of elements in the input array `x`. Must be a power of `radix`.
 * @param s The stride between successive elements to be transformed, allowing in-place transformation
 *          of subarrays within a larger dataset.
 * @param eo A boolean flag indicating whether the transformation is even or odd, which can
 *          influence certain optimization paths or conditional transformations.
 * @param radix The radix to be used for the FFT computation, typically a small power of 2 (e.g., 2, 4, 8).
 *              This parameter dictates the size of the decomposition at each recursive step.
 * @param twiddle Pointer to an array containing precomputed twiddle factors used in the FFT computations.
 * @param x Pointer to the input array of complex numbers to be transformed. The array size should be at least `n * s`.
 * @param y Pointer to the output array where the transformed data will be stored. Must have the same size as `x`.
 *
 * @note The function modifies the array `y` directly, and it is the caller's responsibility to ensure
 *       that `x` and `y` are appropriately sized and allocated before calling this function.
 */
void fft_radix(const int n, const int s, bool eo, const int radix, complex_t* twiddle, complex_t* x, complex_t *y)
{
   if (n == 1) { if (eo) for (int q = 0; q < s; q++) y[q] = x[q]; }
   else 
   {
      const int m = n/radix;
      complex_t Atwiddle[radix*radix];

      // Precompute twiddle factors for matrix multiplication
      for (int r2=0; r2<radix; ++r2)
      for (int r1=0; r1<radix; ++r1)
         Atwiddle[r1+r2*radix] = twiddle[(r1*r2)%radix];

      for (int p=0; p<m; ++p) 
      {
         for (int q=0; q<s; ++q) 
         {
            complex_t* x1 = x+q+s*p;
            complex_t* y1 = y+q+s*radix*p;

            // Initialize the output y1 vector
            for (int r1 = 0; r1 < radix; ++r1) 
               y1[s * r1] = 0.0;

            // Matrix-vector multiplication
            for (int r2 = 0; r2 < radix; ++r2)
            for (int r1 = 0; r1 < radix; ++r1) 
                    y1[s*r1] += x1[s*m*r2] * Atwiddle[r1 + r2 * radix];

            // Apply phase factor to each result
            for (int r1 = 0; r1 < radix; ++r1) 
            {
               complex_t fac = complex_t(1.0, 0.0);
               for (int ps = 0; ps < p; ++ps)
                  fac *= twiddle[radix + r1];

               y1[s * r1] *= fac;
            }
         }
      }
   }
}


/*****************************************
 *                                       *
 *            fft_twiddle                *
 *                                       *
 *****************************************/
/**
 * @brief Performs a mixed-radix Fast Fourier Transform (FFT) on a complex array.
 *
 * This function implements an FFT using a mixed-radix algorithm, which supports multiple
 * radix values (2, 3, 4, ..., 17). It utilizes the `fft_radix` function to handle each 
 * specific radix. The function decomposes the input size recursively based on available 
 * radix factors until all elements are processed. The FFT processes the array in-place 
 * with the help of a secondary buffer `y` for intermediate results.
 *
 * @param n The total number of elements in the input array `x`. It should ideally be
 *          a composite number that can be factored into the radix values supported by
 *          this function (2 through 17).
 * @param twiddle Pointer to an array containing precomputed twiddle factors specific to
 *                each radix operation. This array should be sized appropriately to 
 *                accommodate all necessary twiddle factors for the given `n`.
 * @param x Pointer to the input/output array of complex numbers. The array is modified
 *          in-place to contain the Fourier coefficients after transformation. The array
 *          size must be at least `n`.
 *
 * @note The function dynamically allocates memory for the secondary buffer `y`, which is
 *       used during the computation and freed before exiting. Users should ensure that
 *       the input array `x` and the twiddle factor array `twiddle` are allocated and
 *       populated before calling this function. Post-processing (like normalizing the
 *       FFT output) is not performed inside this function and should be handled externally
 *       if required.
 */
void fft_twiddle(int n, complex_t* twiddle, complex_t* x) // Fourier transform
{
   complex_t* y = new complex_t[n];
   int eo = 0;
   int s  = 1;
   int nn = n;
   int nsize = 0;
   while (s<=n) 
   {
      std::cout << "nn=" << nn << " s=" << s << " eo=" << eo  << std::endl;

      // Identify the largest radix applicable for current nn
      int radix = 2;  // Default to radix-2
      for (int r : radix_values) {
         if (nn % r == 0) {
            radix = r;
            break;
         }
      }

      std::cout << " - radix-" << radix << std::endl;

      // Perform FFT with the determined radix
      if (eo)
         fft_radix(nn, s, eo, radix, twiddle + nsize, y, x);
      else
         fft_radix(nn, s, eo, radix, twiddle + nsize, x, y);

      nsize += 2*radix;
      nn /= radix;
      s *= radix;
      eo = !eo;  // Toggle the 'even-odd' flag

      std::cout << " - final nn=" << nn << " s=" << s << std::endl;
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
   const double theta0 = 2*M_PI/((double) n);
   const double theta_radix = 2*M_PI/((double) radix);

   // Calculate radix-specific twiddle factors
   for (int r=0; r<radix; ++r)
      twiddle[r] = complex_t(cos(r*theta_radix), isgn*sin(r*theta_radix));

   // Calculate the main twiddle factors for the FFT
   for (int r=0; r<radix; ++r)
      twiddle[radix+r] = complex_t(cos(r*theta0), isgn*sin(r*theta0));

}


/*****************************************
 *                                       *
 *            size_fft_twiddle           *
 *                                       *
 *****************************************/
/**
 * @brief Calculate the total storage size required for FFT twiddle factors.
 *
 * This function iterates through the provided radix values to decompose the FFT size (n)
 * into its factors. For each factor found, it accumulates the storage size required for 
 * the twiddle factors of the FFT, considering both real and imaginary parts. The function 
 * ensures that the entire FFT size can be factorized using the given radices and throws an 
 * exception if any part of the size remains unfactorized.
 *
 * @param n The total size of the FFT for which twiddle factors are required.
 * @param radix_values A vector of integers representing the radix values that can be used 
 *        to factorize the FFT size.
 * @return int The total storage size required for the twiddle factors, accounting for 
 *         all levels of factorization.
 * @throws std::runtime_error If the FFT size cannot be completely factorized by the provided 
 *         radices.
 */
int size_fft_twiddle(const int n) 
{
   int nsize = 0;
   int s = 1;
   int nn = n;

   while (s <= n) 
   {
      bool found = false;
      // Loop through possible radix values to find the largest factor of nn
      for (int radix : radix_values) 
      {
         if (nn % radix == 0) 
         {
            nsize += 2 * radix;
            nn /= radix;
            s *= radix;
            found = true;
            break;
         }
      }
      if (!found) break;
   }

   return nsize;
}


/*****************************************
 *                                       *
 *            set_fft_twiddle            *
 *                                       *
 *****************************************/
/**
 * Initializes the FFT twiddle factors for different radix levels based on the input size.
 * 
 * @param isgn The sign indicator for the FFT (usually +1 or -1).
 * @param n The size of the input to the FFT. Must be decomposable by the radices used within.
 * @param twiddle Pointer to the array where the twiddle factors are to be stored.
 */
void set_fft_twiddle(const int isgn, const int n, complex_t *twiddle) 
{
   int nsize = 0;
   int s = 1;
   int nn = n;

   while (s <= n) 
   {
      bool found = false;
      // Loop through possible radix values to find the largest factor of nn
      for (int radix : radix_values) 
      {
         if (nn % radix == 0) 
         {
            set_sub_fft_twiddle(isgn, nn, radix, twiddle + nsize);
            nsize += 2 * radix;
            nn /= radix;
            s *= radix;
            found = true;
            break;
         }
      }
      if (!found) break;
   }
}


int main()
{
   complex_t x[Nbig];

   for (int i=0; i<Nbig; ++i)
   {
      if (i>=N)
         x[i] = 0.0;
      else
         x[i] = (i+1)*1.0;
      std::cout << x[i] << " ";
   }
   std::cout << std::endl;

   int nsize = size_fft_twiddle(Nbig);
   std::cout << "nsize=" << nsize << std::endl;

   complex_t forward_twiddle[nsize];
   complex_t backward_twiddle[nsize];
   set_fft_twiddle(-1,Nbig,forward_twiddle);
   set_fft_twiddle(1,Nbig,backward_twiddle);

   // forward fft
   //fft(-1,Nbig,x);
   // backward fft
   //fft(1,N,x);

   fft_twiddle(Nbig,forward_twiddle,x);
   std::cout << "x =" ;
   for (int i=0; i<Nbig; ++i)
      std::cout << x[i] << " ";
   std::cout << std::endl;
   std::cout << std::endl;

   fft_twiddle(Nbig,backward_twiddle,x);
   std::cout << "x2 =" ;
   for (int i=0; i<Nbig; ++i)
      std::cout << x[i]/((double) Nbig) << " ";
   std::cout << std::endl;
}
