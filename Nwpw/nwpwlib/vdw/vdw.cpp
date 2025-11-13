


#include <cmath>

#include "Parallel.hpp"
#include "util.hpp"
#include "vdw.hpp"
#include "Control2.hpp"
#include "compressed_io.hpp"

#include <fstream>
#include <iostream>
#include "iofmt.hpp"
#include <cstring>



namespace pwdft {


/**********************************************
 *                                            *
 *        vdw_DF_kernel_phi_value             *
 *                                            *
 **********************************************
 *
 *  Purpose:
 *  --------
 *  Evaluate the dimensionless van der Waals kernel function
 *  φ(d₁, d₂) for a given pair of scaled distances d₁, d₂
 *  within the nonlocal correlation functional (vdW-DF).
 *
 *  This routine computes a single value of the kernel φ by
 *  performing the double summation over the tabulated angular
 *  integration grid {aᵢ}, using precomputed coupling weights
 *  W_ab(i,j) and dynamically updated “ν” parameters ν(aᵢ,d₁)
 *  and ν₁(aᵢ,d₂).
 *
 *  Mathematical Form:
 *  ------------------
 *       φ(d₁,d₂) = (1/π²) Σᵢ Σⱼ  W_ab(i,j) ·
 *                    [ (1/(νᵢ+νⱼ) + 1/(ν₁ᵢ+ν₁ⱼ)) ×
 *                      ( 1/((νᵢ+ν₁ᵢ)(νⱼ+ν₁ⱼ))
 *                      + 1/((νᵢ+ν₁ⱼ)(ν₁ᵢ+νⱼ)) ) ]
 *
 *  where ν and ν₁ are smooth functions of a and d:
 *
 *       ν(a,d)  = a² / [ 2(1 − exp(−a²γ / d²)) ]
 *       γ       = 4π/9
 *
 *  and asymptotic limits are handled analytically to prevent
 *  numerical overflow when d → 0 or a²γ/d² → ∞.
 *
 *  Parameters:
 *  -----------
 *    Na        : number of angular grid points
 *    a[Na]     : transformed abscissas  (tan θ)
 *    a2[Na]    : squares of abscissas   (aᵢ²)
 *    nu[Na]    : working array for ν(aᵢ,d₁)   (output/temporary)
 *    nu1[Na]   : working array for ν(aᵢ,d₂)   (output/temporary)
 *    Wab[]     : flattened row-major Na×Na matrix of coupling weights
 *    d1, d2    : scaled dimensionless distances
 *
 *  Returns:
 *  --------
 *    φ(d₁,d₂)  : kernel value suitable for real-space integration
 *                over r = |r₁ − r₂|.
 *
 *  Notes:
 *  ------
 *  • For d₁ = d₂ = 0 the kernel is set to zero.
 *  • A small-value cutoff (1×10⁻¹²) prevents division by zero.
 *  • Exponentials > exp(700) are clamped to avoid overflow.
 *  • The arrays `nu` and `nu1` are reused to minimize allocations.
 *  • This function is typically called inside
 *        vdw_DF_kernel_gen_phir()
 *    to populate φ_r(r) on a radial grid.
 *
 *  References:
 *  -----------
 *    M. Dion et al., *Phys. Rev. Lett.* **92**, 246401 (2004)
 *    K. Lee et al., *Phys. Rev. B* **82**, 081101(R) (2010)
 *    J. Klimeš, D. R. Bowler, A. Michaelides,
 *       *J. Phys.: Condens. Matter* **22**, 022201 (2010)
 *
 **********************************************/
double vdw_DF_kernel_phi_value(int Na,
                               const double* a,
                               const double* a2,
                               double* nu,
                               double* nu1,
                               const double* Wab,  // flattened Na×Na array, row-major
                               double d1,
                               double d2)
{
    const double small = 1.0e-12;
    const double pi = 4.0 * std::atan(1.0);
    const double gamma = 4.0 * pi / 9.0;

    double d1s = d1 * d1;
    double d2s = d2 * d2;

    double phi_value = 0.0;

    // Handle trivial case
    if (d1 == 0.0 && d2 == 0.0)
        return 0.0;

    // Compute nu(i) and nu1(i)
    for (int i = 0; i < Na; ++i)
    {
        // ---- nu(i)
        if (a[i] <= small && d1 > small)
            nu[i] = (9.0 / 8.0) * d1s / pi;
        else if (d1 <= small || (a2[i] * gamma / d1s) > 700.0)
            nu[i] = a2[i] / 2.0;
        else
            nu[i] = a2[i] / ((1.0 - std::exp(-(a2[i] * gamma) / d1s)) * 2.0);

        // ---- nu1(i)
        if (a[i] <= small && d2 > small)
            nu1[i] = (9.0 / 8.0) * d2s / pi;
        else if (d2 <= small || (a2[i] * gamma / d2s) > 700.0)
            nu1[i] = a2[i] / 2.0;
        else
            nu1[i] = a2[i] / ((1.0 - std::exp(-(a2[i] * gamma) / d2s)) * 2.0);
    }

    // Double summation over i,j
    for (int i = 0; i < Na; ++i)
    {
        for (int j = 0; j < Na; ++j)
        {
            double w = nu[i];
            double x = nu[j];
            double y = nu1[i];
            double z = nu1[j];

            if (w > small && x > small && y > small && z > small)
            {
                double T = (1.0 / (w + x) + 1.0 / (y + z)) *
                           (1.0 / ((w + y) * (x + z)) +
                            1.0 / ((w + z) * (y + x)));

                phi_value += T * Wab[i * Na + j];  // row-major indexing
            }
        }
    }

    phi_value /= (pi * pi);
    return phi_value;
}



/**********************************************
 *                                            *
 *         vdw_DF_kernel_gen_phir             *
 *                                            *
 **********************************************
 *
 *  Purpose:
 *  --------
 *  Generate the real-space radial kernel function φ_r(r)
 *  for the van der Waals density functional (vdW-DF)
 *  nonlocal correlation term, evaluated on a uniform
 *  radial grid.
 *
 *  This routine forms φ_r(r) by integrating over
 *  precomputed tables (W_ab, ν, ν₁, etc.) for each
 *  pair of q₁ and q₂ parameters describing the
 *  local response functions.
 *
 *  Mathematical Context:
 *  ---------------------
 *  The vdW-DF kernel φ depends on the reduced response
 *  parameters q₁, q₂ and distance r through the
 *  symmetrized dimensionless variables:
 *
 *       d₁ = r (1 + δq)
 *       d₂ = r (1 − δq)
 *       δq = (q₁ − q₂) / (q₁ + q₂)
 *
 *  For each grid point rᵢ, φ_r(rᵢ) is obtained by:
 *
 *       φ_r(rᵢ) = φ_value(Na, a, a², ν, ν₁, W_ab, d₁, d₂)
 *
 *  where `vdw_DF_kernel_phi_value` interpolates and
 *  evaluates the full tabulated kernel.
 *
 *  Parameters:
 *  -----------
 *    Na        : number of tabulated angular points (a-grid)
 *    a[Na]     : a_i = tan(θ_i) abscissas
 *    a2[Na]    : squares of a_i
 *    nu[Na]    : precomputed kernel parameter array ν(a_i)
 *    nu1[Na]   : precomputed kernel derivative array ν₁(a_i)
 *    Wab[]     : flattened (row-major) Na×Na coupling matrix
 *    q1, q2    : response parameters at r₁, r₂
 *    nr, dr    : number of radial points and radial spacing
 *    phir[]    : output array of φ_r(rᵢ), size nr+1
 *
 *  Notes:
 *  ------
 *  • The function assumes a uniform radial grid: rᵢ = i·dr.
 *  • Fortran convention i=1..nr is extended to 0..nr for safety.
 *  • Calls external routine:
 *        vdw_DF_kernel_phi_value(...)
 *    which performs the core evaluation of the kernel function.
 *  • δq controls the asymmetry between q₁ and q₂; δq = 0
 *    corresponds to the symmetric case q₁ = q₂.
 *
 *  References:
 *  -----------
 *    M. Dion et al., Phys. Rev. Lett. 92, 246401 (2004)
 *    K. Lee et al., Phys. Rev. B 82, 081101(R) (2010)
 *    J. Klimeš et al., J. Phys.: Condens. Matter 22, 022201 (2010)
 *
 **********************************************/
void vdw_DF_kernel_gen_phir(int Na,
                            const double* a,
                            const double* a2,
                            double* nu,
                            double* nu1,
                            const double* Wab,  // flattened Na×Na array (row-major)
                            double q1,
                            double q2,
                            int nr,
                            double dr,
                            double* phir)      // output: size nr+1
{
    double qdelta = (q1 - q2) / (q1 + q2);

    for (int i = 0; i <= nr; ++i)  // Fortran loop: i=1,nr -> includes 0..nr for safety
    {
        double r  = i * dr;
        double d1 = r * (1.0 + qdelta);
        double d2 = r * (1.0 - qdelta);

        phir[i] = vdw_DF_kernel_phi_value(Na, a, a2, nu, nu1, Wab, d1, d2);
    }
}




/**********************************************
 *                                            *
 *            vdw_GaussLegendre               *
 *                                            *
 **********************************************
 *
 *  Purpose:
 *  --------
 *  Generate Gauss–Legendre quadrature abscissas and weights
 *  on an arbitrary interval [amin, amax] for use in the
 *  vdW-DF kernel integration routines.
 *
 *  The points {a_i} and weights {w_i} satisfy:
 *
 *       ∫ₐₘᵢₙ^ₐₘₐₓ f(x) dx ≈ Σ w_i f(a_i)
 *
 *  where the abscissas are the roots of the Legendre polynomial
 *  P_N(x) of order Npoints, and the weights are obtained from
 *  the standard Legendre–Gauss formula.
 *
 *  Method:
 *  -------
 *  1. Use Chebyshev-based initial guesses for P_N(x) roots:
 *         x₀ ≈ cos[π (i − ¼) / (N + ½)]
 *  2. Refine each root using Newton’s method applied to P_N(x).
 *  3. Compute the corresponding quadrature weight:
 *         w_i = 2 / [(1 − x_i²) (P′_N(x_i))²]
 *  4. Map {x_i} and {w_i} from [−1,1] to [amin, amax].
 *
 *  Arguments:
 *  ----------
 *    amin, amax : integration limits
 *    Npoints    : number of quadrature points
 *    a[]        : output array of mapped abscissas (size Npoints)
 *    weights[]  : output array of mapped weights  (size Npoints)
 *
 *  Notes:
 *  ------
 *  • Symmetry is used: only half of the roots are computed and
 *    mirrored for efficiency.
 *  • Convergence tolerance is 1×10⁻¹⁴ for Newton refinement.
 *  • This routine provides the numerical backbone for constructing
 *    integration grids in routines such as `vdw_DF_kernel_gen_Wab`.
 *
 *  References:
 *  -----------
 *    W. H. Press et al., *Numerical Recipes in C*, 2nd Ed.,
 *        Cambridge Univ. Press (1992), §4.5 “Gauss–Legendre Quadrature.”
 *    M. Dion et al., *Phys. Rev. Lett.* 92, 246401 (2004)
 *    K. Lee et al., *Phys. Rev. B* 82, 081101(R) (2010)
 *
 **********************************************/
void vdw_GaussLegendre(double amin, double amax, int Npoints,
                   double *a, 
                   double *weights)
{
    const double pi = 4.0 * std::atan(1.0);

    //a.resize(Npoints);
    //weights.resize(Npoints);

    int N = (Npoints + 1) / 2;
    double midpoint = 0.5 * (amin + amax);
    double length   = 0.5 * (amax - amin);

    for (int i = 1; i <= N; ++i)
    {
        // Initial guess for root using the Chebyshev nodes
        double root = std::cos(pi * (i - 0.25) / (Npoints + 0.5));
        double last_root, poly1, poly2, poly3, dpdx;

        bool done = false;

        // Newton iteration to refine root
        while (!done)
        {
            poly1 = 1.0;
            poly2 = 0.0;
            for (int j = 1; j <= Npoints; ++j)
            {
                poly3 = poly2;
                poly2 = poly1;
                poly1 = ((2.0 * j - 1.0) * root * poly2 - (j - 1.0) * poly3) / j;
            }
            dpdx = Npoints * (root * poly1 - poly2) / (root * root - 1.0);

            last_root = root;
            root = last_root - poly1 / dpdx;

            if (std::fabs(root - last_root) <= 1.0e-14)
                done = true;
        }

        // Compute points mapped to [amin, amax]
        a[i - 1]             = midpoint - length * root;
        a[Npoints - i]       = midpoint + length * root;

        // Compute weights
        double w = 2.0 * length / ((1.0 - root * root) * dpdx * dpdx);
        weights[i - 1]       = w;
        weights[Npoints - i] = w;
    }
}


/**********************************************
 *                                            *
 *               vdw_DF_Fsin                  *
 *                                            *
 **********************************************
 *
 *  Purpose:
 *  --------
 *  Evaluate the sine–weighted integral of a cubic spline segment
 *  over an interval [x₀, x₁] for use in Filon-type quadrature of
 *  oscillatory kernels in the vdW-DF method.
 *
 *  Specifically, this routine computes the contribution of one
 *  cubic spline segment y(x) between x₀ and x₁ to an integral of
 *  the form:
 *
 *        I = ∫ₓ₀ˣ₁ y(x) · sin(gx) dx
 *
 *  where y(x) is represented by its values and second derivatives
 *  at the endpoints (y₀, y₁, y″₀, y″₁).  Analytic expressions for
 *  ∫ φ(x)sin(gx)dx of cubic basis functions are combined to yield
 *  coefficients A–D such that:
 *
 *        I = A·y₀ + B·y₁ + C·y″₀ + D·y″₁
 *
 *  Parameters:
 *  -----------
 *    x0, x1  : interval endpoints
 *    y0, y1  : spline function values at x₀, x₁
 *    ypp0,
 *    ypp1    : second derivatives y″(x₀), y″(x₁)
 *    g       : oscillation frequency (sin argument multiplier)
 *
 *  Returns:
 *  --------
 *    The integral value  ∫ₓ₀ˣ₁ y(x) sin(gx) dx
 *    evaluated analytically via cubic-spline coefficients.
 *
 *  Method:
 *  -------
 *    1. Construct analytic coefficients (A,B,C,D) for each term
 *       in the cubic polynomial representation of y(x).
 *    2. Combine with boundary values and second derivatives.
 *
 *  Usage:
 *  ------
 *    Called within Filon-type integrators such as
 *    `vdw_DF_kernel_bessel()` to efficiently integrate spline-
 *    interpolated functions multiplied by sin(ωx).
 *
 *  Notes:
 *  ------
 *    • Avoid calling with g = 0 (singular denominators).
 *    • x₀ ≠ x₁ required.
 *    • Equivalent to Fortran routine `vdw_DF_Fsin` from the
 *      original NWPW/Quantum vdW-DF implementation.
 *
 *  References:
 *  -----------
 *    M. Dion et al., Phys. Rev. Lett. 92, 246401 (2004)
 *    K. Lee et al., Phys. Rev. B 82, 081101(R) (2010)
 *    P. M. Morse & H. Feshbach, *Methods of Theoretical Physics*,
 *    McGraw–Hill (1953), Sec. 7.3 (Filon integration methods)
 *
 **********************************************/
double vdw_DF_Fsin(double x0, double x1,
                   double y0, double y1,
                   double ypp0, double ypp1,
                   double g)
{
    double Asin0, Bsin0, Csin0, Dsin0;

    Asin0 = (g * (x0 - x1) * std::cos(g * x0)
            - std::sin(g * x0)
            + std::sin(g * x1))
            / (std::pow(g, 2.0) * (x0 - x1));

    Bsin0 = (-(g * (x0 - x1) * std::cos(g * x1))
            + std::sin(g * x0)
            - std::sin(g * x1))
            / (std::pow(g, 2.0) * (x0 - x1));

    Csin0 = -(6.0 * g * (x0 - x1) * std::cos(g * x0)
            + 2.0 * (-3.0 + std::pow(g, 2.0) * std::pow(x0 - x1, 2.0)) * std::sin(g * x0)
            + (6.0 + std::pow(g, 2.0) * std::pow(x0 - x1, 2.0)) * std::sin(g * x1))
            / (6.0 * std::pow(g, 4.0) * (x0 - x1));

    Dsin0 = (-6.0 * std::sin(g * x0)
            + 6.0 * std::sin(g * x1)
            + g * (x0 - x1)
              * (6.0 * std::cos(g * x1)
              - g * (x0 - x1)
                * (std::sin(g * x0) + 2.0 * std::sin(g * x1))))
            / (6.0 * std::pow(g, 4.0) * (x0 - x1));

    double vdw_DF_Fsin = Asin0 * y0 + Bsin0 * y1 + Csin0 * ypp0 + Dsin0 * ypp1;
    return vdw_DF_Fsin;
}



/**********************************************
 *                                            *
 *          vdw_DF_kernel_gen_Wab             *
 *                                            *
 **********************************************
 *
 *  Purpose:
 *  --------
 *  Generate the W_ab(Na×Na) coupling matrix and associated
 *  integration arrays for the van der Waals density functional
 *  (vdW-DF) nonlocal correlation kernel.
 *
 *  The matrix elements W_ab correspond to the angular part of
 *  the double integral appearing in the kernel:
 *
 *        E_c^nl = 1/2 ∬ n(r) φ(|r−r'|, q₁, q₂) n(r') dr dr'
 *
 *  where φ depends on pretabulated functions of parameters
 *  (a_i, a_j) defined via the variable transformation
 *  a = tan(θ) over θ ∈ [atan(0), atan(64)].
 *
 *  Method:
 *  -------
 *  1. Obtain Gauss–Legendre quadrature points and weights for
 *     θ ∈ [atan(0), atan(64)] via `vdw_GaussLegendre`.
 *  2. Transform integration variable from θ → a = tan(θ);
 *     apply Jacobian factor (1 + a²) to weights.
 *  3. Precompute trigonometric factors sin(a), cos(a).
 *  4. Construct symmetric matrix W_ab:
 *
 *        W_ab(i,j) = 2 w_i w_j
 *           [ (3−a_i²)a_j cos a_j sin a_i
 *           + (3−a_j²)a_i cos a_i sin a_j
 *           + (a_i²+a_j²−3) sin a_i sin a_j
 *           − 3 a_i a_j cos a_i cos a_j ] / (a_i a_j)
 *
 *  Arguments:
 *  ----------
 *    Na         : number of quadrature points
 *    a[Na]      : transformed abscissas (tan θ)
 *    a2[Na]     : squares of a (a_i²)
 *    aweights[] : modified Gauss–Legendre weights (w_i (1+a_i²))
 *    cos_a[]    : cos(a_i) precomputed
 *    sin_a[]    : sin(a_i) precomputed
 *    Wab[]      : output coupling matrix, flattened row-major (Na×Na)
 *
 *  Notes:
 *  ------
 *  • Integration limits correspond to a ∈ [0, 64].
 *  • The output W_ab matrix is dimensionless and symmetric.
 *  • The prefactor 4π is prepared for subsequent normalization
 *    in kernel tabulation routines.
 *  • Requires external function:
 *        vdw_GaussLegendre(...) – Gauss–Legendre node/weight generator.
 *
 *  References:
 *  -----------
 *    M. Dion, H. Rydberg, E. Schröder, D. C. Langreth, and B. I. Lundqvist,
 *        “van der Waals density functional for general geometries,”
 *        Phys. Rev. Lett. 92, 246401 (2004)
 *    K. Lee, E. D. Murray, L. Kong, B. I. Lundqvist, and D. C. Langreth,
 *        “Higher-accuracy van der Waals density functional,”
 *        Phys. Rev. B 82, 081101(R) (2010)
 *
 **********************************************/
void vdw_DF_kernel_gen_Wab(int Na,
                           double* a,
                           double* a2,
                           double* aweights,
                           double* cos_a,
                           double* sin_a,
                           double* Wab)
{
    const double fourpi = 16.0 * std::atan(1.0);  // 4π
    const double amin = 0.0;
    const double amax = 64.0;

    // --- Step 1: Gauss-Legendre integration setup
    vdw_GaussLegendre(std::atan(amin), std::atan(amax), Na, a, aweights);

    // --- Step 2: Transform integration variables and compute trigonometric arrays
    for (int i = 0; i < Na; ++i)
    {
        a[i] = std::tan(a[i]);
        a2[i] = a[i] * a[i];
        aweights[i] = aweights[i] * (1.0 + a2[i]);
        cos_a[i] = std::cos(a[i]);
        sin_a[i] = std::sin(a[i]);
    }

    // --- Step 3: Construct Wab matrix (Na×Na)
    for (int i = 0; i < Na; ++i)
    {
        for (int j = 0; j < Na; ++j)
        {
            double numerator =
                (3.0 - a2[i]) * a[j] * cos_a[j] * sin_a[i] +
                (3.0 - a2[j]) * a[i] * cos_a[i] * sin_a[j] +
                (a2[i] + a2[j] - 3.0) * sin_a[i] * sin_a[j] -
                3.0 * a[i] * a[j] * cos_a[i] * cos_a[j];

            Wab[i * Na + j] =
                2.0 * aweights[i] * aweights[j] * numerator / (a[i] * a[j]);
        }
    }
}





/**********************************************
 *                                            *
 *             vdw_DF_kernel_bessel           *
 *                                            *
 **********************************************
 *
 *  Purpose:
 *  --------
 *  Compute the Fourier–Bessel transform of a radial kernel φ(r)
 *  used in the van der Waals density functional (vdW-DF) method.
 *  The transform is evaluated numerically on a uniform radial grid
 *  using a Filon-type quadrature with cubic-spline interpolation.
 *
 *       φ_k(k) = (4π / (k q²)) ∫₀^∞ r φ(r) sin((k/q) r) dr
 *
 *  This routine generates the transformed kernel φ_k(k) from φ(r)
 *  sampled at discrete radii r = i·dr (i = 0 … nr).
 *
 *  Method:
 *  -------
 *  1. Construct rφ(r) and its second derivatives via `util_spline`.
 *  2. For each k-point, integrate rφ(r)·sin((k/q)r) between
 *     successive spline intervals using `vdw_DF_Fsin`, a
 *     Filon-type subintegrator for cubic splines.
 *  3. Apply the prefactor 4π / (k q²) to obtain φ_k(k).
 *
 *  Arguments:
 *  ----------
 *    q        : scaling parameter (must be non-zero)
 *    dq       : differential step in q (currently unused)
 *    nr, dr   : number of radial points and radial spacing
 *    phir[]   : input radial kernel φ(r) on [0, nr·dr]
 *    nk, dk   : number of k-points and k-spacing
 *    phik[]   : output array for transformed kernel φ_k(k)
 *    x[]      : work array holding r-values
 *    rphi[]   : work array holding r·φ(r)
 *    rphipp[] : work array holding second derivatives from spline
 *    utmp[]   : temporary workspace for spline construction
 *
 *  Notes:
 *  ------
 *  • The k = 0 term is explicitly set to zero.
 *  • End-point slopes for the spline are clamped by finite
 *    differences; set yp1=ypn=0 for a natural spline if desired.
 *  • A tail-fit correction (φ(r) ~ A/r⁶) can be added for better
 *    convergence at large r, though it is currently omitted.
 *  • Requires external functions:
 *        util_spline(...)    – cubic spline generator
 *        vdw_DF_Fsin(...)    – Filon subintegrator for sin(ωr)
 *
 *  References:
 *  -----------
 *    M. Dion et al., Phys. Rev. Lett. 92, 246401 (2004)
 *    K. Lee et al., Phys. Rev. B 82, 081101(R) (2010)
 *
 **********************************************/
void vdw_DF_kernel_bessel(const double q, const double dq, const int nr, const double dr, 
                          double *phir, const int nk, const double dk, 
                          double *phik, double *x, double *rphi, double *rphipp,double *utmp)
{
   // ---- basic sanity checks ----
   if (nr < 1 || nk < 0) return;                 // nothing to do
   if (dr <= 0.0 || dk <= 0.0) return;           // invalid grid
   if (q == 0.0) throw std::runtime_error("vdw_DF_kernel_bessel: q must be nonzero");

   // 4*pi without relying on nonstandard M_PI
   const double fourpi = 16.0 * std::atan(1.0);

   // optional tail-fit amplitude (currently unused, but keep the estimate)
   const double rmax = dr * static_cast<double>(nr);
   //const double A_tail = phir[nr] * std::pow(rmax, 6); // if φ(r) ~ A/r^6 beyond rmax

   // zero output safely/portably
   std::fill(phik, phik + (nk + 1), 0.0);

   // k = 0 is defined as 0 here
   phik[0] = 0.0;

   // prep arrays: x(i)=r, rphi(i)=r*phi(r)
   for (int i=0; i<=nr; ++i) 
   {
      const double r = dr * static_cast<double>(i);
      x[i]    = r;
      rphi[i] = r * phir[i];
   }

   // clamped slopes from one-sided finite differences (switch to 0,0 for natural spline)
   const double yp1 = (rphi[1]  - rphi[0])     / dr;
   const double ypn = (rphi[nr] - rphi[nr-1])  / dr;

   // spline second derivs into rphipp, workspace utmp
   util_spline(x, rphi, nr+1, yp1, ypn, rphipp, utmp);

   // Filon-like composite integral over [r_{i-1}, r_i] with spline pieces
   for (int k=1; k<=nk; ++k) 
   {
      const double kk        = dk * static_cast<double>(k);
      const double omega     = kk / q;               // frequency in sin(omega r)
      double accum           = 0.0;

      for (int i = 1; i <= nr; ++i) 
      {
         // integrate on [x[i-1], x[i]] the cubic * sin(omega r)
         accum += vdw_DF_Fsin(x[i-1], x[i], rphi[i-1], rphi[i], rphipp[i-1], rphipp[i], omega);
      }

       // scale factor 4π/(k q^2) = (4π/kk) * (1/q^2)
       phik[k] = (fourpi/(kk*q*q))*accum;
   }
}


/*
*     *****************************************
*     *                                       *
*     *         vdw_DF_get_qmesh_filename     *
*     *                                       *
*     *****************************************
*
*     This function returns the filename of the qmesh datafile.
*
*     *** Order of precedence for choosing name                     ***
*     *** 1) value of NWCHEM_QMESH_DATA environment variable        ***
*     *** 2) value of NWCHEM_QMESH_DATA set in $HOME/.nwchemrc file ***
*     *** 3) value of the compiled in library name                  ***
*
*     This is a serial io routine
*
      subroutine vdw_DF_get_qmesh_filename(qmesh_data_name)
      implicit none
      character*(*) qmesh_data_name

#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "bafdecls.fh"

*     **** local variables ****
      logical mprint,hprint,debug,does_it_exist
      logical from_environment,from_compile,from_nwchemrc
      integer iop,lgth,unitf,print_level,i,j
      character*255 qmesh_library

*     **** external functions ****
      logical  util_find_dir
      external util_find_dir

      call util_print_get_level(print_level)
      mprint = print_medium.le.print_level
      hprint = print_high  .le.print_level
      debug  = print_debug .le.print_level
      from_environment = .false.
      from_nwchemrc    = .false.
      from_compile     = .false.

*     **** Try to get from NWCHEM_QMESH_DATA environment variable ****
      call util_getenv('NWCHEM_QMESH_DATA',qmesh_data_name)
      lgth=inp_strlen(qmesh_data_name)
      if (lgth.gt.0) then
         if (util_find_dir(qmesh_data_name)) then
            from_environment = .true.
            goto 99
         else
            write(luout,*)' warning:::::::::::::: from_environment'
            write(luout,*)' NWCHEM_QMESH_DATA set to: <',
     &       qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
            write(luout,*)' but file does not exist !'
            write(luout,*)' using compiled library'
         end if
      end if


*     **** Try to get from NWCHEM_QMESH_DATA nwchemrc ****
*2:   check for NWCHEM_QMESH_DATA defined in users .nwchemrc file
*     assumed structure in .nwchemrc file is variable [whitespace] value
*     one setting per line
*
      qmesh_library='nwchem_qmesh_data'
      call inp_save_state() ! save state of any inp unit
      if(.not.util_nwchemrc_get(qmesh_library,qmesh_data_name)) then
        if (debug) then
          write(luout,*)'util_nwchemrc_get failed'
        endif
      else
        does_it_exist=util_find_dir(qmesh_data_name)
        if (does_it_exist)then
          from_nwchemrc = .true.
          call inp_restore_state() ! restore state of any inp unit
          goto 99
        else
          write(luout,*)' warning:::::::::::::: from_nwchemrc'
          write(luout,*)' NWCHEM_QMESH_DATA set to: <',
     &     qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
          write(luout,*)' but file does not exist !'
          write(luout,*)' using compiled in library'
        endif
      endif
      call inp_restore_state() ! restore state of any inp unit



*     **** Try to get from compile ****
      from_compile = .true.
      call util_nwchem_srcdir(qmesh_data_name)
      qmesh_data_name
     >     =qmesh_data_name(1:inp_strlen(qmesh_data_name))
     >     //"/nwpw/pspw/lib/exchange-correlation/vdw-DF/vdw_qmesh.dat"


 99   continue

      if (from_environment) then
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: environment'
      else if (from_nwchemrc) then
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: .nwchemrc'
      else
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: compiled reference'
      endif
      if (debug) then
         write(luout,*) ' NWCHEM_QMESH_DATA set to: <',
     >    qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
      end if

      return
      end


*/

/**********************************************
 *                                            *
 *        vdw_DF_kernel_gen_data (C++)        *
 *                                            *
 **********************************************
 *
 *  C++17 / MPI port of the Fortran driver.
 *  - Reads qmesh on rank 0, broadcasts to all ranks.
 *  - Builds Wab and grids (r, k).
 *  - Distributes (qa,qb) pairs round-robin across ranks.
 *  - Each worker computes phir -> phik, sends to rank 0.
 *  - Rank 0 spline-postprocesses phik to phik2 and writes both.
 *
 *  Output file binary layout:
 *    int32  Nqs
 *    int32  nk
 *    double kmax
 *    double qmesh[Nqs]
 *    // For each pair (j=1..Nqs, i=1..j) in that order:
 *    double phik[nk+1]
 *    double phik2[nk+1]
 *
 **********************************************/

// ---- Driver ----
void vdw_DF_kernel_gen_data(Parallel *myparall, 
                            const char *out_filename,
                            const char *qmesh_filename)
{
    std::string out_file(out_filename);
    std::string qmesh_file(qmesh_filename);
    //int rank = 0, nprocs = 1;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // ---- read qmesh on rank 0 ----
    int Nqs = 0;
    std::vector<double> qmesh;
    if (myparall->is_master())
    {
        //std::cout << "Generating VDW kernel filename:" << outfile << std::endl;
        //std::cout << "reading qmesh file: " << qmesh_file << std::endl;
        std::ifstream fin(qmesh_file);
        if (!fin) {
            throw std::runtime_error("Error opening qmesh file: " + qmesh_file);
        }
        fin >> Nqs;
        if (!fin || Nqs <= 0) {
            throw std::runtime_error("Error reading Nqs from qmesh file");
        }
        qmesh.resize(Nqs);
        for (int i = 0; i < Nqs; ++i) {
            fin >> qmesh[i];
            if (!fin) throw std::runtime_error("Error reading qmesh value");
        }
    }

   // ---- broadcast Nqs and qmesh ----
   myparall->Brdcst_iValue(0,MASTER,&Nqs);
   if (!(myparall->is_master())) qmesh.resize(Nqs);
   myparall->Brdcst_Values(0,MASTER,Nqs,qmesh.data());


   // ---- constants / grids ----
   const int Na = 256;

   // r-grid
   int    nr   = 2048;                 // 32768 if you want very fine
   double rmax = 100.0;
   double dr   = rmax / static_cast<double>(nr);
   const int nrp1 = nr + 1;

   // k-grid
   int    nk   = 1024;                 // 16384 if you want very fine
   double kmax = 64.0;
   double dk   = kmax / static_cast<double>(nk);
   const int nkp1 = nk + 1;

   // ---- buffers ----
   std::vector<double> a(Na), a2(Na), aweights(Na), cos_a(Na), sin_a(Na);
   std::vector<double> nu(Na), nu1(Na);
   std::vector<double> Wab(static_cast<size_t>(Na) * Na);
   std::vector<double> g(nkp1);

   std::vector<double> phir(nrp1), phik(nkp1), phik0(nkp1), phik2(nkp1);
   std::vector<double> xtmp(nrp1), rphi(nrp1), sphi(nrp1), utmp(nrp1);

   // ---- build Wab and related tables ----
   // (gen_Wab itself calls vdw_GaussLegendre + transforms)
   vdw_DF_kernel_gen_Wab(Na,
                         a.data(), a2.data(),
                         aweights.data(),
                         cos_a.data(), sin_a.data(),
                         Wab.data());

    // ---- |g| grid ----
    for (int i = 0; i <= nk; ++i) g[i] = dk * static_cast<double>(i);


    // ---- open output on rank 0, write header ----
    if (myparall->is_master())
    {
       openfile(6, out_filename, "w");
       iwrite(6, &Nqs, 1);
       iwrite(6, &nk, 1);
       dwrite(6, &kmax, 1);
       dwrite(6, qmesh.data(), Nqs);

       // (Fortran had g write commented; we mirror that and don't write g)
       // fout.write(reinterpret_cast<const char*>(g.data()), sizeof(double)*nkp1);
    }


    // ---- distribute work: pairs (j=1..Nqs, i=1..j) in Fortran order ----
    // Map to 0-based: j=0..Nqs-1, i=0..j
    int nprocs = myparall->np();
    int pcount = 0;
    //const int MASTER = 0;
    int rank = myparall->taskid();

    for (int j = 0; j < Nqs; ++j)
    {
        const double qb = qmesh[j];
        for (int i = 0; i <= j; ++i)
        {
            const double qa = qmesh[i];
            const int owner = static_cast<int>(pcount % nprocs);

            if (rank == owner)
            {
                // ----- worker computes phir -> phik -----
                vdw_DF_kernel_gen_phir(Na,
                                       a.data(), a2.data(),
                                       nu.data(), nu1.data(), Wab.data(),
                                       qa, qb, nr, dr, phir.data());

                const double qave   = 0.5 * (qa + qb);
                const double qdelta = (qb - qa) / (qa + qb);

                vdw_DF_kernel_bessel(qave, qdelta,
                                     nr, dr, phir.data(),
                                     nk, dk, phik.data(),
                                     xtmp.data(), rphi.data(),
                                     sphi.data(), utmp.data());

                if (rank != MASTER) {
                    // send phik to master
                   myparall->dsend(0, 900 + (owner & 0x7FFF), MASTER, nkp1, phik.data()); 
                }
            }

            if (rank == MASTER)
            {
                if (owner != MASTER)
                {
                    // recv into phik0
                    myparall->dreceive(0, 900 + (owner & 0x7FFF), owner, nkp1, phik0.data());

                    // spline to get phik2
                    const double yp1 = (phik0[1]   - phik0[0])     / dk;
                    const double ypn = (phik0[nk]  - phik0[nk-1])  / dk;
                    util_spline(g.data(), phik0.data(), nkp1, yp1, ypn, phik2.data(), utmp.data());

                    // write phik0, phik2
                    dwrite(6, phik0.data(), nkp1);
                    dwrite(6, phik2.data(), nkp1);
                }
                else
                {
                    // locally computed phik
                    const double yp1 = (phik[1]   - phik[0])     / dk;
                    const double ypn = (phik[nk]  - phik[nk-1])  / dk;
                    util_spline(g.data(), phik.data(), nkp1, yp1, ypn, phik2.data(), utmp.data());

                    dwrite(6, phik.data(), nkp1);
                    dwrite(6, phik2.data(), nkp1);
                }
            }

            ++pcount;
        }
    }

    // ---- finalize ----
    if (myparall->is_master()) closefile(6);
    myparall->Barrier();


}


}

