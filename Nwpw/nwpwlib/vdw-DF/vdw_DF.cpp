


#include <cmath>

#include "vdw_DF.hpp"

#include <iostream>
#include "iofmt.hpp"
#include <cstring>



namespace pwdft {


/********************************
 *                              *
 *    vdw_DF_kernel_phi_value   *
 *                              *
 ********************************/
// -----------------------------------------------------------------------------
// vdw_DF_kernel_phi_value
//   C++ translation of the Fortran routine using raw arrays.
// -----------------------------------------------------------------------------
static double vdw_DF_kernel_phi_value(int Na,
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


/********************************
 *                              *
 *   vdw_DF_kernel_gen_phir     *
 *                              *
 ********************************/
// -----------------------------------------------------------------------------
// vdw_DF_kernel_gen_phir
//   C++ translation of the Fortran subroutine
//   Computes phir(i) for radial grid points based on d1,d2 parameters
// -----------------------------------------------------------------------------
static void vdw_DF_kernel_gen_phir(int Na,
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




/********************************
 *                              *
 *       GaussLegendre          *
 *                              *
 ********************************/
static void GaussLegendre(double amin, double amax, int Npoints,
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

/********************************
 *                              *
 *         vdw_DF_Fsin          *
 *                              *
 ********************************/
// -----------------------------------------------------------
//  vdw_DF_Fsin
//  C++ translation of the Fortran function
//  Computes interpolation coefficients involving sin/cos terms
// -----------------------------------------------------------
static double vdw_DF_Fsin(double x0, double x1,
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


/********************************
 *                              *
 *    vdw_DF_kernel_gen_Wab     *
 *                              *
 ********************************/
// -----------------------------------------------------------------------------
// vdw_DF_kernel_gen_Wab
//   C++ translation of the Fortran subroutine.
//
//   Generates Wab(Na×Na) matrix and associated arrays for vdW-DF kernel.
//
//   Parameters:
//     Na        : number of integration points
//     a, a2     : arrays of abscissas and squared abscissas (size Na)
//     aweights  : Gauss–Legendre weights (modified) (size Na)
//     cos_a, sin_a : trigonometric precomputations (size Na)
//     Wab       : 2D matrix (flattened Na×Na, row-major)
// -----------------------------------------------------------------------------

static void vdw_DF_kernel_gen_Wab(int Na,
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
    GaussLegendre(std::atan(amin), std::atan(amax), Na, a, aweights);

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



/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/**
 * @brief Constructor for the `d3db` class.
 *
 * This constructor initializes an instance of the `d3db` class with the given parameters.
 *
 * @param inparall A pointer to a `Parallel` object.
 * @param inmaptype An integer specifying the mapping type.
 * @param nx The number of grid points in the x-direction.
 * @param ny The number of grid points in the y-direction.
 * @param nz The number of grid points in the z-direction.
 */
vdw_DF::vdw_DF()
{


*     **** read and allocate qmesh data ****
      if (taskid.eq.MASTER) THEN
         call vdw_DF_get_qmesh_filename(qmesh_data_name)

         if(.not.util_io_unit(80,90,unitf))
     >     call errquit("vdw-DF cannot get io unit",0,DISK_ERR)

         open(unit=unitf,file=qmesh_data_name,status='old',
     >     form='formatted',ERR=999)

         read(unitf,*,ERR=999,END=999) Nqs
      end if
      call Parallel_Brdcst_ivalue(MASTER,Nqs)



c     **** Langreth kernel data ****
      integer Na
      parameter (Na=256)

   // initialize r and k grid
   // r-grid 
   int    nr  = 2048;
   int    nr1 = nr+1;
   double rmax = 100.0;
   double dr = rmax/dble(nr)

   // kgrid - maximum g=64 and gg=4096 ... 
   int    nk   = 1024
   int    nk1  = nk+1
   double kmax = 64.0
   double dk = kmax/dble(nk)    

   double *qmesh = new (std::nothrow) double[Nqs]();
   double *g     = new (std::nothrow) double[nk1]();
   double *phir  = new  (std::nothrow) double[nr1]();
   double *phir0 = new  (std::nothrow) double[nr1]();
   double *rphi  = new  (std::nothrow) double[nr1]();
   double *sphi  = new  (std::nothrow) double[nr1]();
   double *utmp  = new  (std::nothrow) double[nr1]();
   double *xtmp  = new  (std::nothrow) double[nr1]();

   double *phik  = new  (std::nothrow) double[nk1]();
   double *phik0 = new  (std::nothrow) double[nk1]();
   double *phik2 = new  (std::nothrow) double[nk1]();
   double *a     = new  (std::nothrow) double[Na]();
   double *a2     = new  (std::nothrow) double[Na]();
   double *aweights = new  (std::nothrow) double[Na]();
   double *cos_a = new  (std::nothrow) double[Na]();
   double *sin_a = new  (std::nothrow) double[Na]();
   double *nu = new  (std::nothrow) double[Na]();
   double *nu1 = new  (std::nothrow) double[Na]();
   double *Wab = new  (std::nothrow) double[Na*Na]();



   // delete data
   delete [] Wab;
   delete [] nu1;
   delete [] nu;
   delete [] sin_a;
   delete [] cos_a;
   delete [] aweights;
   delete [] a2;
   delete [] a;
   delete [] phik;
   delete [] phik0;
   delete [] xtmp;
   delete [] utmp;
   delete [] sphi;
   delete [] rphi;
   delete [] phir0;
   delete [] phir;
   delete [] g;
   delete [] qmesh;



      value = BA_pop_stack(Wab(2))
      value = value.and.BA_pop_stack(nu1(2))
      value = value.and.BA_pop_stack(nu(2))
      value = value.and.BA_pop_stack(sin_a(2))
      value = value.and.BA_pop_stack(cos_a(2))
      value = value.and.BA_pop_stack(aweights(2))
      value = value.and.BA_pop_stack(a2(2))
      value = value.and.BA_pop_stack(a(2))
      value = value.and.BA_pop_stack(phik2(2))
      value = value.and.BA_pop_stack(phik0(2))
      value = value.and.BA_pop_stack(phik(2))
      value = value.and.BA_pop_stack(xtmp(2))
      value = value.and.BA_pop_stack(utmp(2))
      value = value.and.BA_pop_stack(sphi(2))
      value = value.and.BA_pop_stack(rphi(2))
      value = value.and.BA_pop_stack(phir0(2))
      value = value.and.BA_pop_stack(phir(2))
      value = value.and.BA_pop_stack(g(2))
      value = value.and.BA_pop_stack(qmesh(2))
      if (.not.value) call errquit('vdw_DF_gen_data:pop stack',0,MA_ERR)
}


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



}



