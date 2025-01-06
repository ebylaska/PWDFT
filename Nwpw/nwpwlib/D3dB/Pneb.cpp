/* Pneb.cpp
   Author - Eric Bylaska

        this class is used for defining 3d parallel maps
*/

/**
 * @brief Class representing the PNEB (Parallel NEB) calculation.
 *
 * The `Pneb` class is responsible for managing the Parallel NEB (PNEB) calculation.
 * It inherits properties and methods from several base classes such as `PGrid` and `d1db`.
 * The PNEB calculation involves complex operations related to parallelization,
 * matrix manipulations, and more.
 */


#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring> //memset()
#include <iostream>
#include <sstream>
#include <stdexcept> // runtime_error()

#include "Pneb.hpp"

#include "blas.h"

namespace pwdft {


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
// Pneb::Pneb(Parallel *inparall, Lattice *inlattice, Control2& control, int
// ispin, int *ne)
//   : PGrid(inparall, inlattice, control),
//     d1db(inparall,control.mapping1d(),ispin,ne)

Pneb::Pneb(Parallel *inparall, Lattice *inlattice, Control2 &control, int ispin, int *ne)
     : d1db(inparall, control.mapping1d(), ispin, ne), PGrid(inparall, inlattice, control) 
{

   int np_i = d1db::parall->np_i();
   int np_j = d1db::parall->np_j();
   int taskid_i = d1db::parall->taskid_i();
   int taskid_j = d1db::parall->taskid_j();
 
   parallelized = (np_j > 1);
 
   s22 = new (std::nothrow) double[8 * ne[0] * ne[0]]();
   s21 = &s22[1 * ne[0] * ne[0]];
   s12 = &s22[2 * ne[0] * ne[0]];
   s11 = &s22[3 * ne[0] * ne[0]];
   sa1 = &s22[4 * ne[0] * ne[0]];
   sa0 = &s22[5 * ne[0] * ne[0]];
   st1 = &s22[6 * ne[0] * ne[0]];
   st2 = &s22[7 * ne[0] * ne[0]];

   if (parallelized) 
   {
      mcq[0] = mcq[1] = 0;
      ncq[0] = ncq[1] = 0;
      for (auto ms=0; ms<ispin; ++ms) 
      {
         ma[ms]  = new (std::nothrow) int[np_i];
         ma1[ms] = new (std::nothrow) int[np_i];
         ma2[ms] = new (std::nothrow) int[np_i];
         mc[ms] = new (std::nothrow) int[np_i];
         na[ms] = new (std::nothrow) int[np_j];
         nc[ms] = new (std::nothrow) int[np_j];
         std::memset(ma[ms], 0,np_i*sizeof(int));
         std::memset(ma1[ms],0,np_i*sizeof(int));
         std::memset(ma2[ms],0,np_i*sizeof(int));
         std::memset(mc[ms], 0,np_i*sizeof(int));

         std::memset(na[ms],0,np_j*sizeof(int));
         std::memset(nc[ms],0,np_j*sizeof(int));
 
         {  auto i=0; auto j=0;
            for (auto k=0; k<ne[ms]; ++k)
            {
               mc[ms][i] += 1;
               nc[ms][j] += 1;
               na[ms][j] += 1;
               i = (i+1)%np_i;
               j = (j+1)%np_j;
            }
         }
 
         ma[ms][taskid_i]  = 2*PGrid::npack(1);
         ma1[ms][taskid_i] = 2*PGrid::nzero(1);
         ma2[ms][taskid_i] = n2ft3d;
         d1db::parall->Vector_ISumAll(1,np_i,ma[ms]);
         d1db::parall->Vector_ISumAll(1,np_i,ma1[ms]);
         d1db::parall->Vector_ISumAll(1,np_i,ma2[ms]);
 
         mcq[ms] = mc[ms][taskid_i];
         ncq[ms] = nc[ms][taskid_j];
 
         mcqmax[ms] = 0;
         for (auto i=0; i<np_i; ++i)
             if (mc[ms][i] > mcqmax[ms])
                mcqmax[ms] = mc[ms][i];
 
         ncqmax[ms] = 0;
         for (auto j=0; j<np_j; ++j)
             if (nc[ms][j] > ncqmax[ms])
                ncqmax[ms] = nc[ms][j];
              
         ncqmax0 = 0;
         for (auto j=0; j<np_j; ++j)
            if (nc[0][j] > ncqmax0)
               ncqmax0 = nc[0][j];
 
         npack1_all = 0;
         nida1_all  = 0;
         n2ft3d_all = 0;
         for (auto i=0; i<np_i; ++i)
         {
            npack1_all += ma[ms][i];
            nida1_all  += ma1[ms][i];
            n2ft3d_all += ma2[ms][i];
         }
      }


      work1 = new (std::nothrow) double[2*2*64*ma2[0][taskid_i]]();
      int nework = 2*64*nc[0][taskid_j];
      if (nework<mcq[0]*ncq[0]) nework = mcq[0]*ncq[0];
      work2 = new (std::nothrow) double[3*nework]();
      bcolwork = new (std::nothrow) double[ne[0]*ncqmax0]();
      bwork2   = new (std::nothrow) double[ne[0]*ncqmax0]();
      //rwork1   = new (std::nothrow) double[2*(na[0][0]*ma2[0][taskid_i])]();
      //rwork2   = new (std::nothrow) double[2*(na[0][0]*ma2[0][taskid_i])]();
      rwork1   = new (std::nothrow) double[2*(na[0][0]*ma2[0][0])]();
      rwork2   = new (std::nothrow) double[2*(na[0][0]*ma2[0][0])]();
      mat_tmp  = new (std::nothrow) double[3*(mcq[0]*ncq[0] + mcq[1]*ncq[1])]();

      mindx[0] = new (std::nothrow) int[(mcq[0]*ncq[0] + mcq[1]*ncq[1])]();
      mindx[1] = new (std::nothrow) int[(mcq[0]*ncq[0])]();
      if (ispin>1) mindx[2] = new (std::nothrow) int[(mcq[1]*ncq[1])]();

      mall[0]  = ne[0]*ne[0] + ne[1]*ne[1];
      mall[1]  = ne[0]*ne[0];
      mall[2]  = ne[1]*ne[1];
      mpack[0] = mcq[0]*ncq[0] + mcq[1]*ncq[1];
      mpack[1] = mcq[0]*ncq[0];
      mpack[2] = mcq[1]*ncq[1];

      int indx0=0;
      int indx1=0;
      int indx2=0;
      int jj=0;
      int jcur=0;
      for (auto j=0; j<ne[0]; ++j)
      {
         auto ii=0;
         auto icur=0;
         for (auto i=0; i<ne[0]; ++i)
         {
            if ((icur==taskid_i) && (jcur==taskid_j))
            {
               mindx[0][indx0] = i + j*ne[0];
               mindx[1][indx1] = i + j*ne[0];
               ++indx0;
               ++indx1;
            }
            ++ii;
            if (ii>=mc[0][icur])
            {
               ++icur;
               ii=0;
            }
         }
         ++jj;
         if (jj>=nc[0][jcur])
         {
            ++jcur;
            jj=0;
         }
      }
      if (ispin>1)
      {
         int jj=0;
         int jcur=0;
         for (auto j=0; j<ne[1]; ++j)
         {
            auto ii=0;
            auto icur=0;
            for (auto i=0; i<ne[1]; ++i)
            {
               if ((icur==taskid_i) && (jcur==taskid_j))
               {
                  mindx[0][indx0] = i + j*ne[1] + ne[0]*ne[0];
                  mindx[2][indx2] = i + j*ne[1];
                  ++indx0;
                  ++indx2;
               }
               ++ii;
               if (ii>=mc[1][icur])
               {
                  ++icur;
                  ii=0;
               }
            }
            ++jj;
            if (jj>=nc[1][jcur])
            {
               ++jcur;
               jj=0;
            }
         }
      }
   }
 
   g_rnd_algorithm = control.initial_psi_random_algorithm();

   io_norbs_max = control.io_norbs_max();
   io_buffer    = control.io_buffer();
}

/*************************************
 *                                   *
 *      Pneb::g_generate1_random     *
 *                                   *
 *************************************/
void Pneb::g_generate1_random(double *psi) 
{
   int ms, n, indx, i, pj, qj, taskid_j;
   int filling[4], nfft[3];
   double zvalue[2];
   // double *tmp2 = new (std::nothrow) double[n2ft3d]();
 
   double tmp2[n2ft3d];
 
   nfft[0] = nx;
   nfft[1] = ny;
   nfft[2] = nz;
 
   taskid_j = d1db::parall->taskid_j();
   for (ms = 0; ms < ispin; ++ms)
      for (n = 0; n < ne[ms]; ++n) 
      {
         util_getfilling(n, nfft, filling, zvalue);
        
         qj = msntoindex(ms, n);
         pj = msntop(ms, n);
        
         if (pj == taskid_j) 
         {
            r_zero(tmp2);
            d3db::c_setpw(filling, zvalue, tmp2);
            d3db::c_addrandom(tmp2);
           
            PGrid::c_pack(1, tmp2);
            indx = 2 * PGrid::npack(1) * qj;
            PGrid::cc_pack_copy(1, tmp2, psi + indx);
            PGrid::c_pack_noimagzero(1, psi + indx);
         }
      }
 
   // delete [] tmp2;
}

void Pneb::g_generate2_random(double *psi) 
{
   double tmp2[n2ft3d];
 
   int taskid_j = d1db::parall->taskid_j();
   for (auto ms = 0; ms < ispin; ++ms)
      for (auto n = 0; n < ne[ms]; ++n) 
      {
         int qj = msntoindex(ms, n);
         int pj = msntop(ms, n);
         if (pj == taskid_j) {
           d3db::r_setrandom(tmp2);
           d3db::r_zero_ends(tmp2);
           d3db::rc_fft3d(tmp2);
        
           PGrid::c_pack(1, tmp2);
           int indx = 2 * PGrid::npack(1) * qj;
           PGrid::cc_pack_copy(1, tmp2, psi + indx);
         }
      }
}


/*************************************
 *                                   *
 *  Pneb::g_generate_excited_random  *
 *                                   *
 *************************************/
void Pneb::g_generate_excited_random(const int nex[], double *psi_excited) 
{
   int taskid = d1db::parall->taskid();
   util_random(taskid + 91);
   double tmp2[n2ft3d];

   int taskid_j = d1db::parall->taskid_j();
   for (auto ms = 0; ms < ispin; ++ms)
   for (auto n = 0; n < nex[ms]; ++n) 
   {
      d3db::r_setrandom(tmp2);
      d3db::r_zero_ends(tmp2);
      d3db::rc_fft3d(tmp2);
   
      PGrid::c_pack(1, tmp2);
      int indx = 2 * PGrid::npack(1) * n;
      PGrid::cc_pack_copy(1, tmp2, psi_excited + indx);
   }
}

/*************************************
 *                                   *
 *      Pneb::g_generate_random      *
 *                                   *
 *************************************/
void Pneb::g_generate_random(double *psi) {
   int taskid = d1db::parall->taskid();
   util_random(taskid + 91);
 
   if (g_rnd_algorithm == 1)
     this->g_generate1_random(psi);
   else
     this->g_generate2_random(psi);
}

/*************************************
 *                                   *
 *           Pneb::g_read            *
 *                                   *
 *************************************/
 /**
 * @brief Read data from an input unit within a parallel computing framework.
 *
 * This function is responsible for reading data from a specified input unit based
 * on the input parameters. It operates within a parallel computing framework and uses
 * the PGrid class for data manipulation.
 *
 * @param iunit An integer specifying the input unit to read data from.
 * @param psi A pointer to a double array where the read data will be stored.
 *
 * @return None.
 */
void Pneb::g_read(const int iunit, double *psi) 
{
   int ms, n, indx, i, pj, qj, taskid_j;
   double *tmp2 = new (std::nothrow) double[n2ft3d]();
 
   taskid_j = d1db::parall->taskid_j();
 
   for (ms=0; ms<ispin; ++ms)
      for (n=0; n<ne[ms]; ++n) 
      {
         qj = msntoindex(ms, n);
         pj = msntop(ms, n);
         c_read(iunit, tmp2, pj);
         if (pj == taskid_j) 
         {
            indx = 2*PGrid::npack(1)*qj;
            PGrid::c_pack(1, tmp2);
            PGrid::cc_pack_copy(1, tmp2, psi + indx);
         }
      }
 
   delete[] tmp2;
}


/*************************************
 *                                   *
 *           Pneb::g_read_ne         *
 *                                   *
 *************************************/
void Pneb::g_read_ne(const int iunit, const int *ne0, double *psi) 
{
   int ms, n, indx, i, pj, qj, taskid_j;
   double *tmp2 = new (std::nothrow) double[n2ft3d]();
   //std::unique_ptr<double*> tmp2(new double[n2ft3d]());

 
   taskid_j = d1db::parall->taskid_j();
 
   for (ms=0; ms<ispin; ++ms)
      for (n=0; n<ne0[ms]; ++n) 
      {
         qj = msntoindex(ms, n);
         pj = msntop(ms, n);
         if (n<ne0[ms])
         {
            c_read(iunit, tmp2, pj);
         }
         else
         {
            d3db::r_setrandom(tmp2);
            d3db::r_zero_ends(tmp2);
            d3db::rc_fft3d(tmp2);
         }
         
         if (pj == taskid_j) 
         {
            indx = 2*PGrid::npack(1)*qj;
            PGrid::c_pack(1, tmp2);
            PGrid::cc_pack_copy(1, tmp2, psi + indx);
         }
      }
 
   delete[] tmp2;
}

/*************************************
 *                                   *
 *      Pneb::g_read_ne_reverse      *
 *                                   *
 *************************************/
void Pneb::g_read_ne_reverse(const int iunit, const int *ne0, double *psi)
{        
   int ms, n, indx, i, pj, qj, taskid_j;
   double *tmp2 = new (std::nothrow) double[n2ft3d]();
   //std::unique_ptr<double[]> tmp2(new double[n2ft3d]());
      
      
   taskid_j = d1db::parall->taskid_j();
            
   for (ms=0; ms<ispin; ++ms)
      for (n=ne0[ms]-1; n>=0; --n) 
      {     
         qj = msntoindex(ms, n);
         pj = msntop(ms, n);
         if (n<ne0[ms])
         {  
            c_read(iunit, tmp2, pj);
         }
         else
         {  
            d3db::r_setrandom(tmp2);
            d3db::r_zero_ends(tmp2);
            d3db::rc_fft3d(tmp2);
         }
            
         if (pj == taskid_j) 
         {
            indx = 2*PGrid::npack(1)*qj;
            PGrid::c_pack(1, tmp2);
            PGrid::cc_pack_copy(1, tmp2, psi + indx);
         }
      }

   delete[] tmp2;
}


/*************************************
 *                                   *
 *           Pneb::g_write           *
 *                                   *
 *************************************/
/**
 * @brief Write data to an output unit within a parallel computing framework.
 *
 * This function is responsible for writing data to a specified output unit based on
 * the input parameters. It operates within a parallel computing framework and uses
 * the PGrid class for data manipulation.
 *
 * @param iunit An integer specifying the output unit to write data to.
 * @param psi A pointer to a double array containing the data to be written.
 *
 * @return None.
 */
void Pneb::g_write(const int iunit, double *psi) 
{
   d3db::parall->Barrier();
   int ms, n, indx, i, pj, qj, taskid_j;
   double *tmp2 = new (std::nothrow) double[n2ft3d]();

   taskid_j = d1db::parall->taskid_j();
   for (ms = 0; ms < ispin; ++ms)
   for (n = 0; n < ne[ms]; ++n) 
   {
      qj = msntoindex(ms, n);
      pj = msntop(ms, n);
      if (pj == taskid_j) 
      {
         indx = 2 * PGrid::npack(1) * qj;
         PGrid::cc_pack_copy(1, psi + indx, tmp2);
         PGrid::c_unpack(1, tmp2);
      }
      if (io_buffer)
         c_write_buffer(iunit,tmp2,pj);
      else
         c_write(iunit,tmp2,pj);
   }

   delete[] tmp2;
}

/*************************************
 *                                   *
 *      Pneb::g_write_excited        *
 *                                   *
 *************************************/
void Pneb::g_write_excited(const int iunit, const int nex[], double *psi)
{
   d3db::parall->Barrier();
   int taskid_j = d1db::parall->taskid_j();

   double *tmp2 = new (std::nothrow) double[n2ft3d]();

   for (auto ms = 0; ms < ispin; ++ms)
   for (auto n = 0; n < nex[ms]; ++n)
   {
      int indx = 2 * PGrid::npack(1) * n;
      PGrid::cc_pack_copy(1, psi + indx, tmp2);
      PGrid::c_unpack(1, tmp2);

      if (io_buffer)
         c_write_buffer(iunit,tmp2,taskid_j);
      else
         c_write(iunit,tmp2,taskid_j);
   }

   delete[] tmp2;
}


/*************************************
 *                                   *
 *     Pneb::gg_traceall_excited     *
 *                                   *
 *************************************/
double Pneb::gg_traceall_excited(const int nex[], double *psi1, double *psi2) 
{
   int n, indx;
   double sum = 0.0;
 
   indx = 0;
   for (n = 0; n < (nex[0] + nex[1]); ++n) 
   {
      sum += PGrid::cc_pack_idot(1, psi1 + indx, psi2 + indx);
      indx += 2 * PGrid::npack(1);
   }
   //if (ispin == 1)
   //   sum *= 2.0;
 
   return d3db::parall->SumAll(0, sum);
}

/*************************************
 *                                   *
 *           Pneb::gg_traceall       *
 *                                   *
 *************************************/
/**
 * @brief Calculate the trace of the inner product between two sets of spinor wave functions.
 *
 * This function calculates the trace of the inner product between two sets of spinor wave functions, represented by the arrays `psi1` and `psi2`. The trace is computed by summing the inner product of corresponding elements in the arrays.
 *
 * @param psi1 A pointer to the first set of spinor wave functions.
 * @param psi2 A pointer to the second set of spinor wave functions.
 *
 * @return The trace of the inner product between the two sets of spinor wave functions.
 *
 * @note This function assumes that the sizes of `psi1` and `psi2` are compatible, and it operates on elements in pairs. The `neq` and `ispin` variables are used to determine the sizes and characteristics of the wave functions.
 */
double Pneb::gg_traceall(double *psi1, double *psi2) 
{
   int n, indx;
   double sum = 0.0;
 
   indx = 0;
   for (n = 0; n < (neq[0] + neq[1]); ++n) 
   {
      sum += PGrid::cc_pack_idot(1, psi1 + indx, psi2 + indx);
      indx += 2 * PGrid::npack(1);
   }
   if (ispin == 1)
      sum *= 2.0;
 
   return d3db::parall->SumAll(0, sum);
}

/*************************************
 *                                   *
 *       Pneb::gg_traceall_occ       *
 *                                   *
 *************************************/
double Pneb::gg_traceall_occ(double *psi1, double *psi2, double *occ) 
{
   double sum = 0.0;
 
   int indx = 0;
   for (auto n=0; n<(neq[0]+neq[1]); ++n) 
   {
      sum += occ[n]*PGrid::cc_pack_idot(1, psi1 + indx, psi2 + indx);
      indx += 2*PGrid::npack(1);
   }
   if (ispin == 1)
      sum *= 2.0;
 
   return d3db::parall->SumAll(0, sum);
}

/*************************************
 *                                   *
 *           Pneb::gg_copy           *
 *                                   *
 *************************************/
/**
 * @brief Copy the elements of one array to another.
 *
 * This function copies the elements from the source array `psi1` to the destination
 * array `psi2`.
 *
 * @param psi1  The source array.
 * @param psi2  The destination array where the elements are copied.
 */
void Pneb::gg_copy(double *psi1, double *psi2) 
{
   int one = 1;
   int nsize = 2 * (neq[0] + neq[1]) * PGrid::npack(1);
   std::memcpy(psi2, psi1, nsize * sizeof(double));
   // DCOPY_PWDFT(nsize, psi1, one, psi2, one);
}


/*************************************
 *                                   *
 *           Pneb::gg_SMul           *
 *                                   *
 *************************************/
/**
 * @brief Multiply an array by a scalar: psi2 = alpha * psi1.
 *
 * This function multiplies each element of the `psi1` array by the scalar `alpha`
 * and stores the result in the `psi2` array.
 *
 * @param alpha The scalar value to multiply with each element of `psi1`.
 * @param psi1  The source array.
 * @param psi2  The destination array where the result is stored.
 */
void Pneb::gg_SMul(double alpha, double *psi1, double *psi2) 
{
   int nsize = 2 * (neq[0] + neq[1]) * PGrid::npack(1);
   for (int i = 0; i < nsize; ++i)
      psi2[i] = alpha * psi1[i];
}
// void Pneb::g_SMul1(double alpha,double *psi1)
//{
//    int nsize = 2*(neq[0]+neq[1])*PGrid::npack(1);
//    for (int i=0; i<nsize; ++i)
//       psi1[i] *= alpha;
// }

/*************************************
 *                                   *
 *           Pneb::gg_daxpy          *
 *                                   *
 *************************************/
/**
 * @brief Perform a linear combination of two arrays: psi2 += alpha * psi1.
 *
 * This function computes the linear combination `psi2 += alpha * psi1` for arrays `psi1` and `psi2`.
 *
 * @param alpha The scaling factor to apply to the `psi1` array before adding to `psi2`.
 * @param psi1  The source array.
 * @param psi2  The destination array where the result is stored.
 */
void Pneb::gg_daxpy(double alpha, double *psi1, double *psi2) 
{
   int nsize = 2 * (neq[0] + neq[1]) * PGrid::npack(1);
   for (int i = 0; i < nsize; ++i)
      psi2[i] += alpha * psi1[i];
}


/*************************************
 *                                   *
 *           Pneb::g_Scale           *
 *                                   *
 *************************************/
/**
 * @brief Scale the elements of an array by a given factor.
 *
 * This function scales each element of the input array `psi1` by the specified factor `alpha`.
 *
 * @param alpha The scaling factor to apply to the array elements.
 * @param psi1  The array to be scaled.
 */
void Pneb::g_Scale(double alpha, double *psi1) 
{
   int nsize = 2 * (neq[0] + neq[1]) * PGrid::npack(1);
   for (int i = 0; i < nsize; ++i)
      psi1[i] *= alpha;
}


/*************************************
 *                                   *
 *           Pneb::gg_Sum2           *
 *                                   *
 *************************************/
/**
 * @brief Element-wise addition of two arrays.
 *
 * This function performs element-wise addition, adding the elements of
 * array `psi1` to the elements of array `psi2`, and stores the result in `psi2`.
 *
 * @param psi1  The array to add to `psi2`.
 * @param psi2  The array to which `psi1` is added, and where the result is stored.
 */
void Pneb::gg_Sum2(double *psi1, double *psi2) 
{
   int nsize = 2 * (neq[0] + neq[1]) * PGrid::npack(1);
   for (int i = 0; i < nsize; ++i)
      psi2[i] += psi1[i];
}


/*************************************
 *                                   *
 *         Pneb::gg_Minus2           *
 *                                   *
 *************************************/
/**
 * @brief Element-wise subtraction of one array from another.
 *
 * This function performs element-wise subtraction, subtracting the elements of
 * array `psi1` from the elements of array `psi2`, and stores the result in `psi2`.
 *
 * @param psi1  The array to subtract from `psi2`.
 * @param psi2  The array from which `psi1` is subtracted, and where the result is stored.
 */
void Pneb::gg_Minus2(double *psi1, double *psi2) 
{
   int nsize = 2 * (neq[0] + neq[1]) * PGrid::npack(1);
   for (int i = 0; i < nsize; ++i)
      psi2[i] -= psi1[i];
}


/*************************************
 *                                   *
 *         Pneb::ggg_Minus           *
 *                                   *
 *************************************/
/**
 * @brief Element-wise subtraction of two sets of spinor wave functions.
 *
 * This function performs element-wise subtraction of two sets of spinor wave functions, `psi1` and `psi2`, and stores the result in `psi3`.
 *
 * @param psi1 A pointer to the first set of spinor wave functions.
 * @param psi2 A pointer to the second set of spinor wave functions.
 * @param psi3 A pointer to the resulting set of spinor wave functions after subtraction.
 *
 * @note This function assumes that the sizes of `psi1`, `psi2`, and `psi3` are compatible, and it operates on elements in pairs. The `neq` variable is used to determine the sizes of the wave functions.
 */
void Pneb::ggg_Minus(double *psi1, double *psi2, double *psi3) 
{
   int nsize = 2 * (neq[0] + neq[1]) * PGrid::npack(1);
   for (int i = 0; i < nsize; ++i)
      psi3[i] = psi1[i] - psi2[i];
}


/*************************************
 *                                   *
 *         Pneb::g_zero              *
 *                                   *
 *************************************/
void Pneb::g_zero(double *psi2) 
{
   int one = 1;
   int zero = 0;
   int nsize = 2 * (neq[0] + neq[1]) * PGrid::npack(1);
   double rzero = 0.0;

   // dcopy_(&nsize,&rzero,&zero,psi2,&one);
   std::memset(psi2, 0, nsize * sizeof(double));
}


/*************************************
 *                                   *
 *         Pneb::gh_fftb             *
 *                                   *
 *************************************/
/**
 * @brief Perform a forward Fast Fourier Transform (FFT) of a set of wave functions using the PFFT3B library.
 *
 * This function performs a forward Fast Fourier Transform (FFT) of a set of wave functions, `psi`, using the PFFT3B library and stores the result in `psi_r`. The transform is performed along one dimension.

 * @param psi A pointer to the input set of wave functions to be transformed.
 * @param psi_r A pointer to the resulting transformed wave functions.
 *
 * @note This function uses the cr_pfft3b_queue library to perform the FFT. It queues the input wave functions for transformation and retrieves the transformed results.

 * @warning The behavior of this function depends on the initialization and configuration of the pfft3b_queue library. Ensure that the library is properly initialized and configured before calling this function.
 */
void Pneb::gh_fftb(double *psi, double *psi_r) 
{
   nwpw_timing_function ftimer(1);
   int n, done;
   int indx1, indx1n, shift1;
   int indx2, indx2n, shift2;
   int nffts_pipeline = PGrid::nffts_max;
   //std::cout << "gh_fftb: nffts=" << nffts_pipeline << std::endl;
 
   n = neq[0] + neq[1];
   shift1 = 2 * PGrid::npack(1);
   shift2 = n2ft3d;
   indx1 = indx1n = 0;
   indx2 = indx2n = 0;
   done = 0;
   while (!done) 
   {
      if (indx1 < n) 
      {
         // Assuming n, indx1, and nffts_pipeline are defined
         int idum = std::min(n - indx1, nffts_pipeline);

         cr_pfft3b_queuein(1, idum, psi + indx1n);

         indx1n += idum*shift1;
         indx1  += idum;
      }
      if (cr_pfft3b_queuefilled() || (indx1 >= n)) 
      {
         // Assuming n, indx2, and nffts_pipeline are defined
         int jdum = std::min(n - indx2, nffts_pipeline);

         cr_pfft3b_queueout(1, jdum, psi_r + indx2n);

         indx2n += jdum*shift2;
         indx2  += jdum;;
         //++indx2;
      }
      done = ((indx1 >= n) && (indx2 >= n));
   }
}

/*************************************
 *                                   *
 *         Pneb::hr_aSumSqr          *
 *                                   *
 *************************************/
/**
 * @brief Calculate the sum of squares of components of a set of wave functions with an optional scaling factor.
 *
 * This function calculates the sum of squares of the components of a set of wave functions, `psir`, with an optional scaling factor `alpha`. The result is stored in the `dn` array.

 * @param alpha An optional scaling factor applied to the sum of squares.
 * @param psir A pointer to the input set of wave functions.
 * @param dn A pointer to the resulting sum of squares of components.
 *
 * @note This function performs the calculation for each component of the wave functions in the `psir` array and accumulates the results in the `dn` array.
 */
void Pneb::hr_aSumSqr(const double alpha, double *psir, double *dn) 
{
   int n,ms,k,indx0,indx1;
   int one = 1;
   int zero = 0;
   int nsize = n2ft3d * ispin;
   double rzero = 0.0;
 
   // dcopy_(&nsize,&rzero,&zero,dn,&one);
   std::memset(dn,0,nsize*sizeof(double));
 
   indx0 = 0;
   indx1 = 0;
   for (ms=0; ms<ispin; ++ms) 
   {
      for (n=0; n<(neq[ms]); ++n) 
      {
         for (k=0; k<n2ft3d; ++k)
            dn[indx0+k] += alpha*psir[indx1+k]*psir[indx1+k];
         indx1 += n2ft3d;
      }
      indx0 += n2ft3d;
   }
   d3db::parall->Vector_SumAll(2, ispin*n2ft3d, dn);
}

/*************************************
 *                                   *
 *         Pneb::hr_aSumSqr_occ      *
 *                                   *
 *************************************/
void Pneb::hr_aSumSqr_occ(const double alpha, const double *occ, double *psir, double *dn) 
{
   int one = 1;
   int zero = 0;
   int nsize = n2ft3d * ispin;
   //double rzero = 0.0;

   // dcopy_(&nsize,&rzero,&zero,dn,&one);
   std::memset(dn,0,nsize*sizeof(double));

   int indx0 = 0;
   int indx1 = 0;
   for (auto ms=0; ms<ispin; ++ms)
   {
      for (auto q=0; q<(neq[ms]); ++q)
      {
         double wf = occ[msntoindex(ms,q)];
         for (auto k=0; k<n2ft3d; ++k)
            dn[indx0+k] += alpha*psir[indx1+k]*psir[indx1+k]*wf;
         indx1 += n2ft3d;
      }
      indx0 += n2ft3d;
   }
   d3db::parall->Vector_SumAll(2, ispin*n2ft3d, dn);
}


/*************************************
 *                                   *
 *         Pneb::hhr_aSumMul         *
 *                                   *
 *************************************/
/*
   This routine calculates the perturbation density.

   dn12(r) = Sum_i (psi0_i(r)*psi1_i(r) + psi1_i(r)*psi0_i(r) )

   where

   psi_i(r) = psi0_i(r) + lmbda*psi1_i(r)
   dn12(r) = n(r) = Sum_i (psi_i(r)*psi_i(r))
                  = Sum_i ( psi0_i(r)*psi0_i(r)
                          + lmbda*(psi0_i(r)*psi1_i(r) + psi1_i(r)*psi0_i(r))
                          + lmbda^2*psi1_i(r)*psi1_i(r) )
*/
/**
 * @brief Calculate the perturbation density with a mixed product of two sets of wave functions.
 *
 * This function calculates the perturbation density `dn12(r)` using a mixed product of two sets of wave functions `psi0_i(r)` and `psi1_i(r)` with an optional scaling factor `alpha`. The result is stored in the `dn12` array.

 * @param alpha An optional scaling factor applied to the mixed product.
 * @param psir0 A pointer to the first set of wave functions.
 * @param psir1 A pointer to the second set of wave functions.
 * @param dn12 A pointer to the resulting perturbation density.
 *
 * @note This function performs the calculation for each component of the wave functions in the `psir0` and `psir1` arrays and accumulates the results in the `dn12` array according to the formula:
 * \f[
 * dn12(r) = \sum_i \left( \psi0_i(r) \psi1_i(r) + \psi1_i(r) \psi0_i(r) \right)
 * \f]
 * where \f$\psi_i(r) = \psi0_i(r) + \lambda \psi1_i(r)\f$, and \f$\lambda\f$ is the optional scaling factor.
 */
void Pneb::hhr_aSumMul(const double alpha, const double *psir0,
                       const double *psir1, double *dn12) 
{
   int n, ms, k, indx0, indx1;
   int one = 1;
   int zero = 0;
   int nsize = n2ft3d * ispin;
   double rzero = 0.0;
 
   std::memset(dn12, 0, nsize * sizeof(double));
 
   indx0 = 0;
   indx1 = 0;
   for (ms = 0; ms < ispin; ++ms) 
   {
      for (n = 0; n < (neq[ms]); ++n) 
      {
         for (k=0; k<n2ft3d; ++k)
            dn12[indx0+k] += alpha*( psir0[indx1+k]*psir1[indx1+k] 
                                   + psir1[indx1+k]*psir0[indx1+k]);
         indx1 += n2ft3d;
      }
      indx0 += n2ft3d;
   }
   d3db::parall->Vector_SumAll(2, ispin * n2ft3d, dn12);
}

/*************************************
 *                                   *
 *      Pneb::ggm_sym_Multiply       *
 *                                   *
 *************************************/
/**
 * @brief Multiply two sets of wave functions with symmetry considerations.
 *
 * This function multiplies two sets of wave functions `psi1` and `psi2` and stores the result in the `hml` matrix while considering symmetry properties. The function is parallelized for performance improvements.

 * @param psi1 A pointer to the first set of wave functions.
 * @param psi2 A pointer to the second set of wave functions.
 * @param hml A pointer to the resulting matrix.
 *
 * @note This function performs matrix multiplication of `psi1` and `psi2` with symmetry considerations. The result is stored in the `hml` matrix, which should be pre-allocated with sufficient memory. Symmetry is taken into account to optimize the calculation.
 */
void Pneb::ggm_sym_Multiply(double *psi1, double *psi2, double *hml) 
{
   nwpw_timing_function ftimer(15);
 
   int npack1 = 2*PGrid::npack(1);
   int ng0    = 2*PGrid::nzero(1);
 
   int one = 1;
   double rzero = 0.0;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;
 
   if (parallelized) 
   {
      auto taskid_i = d1db::parall->taskid_i();
      auto taskid_j = d1db::parall->taskid_j();
      auto ishift2 = mcq[0]*ncq[0];
      for (auto ms=0; ms<ispin; ++ms) 
      {
          if (ne[ms]>0)
          {
             auto shift0 = ms*neq[0]*npack1;
             auto shift2 = ms*ishift2;
             d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                           ne[ms],ne[ms],npack1_all,128,
                           psi1+shift0,psi2+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                           mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                           work1,work2);
          }
      }
      std::memset(hml,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,hml);
      d1db::parall->Vector_SumAll(0,mall[0],hml);

      // Symmetrize hml
      for (auto ms=0; ms<ispin; ++ms)
      {
         int n = ne[ms];
         int mshift0 = ms*ne[0]*ne[0];
         for (auto k=0; k<n; ++k)
            for (auto j=k+1; j<n; ++j)
               hml[mshift0+j+k*n] = hml[mshift0+k+j*n];
      }
   } 
   else 
   {
      auto shift0  = 0;
      auto mshift0 = 0;
      for (auto ms=0; ms<ispin; ++ms) 
      {
         auto n = ne[ms];
         d3db::mygdevice.TN1_dgemm(npack1,n,rtwo,psi1+shift0,psi2+shift0,rzero,hml+mshift0);
        
         if (ng0 > 0) 
         {
            auto shift1  = shift0;
            auto mshift1 = mshift0;
            for (auto k=1; k<=n; ++k) 
            {
               DGEMM_PWDFT((char *)"T", (char *)"N",k,one,ng0,
                           rmone,
                           psi1+shift0,npack1, 
                           psi2+shift1,npack1, 
                           rone, 
                           hml+mshift1, k);
               shift1 += npack1;
               mshift1 += n;
            }
         }
         for (auto k=0; k<n; ++k)
            for (auto j=k+1; j<n; ++j)
               hml[mshift0+j+k*n] = hml[mshift0+k+j*n];
        
         shift0  += npack1*ne[0];
         mshift0 += ne[0]*ne[0];
      }
      d3db::parall->Vector_SumAll(1,ne[0]*ne[0]+ne[1]*ne[1],hml);
   }
   //if (d3db::parall->is_master())
   //   std::cout << "hml= [" << Efmt(15,10) 
   //                         << hml[0] << " " << hml[4] << " " << hml[8]  << " " << hml[12]  << std::endl 
   //             << "      " << hml[1] << " " << hml[5] << " " << hml[9]  << " " << hml[13]  << std::endl
   //             << "      " << hml[2] << " " << hml[6] << " " << hml[10] << " " << hml[14]  << std::endl
   //             << "      " << hml[3] << " " << hml[7] << " " << hml[11] << " " << hml[15]  << "]"<< std::endl << std::endl;
}


/*************************************
 *                                   *
 *        Pneb::ggm_Multiply         *
 *                                   *
 *************************************/
/**
 * @brief Multiply two sets of wave functions without considering symmetry.
 *
 * This function multiplies two sets of wave functions `psi1` and `psi2` and stores the result in the `hml` matrix without considering symmetry properties. The function is parallelized for performance improvements.

 * @param psi1 A pointer to the first set of wave functions.
 * @param psi2 A pointer to the second set of wave functions.
 * @param hml A pointer to the resulting matrix.
 *
 * @note This function performs matrix multiplication of `psi1` and `psi2` without taking symmetry into account. The result is stored in the `hml` matrix, which should be pre-allocated with sufficient memory. Parallelization is applied for improved performance.
 */
void Pneb::ggm_Multiply(double *psi1, double *psi2, double *hml) 
{
   nwpw_timing_function ftimer(15);
   int n, shift0, mshift0;
 
   int one = 1;
   int npack1 = 2 * PGrid::npack(1);
   int ng0 = 2 * PGrid::nzero(1);
 
   double rzero = 0.0;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;
 
   if (parallelized) 
   {
      //std::ostringstream msg;
      //msg << "NWPW Error: ggm_Multiply() parallelized is NOT supported\n"
      //    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      //throw(std::runtime_error(msg.str()));
      auto taskid_i = d1db::parall->taskid_i();
      auto taskid_j = d1db::parall->taskid_j();
      auto ishift2 = mcq[0]*ncq[0];
      for (auto ms=0; ms<ispin; ++ms)
      {  
          if (ne[ms]>0) 
          { 
             auto shift0 = ms*neq[0]*npack1;
             auto shift2 = ms*ishift2;
             d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                           ne[ms],ne[ms],npack1_all,128,
                           psi1+shift0,psi2+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                           mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                           work1,work2);
          }
      }
      std::memset(hml,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,hml);
      d1db::parall->Vector_SumAll(0,mall[0],hml);

   } 
   else 
   {
      shift0 = 0;
      mshift0 = 0;
      for (auto ms = 0; ms < ispin; ++ms) 
      {
         n = ne[ms];
         d3db::mygdevice.TN1_dgemm(npack1,n,rtwo,psi1+shift0,psi2+shift0,rzero,hml+mshift0);
        
         if (ng0 > 0)
           DGEMM_PWDFT((char *)"T", (char *)"N", n, n, ng0, rmone, psi1 + shift0,
                       npack1, psi2 + shift0, npack1, rone, hml + mshift0, n);
        
         shift0 += npack1 * ne[0];
         mshift0 += ne[0] * ne[0];
      }
      d3db::parall->Vector_SumAll(1, ne[0] * ne[0] + ne[1] * ne[1], hml);
   }
}

/*************************************
 *                                   *
 *      Pneb::ffm_sym_Multiply       *
 *                                   *
 *************************************/
void Pneb::ffm_sym_Multiply(const int mb, double *psi1, double *psi2, double *hml) 
{
   nwpw_timing_function ftimer(15);
   int ms, ms1, ms2, ishift2, j, k, n, shift0, shift1, mshift0, mshift1, nn;
 
   int one = 1;
   int npack1 = 2 * PGrid::npack(1);
   int ng0 = 2 * PGrid::nzero(1);
 
   double rzero = 0.0;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;
 
   if (parallelized) 
   {
      //std::ostringstream msg;
      //msg << "NWPW Error: ffm_sym_Multiply() parallelized is NOT supported\n"
      //    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      //throw(std::runtime_error(msg.str()));
      if (mb==-1)
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = mcq[0]*ncq[0];
         nn = mcqmax[0]*ncqmax[0]+mcqmax[1]*ncqmax[1];
      }
      else
      {
         ms1 = mb;
         ms2 = mb+1;
         ishift2 = 0;
         nn = mcqmax[mb]*ncqmax[mb];
      }

      auto taskid_i = d1db::parall->taskid_i();
      auto taskid_j = d1db::parall->taskid_j();

      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi1+shift0,psi1+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(hml,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,hml);
      d1db::parall->Vector_SumAll(0,mall[0],hml);

      // Symmetrize hml
      for (auto ms=0; ms<ispin; ++ms)
      {
         mshift0 = ms*ne[0]*ne[0];
         n = ne[ms];
         for (auto k=0; k<n; ++k)
            for (auto j=k+1; j<n; ++j)
               hml[mshift0+j+k*n] = hml[mshift0+k+j*n];
      }

   } 
   else 
   {
      if (mb == -1) 
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = ne[0] * ne[0];
         nn = ne[0] * ne[0] + ne[1] * ne[1];
         shift0 = 0;
         mshift0 = 0;
      } 
      else 
      {
         ms1 = mb;
         ms2 = mb + 1;
         ishift2 = 0;
         nn = ne[mb] * ne[mb];
         shift0 = mb * ne[0] * npack1;
         mshift0 = 0;
      }
      for (ms = ms1; ms < ms2; ++ms) 
      {
         n = ne[ms];
         d3db::mygdevice.TN1_dgemm(npack1,n,rtwo,&psi1[shift0],&psi2[shift0],rzero,&hml[mshift0]);
        
         if (ng0 > 0) 
         {
            shift1 = shift0;
            mshift1 = mshift0;
            for (k = 1; k <= n; ++k) 
            {
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           &psi1[shift0], npack1, &psi2[shift1], npack1, rone, &hml[mshift1], k);
              
               shift1 += npack1;
               mshift1 += n;
            }
         }
         for (k = 0; k < n; ++k)
         for (j = k + 1; j < n; ++j)
            hml[mshift0 + j + k * n] = hml[mshift0 + k + j * n];
        
         shift0  += npack1*ne[0];
         mshift0 += ne[0]*ne[0];
      }
      d3db::parall->Vector_SumAll(1, nn, hml);
   }
}


/*************************************
 *                                   *
 *        Pneb::ffm_Multiply         *
 *                                   *
 *************************************/
void Pneb::ffm_Multiply(const int mb, double *psi1, double *psi2, double *hml) 
{
   nwpw_timing_function ftimer(15);
   int ms, ms1, ms2, ishift2, j, k, n, shift0, shift1, mshift0, mshift1, nn;
 
   int one = 1;
   //int ng = 2 * PGrid::npack(1);
   int npack1 = 2 * PGrid::npack(1);
   int ng0 = 2 * PGrid::nzero(1);
 
   double rzero = 0.0;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;
 
   if (parallelized) 
   {
      //std::ostringstream msg;
      //msg << "NWPW Error: ffm_Multiply() parallelized is NOT supported\n"
      //    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      //throw(std::runtime_error(msg.str()));

      if (mb==-1)
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = mcq[0]*ncq[0];
         nn = mcqmax[0]*ncqmax[0]+mcqmax[1]*ncqmax[1];
      }
      else
      {
         ms1 = mb;
         ms2 = mb+1;
         ishift2 = 0;
         nn = mcqmax[mb]*ncqmax[mb];
      }

      auto taskid_i = d1db::parall->taskid_i();
      auto taskid_j = d1db::parall->taskid_j();

      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi1+shift0,psi1+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(hml,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,hml);
      d1db::parall->Vector_SumAll(0,mall[0],hml);
         
   } 
   else 
   {
      if (mb == -1) 
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = ne[0] * ne[0];
         nn = ne[0] * ne[0] + ne[1] * ne[1];
         shift0 = 0;
         mshift0 = 0;
      } 
      else 
      {
         ms1 = mb;
         ms2 = mb + 1;
         ishift2 = 0;
         nn = ne[mb] * ne[mb];
         shift0 = mb * ne[0] * npack1;
         mshift0 = 0;
      }
      for (ms = ms1; ms < ms2; ++ms) 
      {
         n = ne[ms];
         d3db::mygdevice.TN1_dgemm(npack1,n,rtwo,&psi1[shift0],&psi2[shift0],rzero,hml+mshift0);
        
         if (ng0 > 0)
            DGEMM_PWDFT((char *)"T", (char *)"N", n, n, ng0, rmone, 
                        psi1+shift0, npack1, 
                        psi2+shift0, npack1, 
                        rone, hml+mshift0, n);
         shift0 += npack1 * ne[0];
         mshift0 += ne[0] * ne[0];
      }
      d3db::parall->Vector_SumAll(1, nn, hml);
   }
}

/*************************************
 *                                   *
 *      Pneb::ffm3_sym_Multiply      *
 *                                   *
 *************************************/
 /* 
  Description:                                    
  This function performs symmetric matrix         
  multiplication and manipulation operations. It  
  operates on input arrays psi1 and psi2, and     
  calculates and stores the results in arrays s11,
  s21, and s22. The calculations are performed based 
  on the given parameters and matrix properties.    
                                                   
  Parameters:                                     
  - mb: Index specifying the block of matrices to 
        operate on. If mb is -1, full matrix      
        multiplication is performed. Otherwise, a 
        subset of matrices is used for calculations. 
  - psi1, psi2: Pointers to double arrays containing 
                input matrix data.                  
  - s11, s21, s22: Pointers to double arrays where 
                  the results will be stored.     
*/

void Pneb::ffm3_sym_Multiply(const int mb, double *psi1, double *psi2,
                             double *s11, double *s21, double *s22) 
{
   nwpw_timing_function ftimer(15);
   int ms1, ms2, ishift2, n, shift0, shift1, mshift0, mshift1, nn;
   int one = 1;
   int npack1 = 2*PGrid::npack(1);
   int ng0    = 2*PGrid::nzero(1);
 
   double rzero = 0.0;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;
 
   if (parallelized) 
   {
      //std::ostringstream msg;
      //msg << "NWPW Error: ffm3_sym_Multiply() parallelized is NOT supported\n"
      //    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      //throw(std::runtime_error(msg.str()));
      if (mb==-1)
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = mcq[0]*ncq[0];
         nn = mcqmax[0]*ncqmax[0]+mcqmax[1]*ncqmax[1];
      }
      else
      {
         ms1 = mb;
         ms2 = mb+1;
         ishift2 = 0;
         nn = mcqmax[mb]*ncqmax[mb];
      }

      auto taskid_i = d1db::parall->taskid_i();
      auto taskid_j = d1db::parall->taskid_j();
      //auto ishift2 = mcq[0]*ncq[0];

      //s11
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi1+shift0,psi1+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s11,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s11);
      d1db::parall->Vector_SumAll(0,mall[0],s11);

      //s21
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi2+shift0,psi1+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s21,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s21);
      d1db::parall->Vector_SumAll(0,mall[0],s21);

      //s22
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi2+shift0,psi2+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s22,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s22);
      d1db::parall->Vector_SumAll(0,mall[0],s22);

      // Symmetrize hml
      for (auto ms=0; ms<ispin; ++ms)
      {
         mshift0 = ms*ne[0]*ne[0];
         n = ne[ms];
         for (auto k=0; k<n; ++k)
            for (auto j=k+1; j<n; ++j)
            {
               s11[mshift0+j+k*n] = s11[mshift0+k+j*n];
               s21[mshift0+j+k*n] = s21[mshift0+k+j*n];
               s22[mshift0+j+k*n] = s22[mshift0+k+j*n];
            }
      }
   } 
   else 
   {
      if (mb == -1) 
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = ne[0]*ne[0];
         nn = ne[0]*ne[0] + ne[1]*ne[1];
         shift0 = 0;
         mshift0 = 0;
      } 
      else 
      {
         ms1 = mb;
         ms2 = mb + 1;
         ishift2 = 0;
         nn = ne[mb]*ne[mb];
         //shift0 = mb*ne[0]*ng;
         shift0 = mb*ne[0]*npack1;
         mshift0 = 0;
      }
      for (auto ms=ms1; ms<ms2; ++ms) 
      {
         n = ne[ms];
        
         d3db::mygdevice.TN3_dgemm(npack1,n,rtwo,psi1+shift0,psi2+shift0,rzero,s11+mshift0,s21+mshift0,s22+mshift0);
        
         if (ng0 > 0) 
         {
            shift1  = shift0;
            mshift1 = mshift0;
            for (auto k=1; k<=n; ++k) 
            {
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi1+shift0,npack1, 
                           psi1+shift1,npack1, 
                           rone, 
                           s11+mshift1,k);
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi1+shift0,npack1, 
                           psi2+shift1,npack1, 
                           rone, 
                           s21+mshift1,k);
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi2+shift0,npack1, 
                           psi2+shift1,npack1, 
                           rone, 
                           s22+mshift1,k);
               shift1  += npack1;
               mshift1 += n;
            }
         }
        
         for (auto k=0; k<n; ++k) 
         {
            for (auto j=k+1; j<n; ++j) 
            {
               s11[mshift0+j+k*n] = s11[mshift0+k+j*n];
               s21[mshift0+j+k*n] = s21[mshift0+k+j*n];
               s22[mshift0+j+k*n] = s22[mshift0+k+j*n];
            }
         }
        
         shift0  += npack1*ne[0];
         mshift0 += ne[0]*ne[0];
      }
      d3db::parall->Vector_SumAll(1, nn, s11);
      d3db::parall->Vector_SumAll(1, nn, s21);
      d3db::parall->Vector_SumAll(1, nn, s22);
   }
}

/*************************************
 *                                   *
 *  Pneb::ffm3_Fulls21_sym_Multiply  *
 *                                   *
 *************************************/
void Pneb::ffm3_Fulls21_sym_Multiply(const int mb, double *psi1, double *psi2,
                             double *s11, double *s21, double *s22)
{
   nwpw_timing_function ftimer(15);
   int ms1, ms2, ishift2, n, shift0, shift1, mshift0, mshift1, nn;
   int one = 1;
   int npack1 = 2*PGrid::npack(1);
   int ng0    = 2*PGrid::nzero(1);

   double rzero = 0.0;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;

   if (parallelized) 
   {
      //std::ostringstream msg;
      //msg << "NWPW Error: ffm3_sym_Multiply() parallelized is NOT supported\n"
      //    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      //throw(std::runtime_error(msg.str()));
      if (mb==-1)
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = mcq[0]*ncq[0];
         nn = mcqmax[0]*ncqmax[0]+mcqmax[1]*ncqmax[1];
      }
      else
      {
         ms1 = mb;
         ms2 = mb+1;
         ishift2 = 0;
         nn = mcqmax[mb]*ncqmax[mb];
      }

      auto taskid_i = d1db::parall->taskid_i();
      auto taskid_j = d1db::parall->taskid_j();
      //auto ishift2 = mcq[0]*ncq[0];

      //s11
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi1+shift0,psi1+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s11,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s11);
      d1db::parall->Vector_SumAll(0,mall[0],s11);


      //s21
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi2+shift0,psi1+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s21,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s21);
      d1db::parall->Vector_SumAll(0,mall[0],s21);

      //s22
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi2+shift0,psi2+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s22,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s22);
      d1db::parall->Vector_SumAll(0,mall[0],s22);


      // Symmetrize hml
      for (auto ms=0; ms<ispin; ++ms)
      {
         mshift0 = ms*ne[0]*ne[0];
         n = ne[ms];
         for (auto k=0; k<n; ++k)
            for (auto j=k+1; j<n; ++j)
            {
               s11[mshift0+j+k*n] = s11[mshift0+k+j*n];
               //s21[mshift0+j+k*n] = s21[mshift0+k+j*n];
               s22[mshift0+j+k*n] = s22[mshift0+k+j*n];
            }
      }
   }
   else
   {
      if (mb == -1)
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = ne[0]*ne[0];
         nn = ne[0]*ne[0] + ne[1]*ne[1];
         shift0 = 0;
         mshift0 = 0;
      }
      else
      {
         ms1 = mb;
         ms2 = mb + 1;
         ishift2 = 0;
         nn = ne[mb]*ne[mb];
         //shift0 = mb*ne[0]*ng;
         shift0 = mb*ne[0]*npack1;
         mshift0 = 0;
      }
      for (auto ms=ms1; ms<ms2; ++ms)
      {
         n = ne[ms];

         d3db::mygdevice.TN3_FullCab_dgemm(npack1,n,rtwo,psi1+shift0,psi2+shift0,rzero,s11+mshift0,s21+mshift0,s22+mshift0);

         if (ng0 > 0)
         {
            shift1  = shift0;
            mshift1 = mshift0;
            for (auto k=1; k<=n; ++k)
            {
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi1+shift0,npack1,
                           psi1+shift1,npack1,
                           rone,
                           s11+mshift1,k);
               //DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
               //            psi1+shift0,npack1, 
               //            psi2+shift1,npack1,
               //            rone,
               //            s21+mshift1,k);
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi2+shift0,npack1,
                           psi2+shift1,npack1,
                           rone,
                           s22+mshift1,k);
               shift1  += npack1;
               mshift1 += n;
            }
            DGEMM_PWDFT((char *)"T", (char *)"N", n, n, ng0, rmone,
                        psi1+shift0,npack1, 
                        psi2+shift0,npack1,
                        rone,
                        s21+mshift0,n);
         }

         for (auto k=0; k<n; ++k)
         {
            for (auto j=k+1; j<n; ++j)
            {
               s11[mshift0+j+k*n] = s11[mshift0+k+j*n];
               //s21[mshift0+j+k*n] = s21[mshift0+k+j*n];
               s22[mshift0+j+k*n] = s22[mshift0+k+j*n];
            }
         }

         shift0  += npack1*ne[0];
         mshift0 += ne[0]*ne[0];
      }
      d3db::parall->Vector_SumAll(1, nn, s11);
      d3db::parall->Vector_SumAll(1, nn, s21);
      d3db::parall->Vector_SumAll(1, nn, s22);
   }
}






/*************************************
 *                                   *
 *      Pneb::ffm4_sym_Multiply      *
 *                                   *
 *************************************/
/**
 * @brief Multiply two sets of wave functions with symmetry consideration.
 *
 * This function multiplies two sets of wave functions `psi1` and `psi2` with symmetry consideration and stores the result in the `s11`, `s21`, `s12`, and `s22` matrices. The function is parallelized for performance improvements.
 *
 * @param mb An integer specifying the band index. If `mb` is -1, the function operates on all bands.
 * @param psi1 A pointer to the first set of wave functions.
 * @param psi2 A pointer to the second set of wave functions.
 * @param s11 A pointer to the resulting matrix for symmetry operation 11.
 * @param s21 A pointer to the resulting matrix for symmetry operation 21.
 * @param s12 A pointer to the resulting matrix for symmetry operation 12.
 * @param s22 A pointer to the resulting matrix for symmetry operation 22.
 *
 * @note This function performs matrix multiplication of `psi1` and `psi2` with symmetry consideration and stores the results in four different matrices `s11`, `s21`, `s12`, and `s22`. The matrices should be pre-allocated with sufficient memory. Parallelization is applied for improved performance.
 */
void Pneb::ffm4_sym_Multiply(const int mb, double *psi1, double *psi2,
                             double *s11, double *s21, double *s12, double *s22) 
{
   nwpw_timing_function ftimer(15);
   int ms, ms1, ms2, ishift2, j, k, n, shift0, shift1, mshift0, mshift1, nn;
   int one = 1;
   int npack1 = 2 * PGrid::npack(1);
   int ng0 = 2 * PGrid::nzero(1);
 
   double rzero = 0.0;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;
 
   if (parallelized) 
   {
      //std::ostringstream msg;
      //msg << "NWPW Error: ffm4_sym_Multiply() parallelized is NOT supported\n"
      //    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      //throw(std::runtime_error(msg.str()));
      if (mb==-1)
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = mcq[0]*ncq[0];
         nn = mcqmax[0]*ncqmax[0]+mcqmax[1]*ncqmax[1];
      }
      else
      {
         ms1 = mb;
         ms2 = mb+1;
         ishift2 = 0;
         nn = mcqmax[mb]*ncqmax[mb];
      }

      auto taskid_i = d1db::parall->taskid_i();
      auto taskid_j = d1db::parall->taskid_j();
      //auto ishift2 = mcq[0]*ncq[0];

      //s11
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi1+shift0,psi1+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s11,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s11);
      d1db::parall->Vector_SumAll(0,mall[0],s11);

      //s21 
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi2+shift0,psi1+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s21,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s21);
      d1db::parall->Vector_SumAll(0,mall[0],s21);

      //s12 
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi1+shift0,psi2+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s12,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s12);
      d1db::parall->Vector_SumAll(0,mall[0],s12);

      //s22
      for (auto ms=0; ms<ispin; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift2 = ms*ishift2;
            d1db::DMatrix_dgemm2c(d1db::parall, &mygdevice,
                          ne[ms],ne[ms],npack1_all,128,
                          psi2+shift0,psi2+shift0, ma[ms][taskid_i],ma[ms],ma1[ms],na[ms],
                          mat_tmp+shift2,mc[ms][taskid_i],mc[ms],nc[ms],
                          work1,work2);
         }
      }
      std::memset(s22,0,mall[0]*sizeof(double));
      t_bindexcopy(mpack[0],mindx[0],mat_tmp,s22);
      d1db::parall->Vector_SumAll(0,mall[0],s22);

      // Symmetrize hml
      for (auto ms=0; ms<ispin; ++ms)
      {
         mshift0 = ms*ne[0]*ne[0];
         n = ne[ms];
         for (auto k=0; k<n; ++k)
            for (auto j=k+1; j<n; ++j)
            {
               s11[mshift0+j+k*n] = s11[mshift0+k+j*n];
               s21[mshift0+j+k*n] = s21[mshift0+k+j*n];
               s12[mshift0+j+k*n] = s12[mshift0+k+j*n];
               s22[mshift0+j+k*n] = s22[mshift0+k+j*n];
            }
      }
   } 
   else 
   {
      if (mb == -1) 
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = ne[0] * ne[0];
         nn = ne[0] * ne[0] + ne[1] * ne[1];
         shift0 = 0;
         mshift0 = 0;
      } 
      else 
      {
         ms1 = mb;
         ms2 = mb + 1;
         ishift2 = 0;
         nn = ne[mb] * ne[mb];
         shift0 = mb * ne[0] * npack1;
         mshift0 = 0;
      }
      for (ms = ms1; ms < ms2; ++ms) 
      {
         n = ne[ms];
        
         d3db::mygdevice.TN4_dgemm(npack1,n,rtwo,psi1+shift0,psi2+shift0,rzero,
                                   s11+mshift0,s21+mshift0,s12+mshift0,s22+mshift0);
         if (ng0 > 0) 
         {
            shift1 = shift0;
            mshift1 = mshift0;
            for (k = 1; k <= n; ++k) 
            {
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi1+shift0, npack1, psi1+shift1, npack1, rone, s11+mshift1, k);
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi2+shift0, npack1, psi1+shift1, npack1, rone, s21+mshift1, k);
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi1+shift0, npack1, psi2+shift1, npack1, rone, s12+mshift1, k);
               DGEMM_PWDFT((char *)"T", (char *)"N", k, one, ng0, rmone,
                           psi2+shift0, npack1, psi2+shift1, npack1, rone, s22+mshift1, k);
               shift1 += npack1;
               mshift1 += n;
            }
         }
        
         for (k = 0; k < n; ++k) 
         {
            for (j = k + 1; j < n; ++j) 
            {
               s11[mshift0 + j+k*n] = s11[mshift0 + k+j*n];
               s21[mshift0 + j+k*n] = s21[mshift0 + k+j*n];
               s12[mshift0 + j+k*n] = s12[mshift0 + k+j*n];
               s22[mshift0 + j+k*n] = s22[mshift0 + k+j*n];
            }
         }
        
         shift0 += npack1* ne[0];
         mshift0 += ne[0] * ne[0];
      }
      d3db::parall->Vector_SumAll(1, nn, s11);
      d3db::parall->Vector_SumAll(1, nn, s21);
      d3db::parall->Vector_SumAll(1, nn, s12);
      d3db::parall->Vector_SumAll(1, nn, s22);
   }
}

/*************************************
 *                                   *
 *        Pneb::fmf_Multiply         *
 *                                   *
 *************************************/
/**
 * @brief Multiply wave functions with a Hermitian matrix and scaling.
 *
 * This function multiplies a set of wave functions `psi1` with a Hermitian matrix `hml` and scales the result by `alpha`, then adds it to another set of wave functions `psi2` scaled by `beta`. The function is parallelized for performance improvements.
 *
 * @param mb An integer specifying the band index. If `mb` is -1, the function operates on all bands.
 * @param psi1 A pointer to the first set of wave functions.
 * @param hml A pointer to the Hermitian matrix.
 * @param alpha A scaling factor for the multiplication of `psi1` with `hml`.
 * @param psi2 A pointer to the second set of wave functions.
 * @param beta A scaling factor for the addition of the result to `psi2`.
 *
 * @note This function performs the multiplication of `psi1` with the Hermitian matrix `hml`, scales the result by `alpha`, and then adds it to `psi2` scaled by `beta`. The matrices and vectors should be pre-allocated with sufficient memory. Parallelization is applied for improved performance.
 */
void Pneb::fmf_Multiply(const int mb, double *psi1, double *hml, double alpha,
                        double *psi2, double beta) 
{
   nwpw_timing_function ftimer(16);
   int ms, ms1, ms2, n, nn,shift1, mshift1, ishift2,ishift3;
   int npack1 = 2 * PGrid::npack(1);

/*
   if (d1db::parall->is_master())
      std::cout << "into HML= [" << Efmt(15,10)
                            << hml[0] << " " << hml[4] << " " << hml[8]  << " " << hml[12]  << std::endl
                << "      " << hml[1] << " " << hml[5] << " " << hml[9]  << " " << hml[13]  << std::endl
                << "      " << hml[2] << " " << hml[6] << " " << hml[10] << " " << hml[14]  << std::endl
                << "      " << hml[3] << " " << hml[7] << " " << hml[11] << " " << hml[15]  << "]"<< std::endl << std::endl;
                */

   if (parallelized) 
   {
      if (mb==-1)
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = mcq[0]*ncq[0];
         ishift3 = ne[0]*ne[0];
         nn = mcqmax[0]*ncqmax[0]+mcqmax[1]*ncqmax[1];
      }
      else
      {
         ms1 = mb;
         ms2 = mb+1;
         ishift2 = 0;
         ishift3 = 0;
         nn = mcqmax[mb]*ncqmax[mb];
      }
      auto taskid_i = d1db::parall->taskid_i();
      auto taskid_j = d1db::parall->taskid_j();




      for (auto ms=ms1; ms<ms2; ++ms)
      {
         if (ne[ms]>0)
         {
            auto shift0 = ms*neq[0]*npack1;
            auto shift3 = ms*ishift3;
            d1db::DMatrix_dgemm1_rot2(d1db::parall, &mygdevice,
                          npack1_all,ne[ms],ne[ms],128,
                          alpha,
                          psi1+shift0,ma[ms][taskid_i],ma[ms],na[ms],
                          hml+shift3,mcq[ms],mc[ms],nc[ms],
                          beta,
                          psi2+shift0,ma[ms][taskid_i],ma[ms],na[ms],
                          bcolwork,bwork2,rwork1,rwork2);
         }
      }
   } 
   else 
   {
      if (mb == -1) 
      {
         ms1 = 0;
         ms2 = ispin;
         ishift2 = ne[0]*ne[0];
         shift1 = 0;
         mshift1 = 0;
      } 
      else 
      {
         ms1 = mb;
         ms2 = mb + 1;
         ishift2 = 0;
         shift1 = mb * ne[0]*npack1;
         mshift1 = 0;
      }
      for (ms=ms1; ms<ms2; ++ms) 
      {
         n = ne[ms];
        
         d3db::mygdevice.NN_dgemm(npack1, n, alpha, psi1+shift1, hml+mshift1, beta, psi2+shift1);
         shift1 += ne[0]*npack1;
         mshift1 += ishift2;
      }
   }
}


/*************************************
 *                                   *
 *            Pneb::ggm_SVD          *
 *                                   *
 *************************************/
/**
 * @brief Perform Singular Value Decomposition (SVD) of a matrix.
 *
 * This function performs the Singular Value Decomposition (SVD) of a given matrix 'A'
 * and returns the components 'U', 'S', and 'V' such that A = U * S * V^T.
 *
 * @param[in] A Input matrix to be decomposed.
 * @param[out] U Output matrix 'U' with left singular vectors.
 * @param[out] S Output array 'S' containing singular values.
 * @param[out] V Output matrix 'V' with right singular vectors.
 */
void Pneb::ggm_SVD(double *A, double *U, double *S, double *V) 
{
   int n, indx;
   double *tmp2 = new (std::nothrow) double[neq[0] + neq[1]]();
 
   /* generate V and Sigma^2 */
   //if (d3db::parall->is_master()) std::cout << "generate V and Sigma2:" << std::endl;
   ggm_sym_Multiply(A, A, V);
   m_diagonalize(V, S);
   //if (d3db::parall->is_master()) 
   //    std::cout << "Sigma=" << Efmt(15,10)  << S[0] << " " << S[1] << " " << S[2] << " " << S[3] << std::endl << std::endl;
 
   /* generate U*Sigma */
   fmf_Multiply(-1, A, V, 1.0, U, 0.0);
 
   //if (d3db::parall->is_master()) 
   //    std::cout << "AFTER FMF_Multiply" << std::endl;

   /* normalize U*sigma */
   indx = 0;
   for (n = 0; n < (neq[0] + neq[1]); ++n) 
   {
      tmp2[n] = PGrid::cc_pack_idot(1, U+indx, U+indx);
      indx += 2*PGrid::npack(1);
   }
   d3db::parall->Vector_SumAll(1, neq[0] + neq[1], tmp2);
 
   for (n = 0; n < (neq[0] + neq[1]); ++n)
      tmp2[n] = 1.0 / std::sqrt(tmp2[n]);
 
   indx = 0;
   for (n = 0; n < (neq[0] + neq[1]); ++n) 
   {
      PGrid::c_pack_SMul(1, tmp2[n], U+indx);
      indx += 2 * PGrid::npack(1);
   }
 
   /* calculated sqrt(S^2) */
   for (n = 0; n < (neq[0] + neq[1]); ++n) 
   {
      if (S[n] < 0.0)
         S[n] = std::fabs(S[n]);
      S[n] = std::sqrt(S[n]);
   }
 
   delete[] tmp2;
}


/*************************************
 *                                   *
 *            Pneb::m_scal           *
 *                                   *
 *************************************/
/**
 * @brief Scale the elements of a matrix by a scalar value.
 *
 * This function scales all elements of the input matrix 'hml' by the specified scalar value 'alpha'.
 *
 * @param[in] alpha The scalar value to scale the matrix by.
 * @param[in,out] hml The input matrix to be scaled, and the resulting scaled matrix.
 */
void Pneb::m_scal(double alpha, double *hml) {
  int one = 1;
  int nsize = ne[0] * ne[0] + ne[1] * ne[1];

  DSCAL_PWDFT(nsize, alpha, hml, one);
}

/*************************************
 *                                   *
 *        Pneb::m_diag_scal          *
 *                                   *
 *************************************/
#define eta 1.0e-9
void Pneb::m_diag_scal(const double *occ, double *hml) 
{
   int mshift  = 0;
   int mshift1 = 0;
   for (auto ms=0; ms<ispin; ++ms)
   {
      for (auto i=0; i<ne[ms]; ++i)
         hml[i + i*ne[ms] + mshift] *= (occ[i+mshift1]+eta);

       mshift1 += ne[0];
       mshift  += ne[0]*ne[0];
   }
}

/*************************************
 *                                   *
 *       Pneb::m_diag_scal_inv       *
 *                                   *
 *************************************/
void Pneb::m_diag_scal_inv(const double *occ, double *hml) 
{
   int mshift  = 0;
   int mshift1 = 0;
   for (auto ms=0; ms<ispin; ++ms)
   {
      for (auto i=0; i<ne[ms]; ++i)
         hml[i + i*ne[ms] + mshift] /= (occ[i+mshift1]+eta);

       mshift1 += ne[0];
       mshift  += ne[0]*ne[0];
   }
}


/*************************************
 *                                   *
 *            Pneb::m_trace          *
 *                                   *
 *************************************/
/**
 * @brief Calculate the trace of a matrix.
 *
 * This function calculates the trace of the input matrix 'hml', which is the sum of the diagonal elements.
 *
 * @param[in] hml The input matrix for which the trace is to be calculated.
 * @return The trace of the matrix.
 */
double Pneb::m_trace(double *hml) 
{
   int ms, i;
   int mshift = 0;
   double sum = 0.0;
   for (ms = 0; ms < ispin; ++ms) 
   {
      for (i = 0; i < ne[ms]; ++i)
         sum += hml[i + i * ne[ms] + mshift];
      mshift += ne[0] * ne[0];
   }
   return sum;
}

/*************************************
 *                                   *
 *            Pneb::m_trace_occ      *
 *                                   *
 *************************************/
double Pneb::m_trace_occ(double *hml, double *occ) 
{
   int ms, i;
   int mshift1 = 0;
   int mshift2 = 0;
   double sum = 0.0;
   for (ms = 0; ms < ispin; ++ms) 
   {
      for (i=0; i<ne[ms]; ++i)
         sum += hml[i + i*ne[ms] + mshift2]*occ[i+mshift1];
      mshift1 += ne[0];
      mshift2 += ne[0]*ne[0];
   }
   return sum;
}


/*************************************
 *                                   *
 *        Pneb::m_diagonalize        *
 *                                   *
 *************************************/
/**
 * @brief Diagonalize a matrix and compute its eigenvalues.
 *
 * This function diagonalizes the input matrix 'hml' and computes its eigenvalues, which are stored in the 'eig' array.
 *
 * @param[in,out] hml The input matrix to be diagonalized.
 * @param[out] eig An array to store the computed eigenvalues.
 */
void Pneb::m_diagonalize(double *hml, double *eig) 
{
   nwpw_timing_function ftimer(17);
 
   if (mparallelized)
   {
      std::ostringstream msg;
      msg << "NWPW Error: m_diagonalize() mparallelized is NOT supported\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
   } 
   else 
   {
      int n = ne[0] + ne[1];
      int nn = ne[0]*ne[0] + ne[1]*ne[1];
 
      if (d1db::parall->is_master())
         d3db::mygdevice.NN_eigensolver(ispin, ne, hml, eig);
      d1db::parall->Brdcst_Values(0, 0, nn, hml);
      d1db::parall->Brdcst_Values(0, 0, n, eig);
   }
}


void Pneb::m_diagonalize(const int n, double *hml, double *eig) 
{
   nwpw_timing_function ftimer(17);

   if (mparallelized)
   {
      std::ostringstream msg;
      msg << "NWPW Error: m_diagonalize() mparallelized is NOT supported\n"
          << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      throw(std::runtime_error(msg.str()));
   }
   else
   {
      int nn = n*n;

      if (d1db::parall->is_master())
         d3db::mygdevice.NN_eigensolver0(n, hml, eig);
      d1db::parall->Brdcst_Values(0, 0, nn, hml);
      d1db::parall->Brdcst_Values(0, 0, n, eig);
   }
}


/*************************************
 *                                   *
 *     Pneb::m_scale_s22_s21_s11     *
 *                                   *
 *************************************/
/**
 * @brief Scale matrices S22, S21, and S11 for NEB calculations.
 *
 * This function scales the matrices S22, S21, and S11, which are used in the context of NEB (Nudged Elastic Band) calculations.
 * The scaling depends on the parameters 'mb' and 'dte'. If 'mb' is -1, it scales all matrices for all spin channels; otherwise, it scales matrices for the specified spin channel.
 *
 * @param[in] mb The spin channel index to which the scaling is applied. If -1, scaling is applied to all spin channels.
 * @param[in] dte The scaling factor for the matrices.
 * @param[in,out] s22 The S22 matrix to be scaled.
 * @param[in,out] s21 The S21 matrix to be scaled.
 * @param[in,out] s11 The S11 matrix to be scaled.
 */
void Pneb::m_scale_s22_s21_s11(const int mb, const double dte, double *s22, double *s21, double *s11) 
{
   int j, k, ms, ms1, ms2, ishift2, indx0, indx, indxt;
 
   if (mb == -1) {
     ms1 = 0;
     ms2 = ispin;
     ishift2 = ne[0] * ne[0];
   } else {
     ms1 = mb;
     ms2 = mb + 1;
     ishift2 = 0;
   }
 
   for (ms = ms1; ms < ms2; ++ms) {
     indx0 = ms * ishift2;
     for (k = 0; k < ne[ms]; ++k) {
       s22[indx0] = (1.0 - s22[indx0]) * (0.5 / dte);
       s21[indx0] = (1.0 - s21[indx0]) * (0.5);
       s11[indx0] *= -0.5 * dte;
 
       indx = indx0 + 1;
       indxt = indx0 + ne[ms];
       for (j = (k + 1); j < ne[ms]; ++j) {
         s22[indx] *= (-0.5 / dte);
         s22[indxt] *= (-0.5 / dte);
         s21[indx] *= -0.5;
         s21[indxt] *= -0.5;
         s11[indx] *= -0.5 * dte;
         s11[indxt] *= -0.5 * dte;
 
         indx += 1;
         indxt += ne[ms];
       }
       indx0 += (ne[ms] + 1);
     }
   }
}


/*************************************
 *                                   *
 *   Pneb::m_scale_s22_s21_s12_s11   *
 *                                   *
 *************************************/
/**
 * @brief Scale matrices S22, S21, S12, and S11 for NEB calculations.
 *
 * This function scales the matrices S22, S21, S12, and S11, which are used in the context of NEB (Nudged Elastic Band) calculations.
 * The scaling depends on the parameters 'mb' and 'dte'. If 'mb' is -1, it scales all matrices for all spin channels; otherwise, it scales matrices for the specified spin channel.
 *
 * @param[in] mb The spin channel index to which the scaling is applied. If -1, scaling is applied to all spin channels.
 * @param[in] dte The scaling factor for the matrices.
 * @param[in,out] s22 The S22 matrix to be scaled.
 * @param[in,out] s21 The S21 matrix to be scaled.
 * @param[in,out] s12 The S12 matrix to be scaled.
 * @param[in,out] s11 The S11 matrix to be scaled.
 */
void Pneb::m_scale_s22_s21_s12_s11(const int mb, const double dte, double *s22,
                                   double *s21, double *s12, double *s11) 
{
   int j, k, ms, ms1, ms2, ishift2, indx0, indx, indxt;
 
   if (mb == -1) 
   {
      ms1 = 0;
      ms2 = ispin;
      ishift2 = ne[0] * ne[0];
   } 
   else 
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift2 = 0;
   }
 
   for (ms = ms1; ms < ms2; ++ms) 
   {
      indx0 = ms * ishift2;
      for (k = 0; k < ne[ms]; ++k) 
      {
         s22[indx0] = (1.0 - s22[indx0]) * (0.5 / dte);
         s21[indx0] = (1.0 - s21[indx0]) * (0.5);
         s12[indx0] = (1.0 - s12[indx0]) * (0.5);
         s11[indx0] *= -0.5 * dte;
        
         indx = indx0 + 1;
         indxt = indx0 + ne[ms];
         for (j = (k + 1); j < ne[ms]; ++j) 
         {
            s22[indx] *= (-0.5 / dte);
            s22[indxt] *= (-0.5 / dte);
            s21[indx] *= -0.5;
            s21[indxt] *= -0.5;
            s12[indx] *= -0.5;
            s12[indxt] *= -0.5;
            s11[indx] *= -0.5 * dte;
            s11[indxt] *= -0.5 * dte;
           
            indx += 1;
            indxt += ne[ms];
         }
         indx0 += (ne[ms] + 1);
      }
   }
}


/*************************************
 *                                   *
 *      Pneb::mm_SCtimesVtrans       *
 *                                   *
 *************************************/
/**
 * @brief Multiply a matrix S by the transpose of a matrix Vt, with scaling.
 *
 * This function multiplies a matrix S by the transpose of a matrix Vt while applying scaling factors based on the parameters 'mb' and 't'. If 'mb' is -1, the operation is performed for all spin channels; otherwise, it is performed for the specified spin channel.
 *
 * @param[in] mb The spin channel index to which the operation is applied. If -1, the operation is applied to all spin channels.
 * @param[in] t The scaling factor used for the operation.
 * @param[in] S The input matrix to be multiplied by the transpose of Vt.
 * @param[in] Vt The transpose of the matrix V.
 * @param[out] A The result matrix A.
 * @param[out] B The result matrix B.
 * @param[out] SA The result matrix SA, after applying scaling factors.
 * @param[out] SB The result matrix SB, after applying scaling factors.
 */
void Pneb::mm_SCtimesVtrans(const int mb, const double t, double *S, double *Vt,
                            double *A, double *B, double *SA, double *SB) 
{
   nwpw_timing_function ftimer(19);

   int ms1,ms2,ishift2,ishift1,nj;
   if (mb == -1) 
   {
      ms1 = 0;
      ms2 = ispin;
      ishift2 = ne[0] * ne[0];
      ishift1 = ne[0];
      nj = ne[0] + ne[1];
   } 
   else 
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift2 = 0;
      ishift1 = 0;
      nj = ne[mb];
   }
   for (auto j=0; j<nj; ++j) 
   {
      SA[j] = cos(S[j] * t);
      SB[j] = sin(S[j] * t);
   }
   for (auto ms=ms1; ms<ms2; ++ms) 
   {
      auto shift1 = ms * ishift1;
      auto shift2 = ms * ishift2;
      for (auto k=0; k<ne[ms]; ++k) 
      {
         auto indx1 = shift1;
         auto indx2 = shift2 + k * ne[ms];
         for (auto j=0; j<ne[ms]; ++j) 
         {
            A[indx2] = SA[indx1] * Vt[indx2];
            B[indx2] = SB[indx1] * Vt[indx2];
            ++indx1;
            ++indx2;
         }
      }
   }
}


/*************************************
 *                                   *
 *      Pneb::mm_SCtimesVtrans2      *
 *                                   *
 *************************************/
/**
 * @brief Multiply a matrix S by the transpose of a matrix Vt, with scaling.
 *
 * This function multiplies a matrix S by the transpose of a matrix Vt while applying scaling factors based on the parameters 'mb' and 't'. If 'mb' is -1, the operation is performed for all spin channels; otherwise, it is performed for the specified spin channel.
 *
 * @param[in] mb The spin channel index to which the operation is applied. If -1, the operation is applied to all spin channels.
 * @param[in] t The scaling factor used for the operation.
 * @param[in] S The input matrix to be multiplied by the transpose of Vt.
 * @param[in] Vt The transpose of the matrix V.
 * @param[out] A The result matrix A.
 * @param[out] B The result matrix B.
 * @param[out] SA The result matrix SA, after applying scaling factors.
 * @param[out] SB The result matrix SB, after applying scaling factors.
 */
void Pneb::mm_SCtimesVtrans2(const int mb, const double t, double *S,
                             double *Vt, double *A, double *B, double *SA, double *SB) 
{
   nwpw_timing_function ftimer(19);
   int ms, n, ms1, ms2, ishift2, shift2, ishift1, shift1, nj;
   int j, k, indx1, indx2;
   if (mb == -1) 
   {
      ms1 = 0;
      ms2 = ispin;
      ishift2 = ne[0] * ne[0];
      ishift1 = ne[0];
      nj = ne[0] + ne[1];
   } 
   else 
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift2 = 0;
      ishift1 = 0;
      nj = ne[mb];
   }
   for (j = 0; j < nj; ++j) 
   {
      SA[j] = S[j] * sin(S[j] * t);
      SB[j] = S[j] * cos(S[j] * t);
   }
   for (ms = ms1; ms < ms2; ++ms) 
   {   
      shift1 = ms * ishift1;
      shift2 = ms * ishift2;
      for (k = 0; k < ne[ms]; ++k) 
      {
         indx1 = shift1;
         indx2 = shift2 + k * ne[ms];
         for (j = 0; j < ne[ms]; ++j) 
         {
            A[indx2] = SA[indx1] * Vt[indx2];
            B[indx2] = SB[indx1] * Vt[indx2];
            ++indx1;
            ++indx2;
         }
      }
   }
}


/*************************************
 *                                   *
 *      Pneb::mm_SCtimesVtrans3      *
 *                                   *
 *************************************/
/**
 * @brief Multiply a matrix S by the transpose of a matrix Vt with scaling.
 *
 * This function multiplies a matrix S by the transpose of a matrix Vt while applying scaling factors based on the parameters 'mb' and 't'. If 'mb' is -1, the operation is performed for all spin channels; otherwise, it is performed for the specified spin channel.
 *
 * @param[in] mb The spin channel index to which the operation is applied. If -1, the operation is applied to all spin channels.
 * @param[in] t The scaling factor used for the operation.
 * @param[in] S The input matrix to be multiplied by the transpose of Vt.
 * @param[in] Vt The transpose of the matrix V.
 * @param[out] A The result matrix A.
 * @param[out] B The result matrix B.
 * @param[out] SA The result matrix SA, after applying scaling factors.
 * @param[out] SB The result matrix SB, after applying scaling factors.
 */
void Pneb::mm_SCtimesVtrans3(const int mb, const double t, double *S,
                             double *Vt, double *A, double *B, double *SA, double *SB) 
{
   nwpw_timing_function ftimer(19);
   int ms, n, ms1, ms2, ishift2, shift2, ishift1, shift1, nj;
   int j, k, indx1, indx2;
   if (mb == -1) 
   {
      ms1 = 0;
      ms2 = ispin;
      ishift2 = ne[0] * ne[0];
      ishift1 = ne[0];
      nj = ne[0] + ne[1];
   } 
   else 
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift2 = 0;
      ishift1 = 0;
      nj = ne[mb];
   }
   for (j = 0; j < nj; ++j) 
   {
      SA[j] = sin(S[j] * t);
      SB[j] = 1.0 - cos(S[j] * t);
   }
   for (ms = ms1; ms < ms2; ++ms) 
   {
      shift1 = ms * ishift1;
      shift2 = ms * ishift2;
      for (k = 0; k < ne[ms]; ++k) 
      {
         indx1 = shift1;
         indx2 = shift2 + k * ne[ms];
         for (j = 0; j < ne[ms]; ++j) 
         {
            A[indx2] = SA[indx1] * Vt[indx2];
            B[indx2] = SB[indx1] * Vt[indx2];
            ++indx1;
            ++indx2;
         }
      }
   }
}


/*************************************
 *                                   *
 *         Pneb::mmm_Multiply        *
 *                                   *
 *************************************/
void Pneb::mmm_Multiply(const int mb, double *a, double *b, double alpha,
                        double *c, double beta) 
{
   nwpw_timing_function ftimer(18);
   int ms, n, ms1, ms2, ishift2, shift2;
   if (mb == -1) 
   {
      ms1 = 0;
      ms2 = ispin;
      ishift2 = ne[0] * ne[0];
   } 
   else 
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift2 = 0;
   }
   for (ms = ms1; ms < ms2; ++ms) 
   {
      n = ne[ms];
      if (n > 0) 
      {
         shift2 = ms * ishift2;
     
         DGEMM_PWDFT((char *)"N", (char *)"N", n, n, n, alpha, a + shift2, n,
                     b + shift2, n, beta, c + shift2, n);
      }
   }
}

/*************************************
 *                                   *
 *         Pneb::mmm_Multiply2       *
 *                                   *
 *************************************/
void Pneb::mmm_Multiply2(const int mb, double *a, double *b, double alpha,
                         double *c, double beta) {
  nwpw_timing_function ftimer(18);
  int ms, n, ms1, ms2, ishift2, shift2;
  if (mb == -1) {
    ms1 = 0;
    ms2 = ispin;
    ishift2 = ne[0] * ne[0];
  } else {
    ms1 = mb;
    ms2 = mb + 1;
    ishift2 = 0;
  }
  for (ms = ms1; ms < ms2; ++ms) {
    n = ne[ms];
    if (n > 0) {
      shift2 = ms * ishift2;

      DGEMM_PWDFT((char *)"T", (char *)"N", n, n, n, alpha, a + shift2, n,
                  b + shift2, n, beta, c + shift2, n);
    }
  }
}

/*************************************
 *                                   *
 *         Pneb::mmm_Multiply3       *
 *                                   *
 *************************************/
void Pneb::mmm_Multiply3(const int mb, double *a, double *b, double alpha,
                         double *c, double beta) {
  nwpw_timing_function ftimer(18);
  int ms, n, ms1, ms2, ishift2, shift2;
  if (mb == -1) {
    ms1 = 0;
    ms2 = ispin;
    ishift2 = ne[0] * ne[0];
  } else {
    ms1 = mb;
    ms2 = mb + 1;
    ishift2 = 0;
  }
  for (ms = ms1; ms < ms2; ++ms) {
    n = ne[ms];
    if (n > 0) {
      shift2 = ms * ishift2;

      DGEMM_PWDFT((char *)"N", (char *)"T", n, n, n, alpha, a + shift2, n,
                  b + shift2, n, beta, c + shift2, n);
    }
  }
}


/********************************
 *                              *
 *       Pneb::mm_transpose     *
 *                              *
 ********************************/
void Pneb::mm_transpose(const int mb, const double *a, double *b) 
{
   int i, j, indx, indxt;
   int ms, n, ms1, ms2, ishift2, shift2;
   if (mb == -1) 
   {
      ms1 = 0;
      ms2 = ispin;
      ishift2 = ne[0] * ne[0];
   } 
   else 
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift2 = 0;
   }
   for (ms = ms1; ms < ms2; ++ms) 
   {
      shift2 = ms * ishift2;
      for (j = 0; j < ne[ms]; ++j)
        for (i = 0; i < ne[ms]; ++i) 
        {
           indx  = i + j * ne[ms] + shift2;
           indxt = j + i * ne[ms] + shift2;
           b[indx] = a[indxt];
        }
   }
}

/********************************
 *                              *
 *   Pneb::mm_Kiril_Btransform  *
 *                              *
 ********************************/
void Pneb::mm_Kiril_Btransform(const int mb, double *a, double *b) {
  int ms1, ms2, ishift2;
  if (mb == -1) {
    ms1 = 0;
    ms2 = ispin;
    ishift2 = ne[0] * ne[0];
  } else {
    ms1 = mb;
    ms2 = mb + 1;
    ishift2 = 0;
  }
  for (auto ms = ms1; ms < ms2; ++ms) {
    int shift2 = ms * ishift2;
    for (auto i = 0; i < ne[ms]; ++i)
      for (auto j = 0; j < i; ++j) {
        int indx = i + j * i + shift2;
        double tmp = 0.5 * (a[indx] + b[indx]);
        a[indx] = tmp;
        b[indx] = tmp;
      }
  }
}

/********************************
 *                              *
 *     Pneb::ggm_lambda         *
 *                              *
 ********************************/

#define ITERLMD 220
#define CONVGLMD 1e-15
#define CONVGLMD2 1e-12

// Lagrange multiplier (expensive method)
void Pneb::ggm_lambda(double dte, double *psi1, double *psi2, double *lmbda) 
{
   nwpw_timing_function ftimer(3);
 
   int one = 1;
   double rmone = -1.0;
   double adiff = 0.0;
 
   for (int ms = 0; ms < ispin; ++ms) {
 
     int nn = m_size(ms);
 
     ffm3_sym_Multiply(ms, psi1, psi2, s11, s21, s22);
     m_scale_s22_s21_s11(ms, dte, s22, s21, s11);
 
     int jj;
     int ii = 0;
     int done = 0;
 
     // DCOPY_PWDFT(nn, s21, one, s12, one);
     // DCOPY_PWDFT(nn, s22, one, sa0, one);
     std::memcpy(s12, s21, nn * sizeof(double));
     std::memcpy(sa0, s22, nn * sizeof(double));
 
     while ((!done) && ((ii++) < ITERLMD)) {
       // DCOPY_PWDFT(nn, s22, one, sa1, one);
       std::memcpy(sa1, s22, nn * sizeof(double));
 
       // mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
       // mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
       // mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
       // mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
       d3db::mygdevice.MM6_dgemm(ne[ms], s12, s12, s11, sa0, sa1, st1);
 
       // DCOPY_PWDFT(nn, sa1, one, st1, one);
       std::memcpy(st1, sa1, nn * sizeof(double));
       DAXPY_PWDFT(nn, rmone, sa0, one, st1, one);
       jj = IDAMAX_PWDFT(nn, st1, one);
       adiff = fabs(st1[IDAMAX_PWDFT(nn, st1, one) - 1]);
 
       if (adiff < CONVGLMD)
         done = 1;
       else
         std::memcpy(sa0, sa1,
                     nn *
                         sizeof(double)); // DCOPY_PWDFT(nn, sa1, one, sa0, one);
     }
     // printf("ierr=10 check nn=%d jj=%d adiff=%le ii=%d
     // done=%d\n",nn,jj,adiff,ii,done);
 
     if (adiff > CONVGLMD2) {
       if (!done)
         printf("ierr=10 adiff=%le\n", adiff);
     }
 
     // DCOPY_PWDFT(nn, sa1, one, &lmbda[ms*ne[0]*ne[0]], one);
     std::memcpy(&lmbda[ms*ne[0]*ne[0]],sa1,nn*sizeof(double));
     // std::memcpy(lmbda+ms*ne[0]*ne[0],sa1,nn*sizeof(double));
 
   } // for loop - ms
 
   /* correction due to contraint */
   fmf_Multiply(-1, psi1, lmbda, dte, psi2, 1.0);
}

/********************************
 *                              *
 *     Pneb::ggm_lambda_sic     *
 *                              *
 ********************************/
// Lagrange multiplier for Stiefel and  SIC and(expensive method)
void Pneb::ggm_lambda_sic(double dte, double *psi1, double *psi2,
                          double *lmbda) {
  nwpw_timing_function ftimer(3);

  int one = 1;
  double rmone = -1.0;
  double adiff = 0.0;

  for (int ms = 0; ms < ispin; ++ms) {

    int nn = m_size(ms);

    ffm4_sym_Multiply(ms, psi1, psi2, s11, s21, s12, s22);
    mm_Kiril_Btransform(ms, s12, s21);

    m_scale_s22_s21_s12_s11(ms, dte, s22, s21, s12, s11);

    int jj;
    int ii = 0;
    int done = 0;

    // DCOPY_PWDFT(nn, s21, one, s12, one);
    // DCOPY_PWDFT(nn, s22, one, sa0, one);
    std::memcpy(s12, s21, nn * sizeof(double));
    std::memcpy(sa0, s22, nn * sizeof(double));

    while ((!done) && ((ii++) < ITERLMD)) {
      // DCOPY_PWDFT(nn, s22, one, sa1, one);
      std::memcpy(sa1, s22, nn * sizeof(double));

      // mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
      // mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
      // mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
      // mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
      d3db::mygdevice.MM6_dgemm(ne[ms], s12, s12, s11, sa0, sa1, st1);

      // DCOPY_PWDFT(nn, sa1, one, st1, one);
      std::memcpy(st1, sa1, nn * sizeof(double));
      DAXPY_PWDFT(nn, rmone, sa0, one, st1, one);
      jj = IDAMAX_PWDFT(nn, st1, one);
      adiff = fabs(st1[IDAMAX_PWDFT(nn, st1, one) - 1]);

      if (adiff < CONVGLMD)
        done = 1;
      else
        std::memcpy(sa0, sa1,
                    nn *
                        sizeof(double)); // DCOPY_PWDFT(nn, sa1, one, sa0, one);
    }
    // printf("ierr=10 check nn=%d jj=%d adiff=%le ii=%d
    // done=%d\n",nn,jj,adiff,ii,done);

    if (adiff > CONVGLMD2) {
      if (!done)
        printf("ierr=10 adiff=%le\n", adiff);
    }

    // DCOPY_PWDFT(nn, sa1, one, &lmbda[ms*ne[0]*ne[0]], one);
    std::memcpy(&lmbda[ms * ne[0] * ne[0]], sa1, nn * sizeof(double));
    // std::memcpy(lmbda+ms*ne[0]*ne[0],sa1,nn*sizeof(double));

  } // for loop - ms

  /* correction due to contraint */
  fmf_Multiply(-1, psi1, lmbda, dte, psi2, 1.0);
}

/********************************
 *                              *
 *     Pneb::gen_Ba_Bs          *
 *                              *
 ********************************/
void Pneb::gen_Ba_Bs(const int ms, double *B, double *Bs, double *Ba)
{
   int one = 1;
   double rone = 1.0;
   double rhalf = 0.5;
   double rmone = -1.0;
   double rmhalf = -0.5;
   int nn = m_size(ms);

   Pneb::mm_transpose(ms,B,Ba);

   std::memcpy(Bs,B,nn*sizeof(double));
   DAXPY_PWDFT(nn,rone,Ba,one,Bs,one);
   DSCAL_PWDFT(nn,rhalf,Bs,one);

   DAXPY_PWDFT(nn,rmone,B,one,Ba,one);
   DSCAL_PWDFT(nn,rmhalf,Ba,one);
}

/********************************
 *                              *
 *     Pneb::gen_UD             *
 *                              *
 ********************************/
void Pneb::gen_UD(const int ms, double *Bs, double *D)
{
    m_diagonalize(ne[ms],Bs,D);
    //m_diagonalize(Bs,D);
}

/********************************
 *                              *
 *     Pneb::gen_X              *
 *                              *
 ********************************/
void Pneb::gen_X(const int ms, double *X1, double *tmp, double *A, double *Ba, double *C,
                 double *U, double *D, double *fnm, double *fweight, bool *failed)
{

   // local variables
   int itrlmd = 120;
   double convg = 1.0e-15;

   int one = 1;
   double rzero = 0.0;
   double rone = 1.0;
   double rhalf = 0.5;
   double rmone = -1.0;
   double rmhalf = -0.5;
   int nn = m_size(ms);

   // A = I-A 
   DSCAL_PWDFT(nn,rmone,A,one);
   m_eye(ms,fnm,rone);
   DAXPY_PWDFT(nn,rone,fnm,one,A,one);

   // fnm = I-A 
   std::memcpy(fnm,A,nn*sizeof(double));

   // solve U*D*Ut*X + X*U*D*Ut = fnm for X
   fnm_to_X(ms,fnm,U,D,fweight,tmp);
   std::memcpy(X1,fnm,nn*sizeof(double));


   int it  = 0;
   *failed = true;
   while (*failed &&  (it<itrlmd))
   {
      ++it;

      // fnm = X*C*X
      mmm_Multiply(ms,C,X1,rone,tmp,rzero);
      mmm_Multiply(ms,X1,tmp,rone,fnm,rzero);

      // fnm = Ba*X - X*C*X 
      mmm_Multiply(ms,Ba,X1,rone,fnm,rmone);


      // fnm = Ba*X - X*Ba - X*C*X
      mmm_Multiply(ms,X1,Ba,rmone,fnm,rone);

      // fnm = I-A + Ba*X - X*Ba - X*C*X 
      DAXPY_PWDFT(nn,rone,A,one,fnm,one);


      // solve U*D*Ut*X + X*U*D*Ut = fnm for X
      fnm_to_X(ms,fnm,U,D,fweight,tmp);

      std::memcpy(tmp,X1,nn*sizeof(double));
      DAXPY_PWDFT(nn,rmone,fnm,one,tmp,one);
      double adiff = fabs(tmp[IDAMAX_PWDFT(nn, tmp, one) - 1]);
      std::memcpy(X1,fnm,nn*sizeof(double));

      if (adiff<convg) *failed = false;
   }

}

void Pneb::printNNMatrix(const std::string& str, const int ms, double *A)
{
   int nn = m_size(ms);
      
   int n = ne[ms];
   std::cout << str << "=" << std::endl;
   for (int i=0; i<n; ++i) 
   {
      for (int j= 0; j<n; ++j) 
         std::cout << Efmt(12,4) << A[i*n + j] << " ";
       std::cout << std::endl;
    }
    std::cout << std::endl;
}


/********************************
 *                              *
 *     Pneb::fnm_to_X           *
 *                              *
 ********************************/
void Pneb::fnm_to_X(const int ms, double *fnm, double *U, double *D, double *fweight, double *tmp)
{
   double rzero = 0.0;
   double rone = 1.0;
   
   // fnm = Ut*fnm*U
   mmm_Multiply(ms,fnm,U,rone,tmp,rzero);
   mmm_Multiply2(ms,U,tmp,rone,fnm,rzero);

   // fnm = (Ut*fnm*U)_nm/(d_n+d_m) 
   m_HmldivideDplusD(ms,fnm,D);


   // fnm = X = U*{(Ut*fnm*U)_nm/(d_n+d_m)}*Ut 
   mmm_Multiply(ms,U,fnm,rone,tmp,rzero);
   mmm_Multiply3(ms,tmp,U,rone,fnm,rzero);

   m_Hmlfweightscale(ms,fnm,fweight);
}

/********************************
 *                              *
 *     Pneb::m_eye              *
 *                              *
 ********************************/
void Pneb::m_eye(const int mb, double *A, double alpha)
{
   int ms1, ms2, ishift2;
   if (mb == -1)
   {
      ms1 = 0;
      ms2 = ispin;
      ishift2 = ne[0]*ne[0];
   }
   else
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift2 = 0;
   }
   for (auto ms=ms1; ms<ms2; ++ms)
   {
      int shift2 = ms*ishift2;
      int n = ne[ms];
      for (auto k=0; k<n; ++k)
         A[shift2+k+k*n] = alpha;
   }
}

/********************************
 *                              *
 *   Pneb::m_HmldivideDplusD    *
 *                              *
 ********************************/
void Pneb::m_HmldivideDplusD(const int mb, double *A, double *D)
{
   int ms1, ms2, ishift1, ishift2;
   if (mb == -1)
   {
      ms1 = 0;
      ms2 = ispin;
      ishift1 = ne[0];
      ishift2 = ne[0]*ne[0];
   }
   else
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift1 = 0;
      ishift2 = 0;
   }
   for (auto ms=ms1; ms<ms2; ++ms)
   {
      int shift1 = ms*ishift1;
      int shift2 = ms*ishift2;
      int n = ne[ms];
      for (auto k=0; k<n; ++k)
      for (auto j=0; j<n; ++j)
         A[shift2+j+k*n] /= (D[shift1+k] + D[shift1+j]);
   }
}

/********************************
 *                              *
 *   Pneb::m_Hmlfweightscale    *
 *                              *
 ********************************/
void Pneb::m_Hmlfweightscale(const int mb, double *A, double *fweight)
{
   int ms1, ms2, ishift1, ishift2;
   if (mb == -1)
   {
      ms1 = 0;
      ms2 = ispin;
      ishift1 = ne[0];
      ishift2 = ne[0]*ne[0];
   }
   else
   {
      ms1 = mb;
      ms2 = mb + 1;
      ishift1 = 0;
      ishift2 = 0;
   }
   for (auto ms=ms1; ms<ms2; ++ms)
   {
      int shift1 = ms*ishift1;
      int shift2 = ms*ishift2;
      int n = ne[ms];
      for (auto k=0; k<n; ++k)
      for (auto j=0; j<n; ++j)
         if ((fweight[shift1+j] + fweight[shift1+k]) > 1.0e-9)
            A[shift2+j+k*n] *= (2.0*fweight[shift1+k]/(fweight[shift1+k]+fweight[shift1+j]));
   }
}



/********************************
 *                              *
 *     Pneb::ggm_occ_lambda     *
 *                              *
 ********************************/

#define ITERLMD 220
#define CONVGLMD 1e-15
#define CONVGLMD2 1e-12

// Lagrange multiplier (expensive method)
void Pneb::ggm_occ_lambda(double dte, double *psi1, double *psi2, double *occ, double *lmbda)
{
   nwpw_timing_function ftimer(3);

   int one = 1;
   double rmone = -1.0;
   double adiff = 0.0;
   double overdte = 1.0/dte;



   double *A = s22;   // 0*nn
   double *B = s21;   // 1*nn 
   double *C = s12;   // 2*nn

   double *Ba = s11;  // 3*nn
   double *Bs = sa1;  // 4*nn
   double *fnm = sa0; // 5*nn
   //st1 = st1;       // 6*nn
   //st2 = st2        // 7*nn
   double *D = B;     // 1*nn
   double *U = Bs;    // 4*nn

   bool failed;


   for (int ms=0; ms<ispin; ++ms) 
   {
      int nn = m_size(ms);

      ffm3_Fulls21_sym_Multiply(ms, psi1, psi2, C, B, A);

      gen_Ba_Bs(ms,B,Bs,Ba);
      gen_UD(ms,Bs,D);
      gen_X(ms,st1,st2,A,Ba,C,U,D,fnm,occ,&failed);

      if (failed) 
      {
         std::cout << "Warning: ierr=10, Lagrange Multiplier generation failed." << std::endl 
                   << "        +Try using a smaller time step" << std::endl
                   << "        +Gram-Schmidt being performed, spin: " << ms << std::endl;
         g_ortho(ms,psi2);
      }
      else
      {
         fmf_Multiply(ms,psi1,st1,1.0,psi2,1.0);

         DSCAL_PWDFT(nn,overdte,st1,one);
         std::memcpy(lmbda + ms*ne[0]*ne[0], st1, nn*sizeof(double));
         //Dneall_mm_Expand(ms,st1,lmbda)
      }

   } // for loop - ms

}



/********************************
 *                              *
 *        Pneb::g_ortho         *
 *                              *
 ********************************/
/*
   Performs a Gram-Schmidt orthogonalization on psi
*/
void Pneb::g_ortho(const int mb, double *psi) 
{
   int indxj, indxk, ishift;
   double w;

   int ms1, ms2;
   if (mb == -1)
   {
      ms1 = 0;
      ms2 = ispin;
   }
   else
   {
      ms1 = mb;
      ms2 = mb + 1;
   }


   if (parallelized) 
   {
      //std::ostringstream msg;
      //msg << "NWPW Error: g_ortho() parallelized is NOT supported\n"
      //    << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
      //throw(std::runtime_error(msg.str()));

      int np_j = d1db::parall->np_j();
      int taskid_j = d1db::parall->taskid_j();
      int npack1 = 2*PGrid::npack(1);
      double *tmp = new (std::nothrow) double[npack1]();

      for (auto ms=ms1; ms<ms2; ++ms)
      {
         auto shift0 = ms*neq[0]*npack1;
         auto kcur = np_j-1;
         auto kk   = na[ms][kcur]-1;

         for (auto k=ne[ms]-1; k>=0; --k)
         {
            if (kcur==taskid_j)
            {
               indxk = npack1*kk + shift0;
               w = PGrid::cc_pack_dot(1, psi+indxk, psi+indxk);
               w = 1.0/std::sqrt(w);
               PGrid::c_pack_SMul(1, w, psi+indxk);
               std::memcpy(tmp,psi+indxk,npack1*sizeof(double));
            }

            if (kcur>0)
                d1db::parall->Brdcst_Values(2,kcur,npack1,tmp);

            //*** set j = k+1 ***
            auto jj   = kk;
            auto jcur = kcur;

            --jj;
            if (jj<0) 
            {
               --jcur;
               jj = na[ms][jcur] - 1;
            }

            for (auto j=k-1; j>=0; --j)
            {
               if (jcur==taskid_j)
               {
                  indxj = npack1*jj + shift0;
                  w = -PGrid::cc_pack_dot(1, tmp, psi+indxj);
                  PGrid::cc_pack_daxpy(1, w, tmp, psi+indxj);

               }

               --jj;
               if (jj<0) 
               {
                  --jcur;
                  jj = na[ms][jcur] - 1;
               }
            }

            --kk;
            if (kk<0) 
            {
               --kcur;
               kk = na[ms][kcur] - 1;
            }
         }
      }

      delete[] tmp;
   } 
   else 
   {
      // npj==1
      for (auto ms=ms1; ms<ms2; ++ms) 
      {
         ishift = ms*ne[0]*2*PGrid::npack(1);
         for (auto k=ne[ms]-1; k>=0; --k) 
         {
            indxk = 2*PGrid::npack(1)*k + ishift;
            w = PGrid::cc_pack_dot(1, psi+indxk, psi+indxk);
            w = 1.0/std::sqrt(w);
            PGrid::c_pack_SMul(1, w, psi+indxk);
           
            for (auto j=k-1; j>=0; --j) 
            {
               indxj = 2 * PGrid::npack(1) * j + ishift;
               w = -PGrid::cc_pack_dot(1, psi+indxk, psi+indxj);
               PGrid::cc_pack_daxpy(1, w, psi+indxk, psi+indxj);
            }
         }
      }
   }
   //if (d1db::parall->is_master())
   //   std::cout << std::endl;
}

/********************************
 *                              *
 *    Pneb::g_ortho_excited     *
 *                              *
 ********************************/
/*
   Performs a Gram-Schmidt orthogonalization on psi
*/
void Pneb::g_ortho_excited(const int mb, double *psi, const int nex[], double *psi_excited)
{
   int ms1, ms2;
   if (mb == -1)
   {           
      ms1 = 0; 
      ms2 = ispin;
   }              
   else           
   {
      ms1 = mb;
      ms2 = mb + 1;
   }           

   for (auto ms=ms1; ms<ms2; ++ms)
   {
      // project out filled and virtual spaces
      int kshift = ms*nex[0]*2*PGrid::npack(1);
      for (auto k=0; k<nex[ms]; ++k)
      {
         int indxk = 2*PGrid::npack(1)*k + kshift;

         // project out filled space
         Pneb::g_project_out_filled(psi, ms, psi_excited + indxk);

         // project out lower virtual space
         Pneb::g_project_out_virtual(ms, nex, k, psi_excited, psi_excited+indxk);

         // normalize 
         Pneb::g_norm(psi_excited + indxk);
      }
   }

   // project out virtual space
}

/********************************
 *                              *
 *  Pneb::g_project_out_filled  *
 *                              *
 ********************************/
/**
 * @brief Projects out the filled (occupied) component of a wavefunction from an excited state.
 *
 * This function iteratively projects out the filled (occupied) components from the `psi_excited`
 * wavefunction using the provided `psi` vector. The projection is performed for each value 
 * of `n` from `0` to `ne[ms] - 1`. The function modifies the `psi_excited` wavefunction in place.
 *
 * @param psi           Pointer to the array containing the filled (occupied) components of the wavefunctions.
 * @param ms            The multiplicity of the state, which is used to calculate the shift in the index.
 * @param psi_excited   Pointer to the array containing the excited state wavefunction to be modified.
 *
 * @note This function assumes that the input arrays are stored in a format compatible 
 *       with the `PGrid` utility functions `cc_pack_dot` and `cc_pack_daxpy`.
 *       The `PGrid::npack(1)` function is used to determine the packing size.
 *
 * @warning The function prints the values of `n` and `w` to `std::cout` for each iteration, 
 *          which may result in a significant amount of output if `ne[ms]` is large.
 */
void Pneb::g_project_out_filled(double *psi, const int ms, double *Horb) 
{
   int ishift = ms*ne[0]*2*PGrid::npack(1);
   for (auto n=0; n<ne[ms]; ++n)
   {
     int indx = 2*PGrid::npack(1)*n + ishift;
     double w = -PGrid::cc_pack_dot(1,psi+indx,Horb);
     //std::cout << "   n=" << n << " w=" << w << std::endl;
     PGrid::cc_pack_daxpy(1,w,psi+indx,Horb);
   }
}

/**************************************
 *                                    *
 *  Pneb::g_project_out_filled_below  *
 *                                    *
 **************************************/
/**
 * @brief Projects out components below a certain index in the provided psi array.
 *
 * This function modifies the Horb array by removing contributions from components 
 * of the psi array that are indexed below the specified value of `k`. It loops over 
 * the elements of `psi` from `k-1` down to `0`, computes the dot product of the 
 * relevant part of `psi` and `Horb`, and applies a scaled update to `Horb` using 
 * daxpy operations. 
 *
 * @param psi  Pointer to an array of doubles, representing the wave function or state vector.
 * @param ms   Integer index, representing a state or iteration.
 * @param k    Integer index, specifying the cutoff for the projection.
 * @param Horb Pointer to an array of doubles, representing orbital data or coefficients to be updated.
 *
 * @details 
 * The function utilizes the PGrid class's `npack`, `cc_pack_dot`, and `cc_pack_daxpy` 
 * methods to perform packed dot product and daxpy operations efficiently.
 * 
 * The `psi` array is accessed via the computed `indx` values, which incorporate 
 * both the loop index `km` and the `ms` shift.
 *
 * @note 
 * This function assumes that the PGrid class provides methods for handling packed 
 * grid data and that `ne[0]` is globally accessible within this class or namespace.
 */
void Pneb::g_project_out_filled_below(double *psi, const int ms, const int k, double *Horb) 
{
   int ishift = ms*ne[0]*2*PGrid::npack(1);
   for (auto km=k-1; km>=0; --km)
   {
     int indx = 2*PGrid::npack(1)*km + ishift;
     double w = -PGrid::cc_pack_dot(1,psi+indx,Horb);
     PGrid::cc_pack_daxpy(1,w,psi+indx,Horb);
   }
}


/**************************************
 *                                    *
 *  Pneb::g_project_out_filled_above  *
 *                                    *
 **************************************/
void Pneb::g_project_out_filled_above(double *psi, const int ms, const int k, double *Horb) 
{
   int ishift = ms*neq[0]*2*PGrid::npack(1);
   for (auto ka=k+1; ka<neq[ms]; ++ka)
   {
     int indx = 2*PGrid::npack(1)*ka + ishift;
     double w = -PGrid::cc_pack_dot(1,psi+indx,Horb);
     PGrid::cc_pack_daxpy(1,w,psi+indx,Horb);
   }
}


/*********************************
 *                               *
 *  Pneb::g_project_out_virtual  *
 *                               *
 *********************************/
/**
 * @brief Projects out the virtual component of a wavefunction from an excited state.
 *
 * This function iteratively projects out the virtual components from the `psi_excited`
 * wavefunction using the provided `psiv` vector. The projection is performed for 
 * each value of `km` from `k-1` down to `0`. The function modifies the `psi_excited`
 * wavefunction in place.
 *
 * @param ms       The multiplicity of the state, which is used to calculate the shift in the index.
 * @param nex      An array representing the number of excitations or related quantum states.
 * @param k        The current state index.
 * @param psiv     Pointer to the array containing the virtual components of the wavefunctions.
 * @param psi_excited Pointer to the array containing the excited state wavefunction to be modified.
 *
 * @note This function assumes that the input arrays are stored in a format compatible 
 *       with the `PGrid` utility functions `cc_pack_dot` and `cc_pack_daxpy`.
 *       The `PGrid::npack(1)` function is used to determine the packing size.
 *
 * @warning The function prints the values of `km`, `k`, and `wkm` to `std::cout` 
 *          for each iteration, which may result in a significant amount of output 
 *          if `k` is large.
 */
void Pneb::g_project_out_virtual(const int ms, const int nex[], const int k,  double *psiv,  double *Horb) 
{
   int kshift = ms*nex[0]*2*PGrid::npack(1);
   for (auto km=k-1; km>=0; --km)
   {
      int indxkm = 2*PGrid::npack(1)*km + kshift;
      double wkm = -PGrid::cc_pack_dot(1, psiv+indxkm, Horb);
      PGrid::cc_pack_daxpy(1, wkm, psiv+indxkm, Horb);
      //std::cout << "    - km=" << km << " k=" << k << " wkm=" << wkm <<  std::endl;
   }
}


/*********************************
 *                               *
 *         Pneb::g_norm          *
 *                               *
 *********************************/
/**
 * @brief Normalizes the given vector `psi_to_norm` in place.
 *
 * This function calculates the L2 norm of the input vector `psi_to_norm`
 * and scales the vector so that it has a unit norm. If the norm is found
 * to be effectively zero (below a small threshold), the function will issue
 * a warning and skip the normalization to avoid division by zero.
 *
 * @param psi_to_norm Pointer to the vector that needs to be normalized.
 *
 * @note The function assumes that the input vector is stored in a format
 *       compatible with the `PGrid` utility functions `cc_pack_dot` and `c_pack_SMul`.
 *       These utility functions handle operations on vectors across a distributed grid.
 *
 * @warning If the norm of the vector is less than or equal to `1.0e-12`, 
 *          the normalization will be skipped and a warning will be printed to `std::cerr`.
 */
void Pneb::g_norm(double *psi_to_norm)
{  
   // Compute the dot product of psi_to_norm with itself, resulting in the squared norm
   double squared_norm = PGrid::cc_pack_dot(1, psi_to_norm, psi_to_norm);

   // Check if the norm is effectively zero (within a small threshold)
   if (squared_norm <= 1.0e-12)
   {
        std::cerr << "Warning: Norm is too small to normalize the vector. Skipping normalization." << std::endl;
        return;  // Exit the function early if the vector is too small to normalize
   }
    
   // Compute the inverse of the square root of the squared norm (1/norm)
   double norm_inverse = 1.0 / std::sqrt(squared_norm);
    
   // Scale psi_to_norm by the computed norm_inverse to normalize it
   PGrid::c_pack_SMul(1, norm_inverse, psi_to_norm);
}


/********************************
 *                              *
 *        Pneb::fm_QR           *
 *                              *
 ********************************/
/*
   Performs a modified Gram-Schmidt QR.
*/
void Pneb::fm_QR(const int mb, double *Q, double *R) {
  int n, ms1, ms2, ishift2, shift2;
  int ishift, indxj, indxk, indxm;
  double w;
  if (mb == -1) {
    ms1 = 0;
    ms2 = ispin;
    ishift2 = ne[0] * ne[0];
    std::memset(R, 0, (ne[0] * ne[0] + ne[1] * ne[1]) * sizeof(double));
  } else {
    ms1 = mb;
    ms2 = mb + 1;
    ishift2 = 0;
    std::memset(R, 0, (ne[mb] * ne[mb]) * sizeof(double));
  }

  // npj>1
  if (parallelized) {
    std::ostringstream msg;
    msg << "NWPW Error: fm_QR parallelized is NOT supported\n"
        << "\t - " << __FILE__ << " : " << __LINE__ << std::endl;
    throw(std::runtime_error(msg.str()));
  }

  // npj==1
  else {
    for (auto ms = ms1; ms < ms2; ++ms) {
      ishift = ms * ne[0] * 2 * PGrid::npack(1);
      for (auto k = 0; k < ne[ms]; ++k) {
        indxk = 2 * PGrid::npack(1) * k + ishift;
        indxm = k + k * ne[ms] + ms * ishift2;
        w = PGrid::cc_pack_dot(1, Q + indxk, Q + indxk);
        w = sqrt(w);
        R[indxm] = w;
        w = 1.0 / w;
        PGrid::c_pack_SMul(1, w, Q + indxk);

        for (auto j = k + 1; j < ne[ms]; ++j) {
          indxj = 2 * PGrid::npack(1) * j + ishift;
          indxm = k + j * ne[ms] + ms * ishift2;
          w = PGrid::cc_pack_dot(1, Q + indxk, Q + indxj);
          R[indxm] = w;
          PGrid::cc_pack_daxpy(1, -w, Q + indxk, Q + indxj);
        }
      }
    }
  }
}

/********************************
 *                              *
 *     Pneb::mmm4_AR_to_T4      *
 *                              *
 ********************************/
void Pneb::mmm4_AR_to_T4(const int mb, const double *A, const double *R, double *T4) 
{
   int ms1, ms2, ishift1, ishift2;
   if (mb == -1) 
   {
      ms1 = 0;
      ms2 = ispin;
      ishift1 = ne[0] * ne[0];
      ishift2 = 4 * ne[0] * ne[0];
      std::memset(T4, 0, (4 * ne[0] * ne[0] + ne[1] * ne[1]) * sizeof(double));
   } 
   else 
   {
      ms1 = mb;
      ms2 = mb;
      ishift1 = 0;
      ishift2 = 0;
      std::memset(T4, 0, (4 * ne[mb] * ne[mb]) * sizeof(double));
   }

   for (auto ms = ms1; ms < ms2; ++ms) 
   {
      const int n = ne[ms];
      const double *Asub = A + ms * ishift1;
      const double *Rsub = R + ms * ishift1;
      double *Tsub = T4 + ms * ishift2;
     
      //**** copy A to upper-left of T ****
      for (auto j = 0; j < n; ++j)
         for (auto i = 0; i < n; ++i)
            Tsub[i + j * (2 * n)] = Asub[i + j * n];
     
      //**** copy R to lower-left of T ****
      for (auto j = 0; j < n; ++j)
         for (auto i = 0; i < n; ++i)
            Tsub[(i + n) + j * (2 * n)] = Rsub[i + j * n];
     
      //**** copy -R^t to upper-right of T ****
      for (auto j = 0; j < n; ++j)
         for (auto i = 0; i < n; ++i)
            Tsub[i + (j + n) * (2 * n)] = -Rsub[j + i * n];
   }
}

/********************************
 *                              *
 *     Pneb::m4_FactorSkew      *
 *                              *
 ********************************/
void Pneb::m4_FactorSkew(const int mb, double *K4, double *V4, double *W4,
                         double *Sigma) {
  int ms1, ms2, ishift1, ishift2;
  if (mb == -1) {
    ms1 = 0;
    ms2 = ispin;
    ishift1 = 2 * ne[0];
    ishift2 = 4 * ne[0] * ne[0];
  } else {
    ms1 = mb;
    ms2 = mb;
    ishift1 = 0;
    ishift2 = 0;
  }

  for (auto ms = ms1; ms < ms2; ++ms) {
    int n = 2 * ne[ms];
    int shift1 = ms * ishift1;
    int shift2 = ms * ishift2;

    FACTOR_SKEW_PWDFT(n, K4 + shift2, V4 + shift2, W4 + shift2, Sigma + shift1);
  }
}

/********************************
 *                              *
 *    mm4_RotationSkew_sub1     *
 *                              *
 ********************************/
 /*
  Description:
    This function calculates the sine and cosine values of the scaled Sigma values multiplied by 't',
    and stores the results in arrays SA and SB. The operation is performed element-wise for the given size 'n'.

  Parameters:
    - n: Size of the arrays.
    - t: Scaling factor for the rotation angle.
    - Sigma: Pointer to a double array containing input matrix data.
    - SA: Pointer to a double array where the cosine results will be stored.
    - SB: Pointer to a double array where the sine results will be stored.
*/
static void mm4_RotateSkew_sub1(const int n, const double t,
                                const double *Sigma, double *SA, double *SB) 
{
   for (auto i = 0; i < n; ++i) 
   {
      SA[i] = cos(Sigma[i] * t);
      SB[i] = sin(Sigma[i] * t);
   }
}

/********************************
 *                              *
 *    mm4_RotationSkew_sub2     *
 *                              *
 ********************************/
 /*
  Description:
    This function performs a rotation and skew transformation on matrices V and W using the precalculated
    sine and cosine values from arrays SA and SB. The resulting matrices A and B are stored accordingly.
    The operation is performed element-wise for the given size 'n'.

  Parameters:
    - n: Size of the arrays.
    - SA, SB: Pointers to double arrays containing precalculated sine and cosine values.
    - V, W: Pointers to double arrays containing input matrix data.
    - A: Pointer to a double array where the rotated matrix A will be stored.
    - B: Pointer to a double array where the skewed matrix B will be stored.
*/
static void mm4_RotateSkew_sub2(const int n, const double *SA, const double *SB,
                                const double *V, const double *W, double *A,
                                double *B) {
  for (auto j = 0; j < n; ++j)
    for (auto i = 0; i < n; ++i) {
      int indx = i + j * n;
      A[indx] = V[indx] * SA[j] + W[indx] * SB[j];
      B[indx] = W[indx] * SA[j] - V[indx] * SB[j];
    }
}

/********************************
 *                              *
 *    Pneb::m4_RotationSkew     *
 *                              *
 ********************************/
 /*
  Description:
    This function performs a rotation and skew transformation on given matrices V4 and W4, storing the
    results in matrices A4, B4, and R4. The transformation is based on the provided parameters and
    matrix properties. The operation is performed within specified matrix blocks determined by 'mb' and
    matrix sizes 'ne'.

  Parameters:
    - mb: Index specifying the block of matrices to operate on. If mb is -1, the entire matrix V4 and W4
          are used for calculations. Otherwise, a subset of matrices is used.
    - t: Rotation angle for the transformation.
    - V4, W4: Pointers to double arrays containing input matrix data.
    - Sigma: Pointer to a double array for storing intermediate transformation values.
    - A4, B4: Pointers to double arrays where the resulting matrices A4 and B4 will be stored.
    - R4: Pointer to a double array where the final transformed matrix R4 will be stored.
*/

void Pneb::m4_RotationSkew(const int mb, const double t, double *V4, double *W4,
                           double *Sigma, double *A4, double *B4, double *R4) {
  int ms1, ms2, ishift1, ishift2, nj;
  if (mb == -1) {
    ms1 = 0;
    ms2 = ispin;
    ishift1 = 2 * ne[0];
    ishift2 = 4 * ne[0] * ne[0];
    nj = 2 * (ne[0] + ne[1]);
  } else {
    ms1 = mb;
    ms2 = mb;
    ishift1 = 0;
    ishift2 = 0;
    nj = 2 * ne[mb];
  }

  double rzero = 0.0;
  double rone = 1.0;
  double SA[nj];
  double SB[nj];

  mm4_RotateSkew_sub1(nj, t, Sigma, SA, SB);

  for (auto ms = ms1; ms < ms2; ++ms) {
    int n = 2 * ne[ms];
    int shift1 = ms * ishift1;
    int shift2 = ms * ishift2;
    mm4_RotateSkew_sub2(n, SA + shift1, SB + shift1, V4 + shift2, W4 + shift2,
                        A4 + shift2, B4 + shift2);

    DGEMM_PWDFT((char *)"N", (char *)"T", n, n, n, rone, V4 + shift2, n,
                A4 + shift2, n, rzero, R4 + shift2, n);
    DGEMM_PWDFT((char *)"N", (char *)"T", n, n, n, rone, W4 + shift2, n,
                B4 + shift2, n, rone, R4 + shift2, n);
  }
}

/********************************
 *                              *
 *     Pneb::m4_R4_to_MN        *
 *                              *
 ********************************/
/*
  Description:
    This function converts a given matrix R4 to two matrices M and N. The matrices M and N are
    extracted from specific subblocks of R4 based on the provided index 'mb' and matrix sizes 'ne'.
    The operation is performed according to the specified parameters and matrix properties.
    
  Parameters:
    - mb: Index specifying the block of matrices to operate on. If mb is -1, the entire matrix R4 is
          used to extract M and N. Otherwise, a subset of matrices is used for calculations.
    - R4: Pointer to a double array containing the input matrix data (R4).
    - M: Pointer to a double array where the resulting matrix M will be stored.
    - N: Pointer to a double array where the resulting matrix N will be stored.
*/
void Pneb::m4_R4_to_MN(const int mb, const double *R4, double *M, double *N) 
{
   int ms1, ms2, ishift1, ishift2;
   if (mb == -1) 
   {
      ms1 = 0;
      ms2 = ispin;
      ishift1 = ne[0] * ne[0];
      ishift2 = 4 * ne[0] * ne[0];
   } 
   else 
   {
      ms1 = mb;
      ms2 = mb;
      ishift1 = 0;
      ishift2 = 0;
   }
 
   for (auto ms = ms1; ms < ms2; ++ms) 
   {
      int n = ne[ms];
      int shift1 = ms * ishift1;
      int shift2 = ms * ishift2;
     
      for (auto j = 0; j < n; ++j)
         for (auto i = 0; i < n; ++i) 
         {
            M[i + j * n] = R4[i + j * (2 * n)];
            N[i + j * n] = R4[(i + n) + j * (2 * n)];
         }
   }
}



/********************************************
 *                                          *
 *        Pneb::m_0define_occupation        *
 *                                          *
 ********************************************
 * Description:
 * This function calculates and defines occupation numbers
 * for electronic states based on eigenvalues and optional
 * Hamiltonian matrix elements. It determines the Fermi
 * level using a bisection method and assigns occupations
 * using a specified smearing scheme.
 *
 * Input Parameters:
 * - initial_alpha (double):
 *     Mixing parameter for updating occupations. If negative,
 *     a default value of 1.0 is used.
 *     Constraint: Must be >= 0.0 if explicitly provided.
 * - use_hml (bool):
 *     Flag to indicate whether to update eigenvalues from the
 *     Hamiltonian matrix (`hml`).
 * - multiplicity (int):
 *     Spin multiplicity of the system.
 *     Constraint: Must be >= 1.
 * - ion_charge (double):
 *     Total charge of the ions in the system.
 * - total_charge (double):
 *     Desired total charge of the system.
 * - eig (double*):
 *     Array of eigenvalues, modified if `use_hml` is true.
 *     Assumption: Must be properly allocated with size matching the
 *     number of states (ne[0] + ne[1]).
 * - hml (double*):
 *     Hamiltonian matrix elements used to update eigenvalues
 *     if `use_hml` is true. Assumption: Properly allocated for the
 *     system's size and structure.
 * - smeartype (int):
 *     Type of smearing to apply. Supported types include:
 *       1 - Fermi-Dirac
 *       2 - Gaussian
 *       4 - Marzari-Vanderbilt
 *     Constraint: Must be 1, 2, or 4.
 * - smearkT (double):
 *     Smearing temperature (related to the broadening of the
 *     Fermi distribution).
 *     Constraint: Must be > 0.
 *
 * Output Parameters:
 * - smearfermi (double*):
 *     Array to store the computed Fermi levels for each spin channel.
 *     Assumption: Must be pre-allocated with at least 2 elements.
 * - smearcorrection (double*):
 *     Pointer to a correction term computed based on the
 *     selected smearing scheme.
 *
 * Input/Output Parameters:
 * - occ (double*):
 *     Array of occupation numbers, updated by the function.
 *     Assumption: Properly allocated for all states.
 *
 * Behavior:
 * - Computes the Fermi level using a robust bisection method, ensuring
 *   charge neutrality.
 * - Updates occupation numbers based on the smearing scheme and applies
 *   mixing controlled by `initial_alpha`.
 * - Handles both spin-polarized and non-spin-polarized systems (based on `ispin`).
 * - Computes smearing corrections for entropy-like terms in the energy.
 *
 * Notes:
 * - The function assumes that `eig`, `hml`, and `occ` arrays are properly
 *   allocated and sized according to the system's states and spins.
 * - Throws a `std::runtime_error` if the Fermi energy cannot be determined.
 * - Smearing corrections are only applied for supported `smeartype` values.
 * - This function is **not thread-safe** if arrays like `occ` or `eig` are
 *   shared across threads without synchronization.
 *
 * Usage Example:
 * ```
 * double eig[100], hml[10000], occ[100], smearfermi[2], smearcorrection = 0.0;
 * Pneb myPneb;
 * myPneb.m_0define_occupation(
 *     -1.0, true, 3, 10.0, 8.0,
 *     eig, hml, occ,
 *     1, 0.01, smearfermi, &smearcorrection
 * );
 * ```
 ********************************************/
void Pneb::m_0define_occupation(const double initial_alpha, const bool use_hml,
                          const int multiplicity,
                          const double ion_charge, const double total_charge,
                          double *eig, double *hml, double *occ,
                          const int smeartype, const double smearkT, double *smearfermi, double *smearcorrection)
{
   double ZZ = ion_charge - total_charge;

   // Initialize smear correction and Fermi levels
   smearfermi[0] = 0.0;
   smearfermi[1] = 0.0;
   *smearcorrection = 0.0;

   // Early exit if charge neutrality is effectively satisfied
   if (std::abs(ZZ) < 1.0e-9) return;

   // Determine alpha (mixing parameter)
   double alpha = (initial_alpha < 0.0) ? 1.0 : initial_alpha;

   // Spin-dependent charge values
   double Z[2] = {0.0, 0.0};
   if (ispin == 2) 
   {
      Z[0] = 0.5 * (ZZ + multiplicity - 1);
      Z[1] = 0.5 * (ZZ - multiplicity + 1);
   } 
   else 
   {
      Z[0] = 0.5 * ZZ;
      Z[1] = 0.0;
   }

   // Loop over spin channels
   for (int ms=0; ms<ispin; ++ms) 
   {
      double elower = 1.0e12, eupper = -1.0e12;

      // Update eigenvalues from HML if required
      if (use_hml) 
      {
         for (int n=0; n<ne[ms]; ++n) 
         {
            int index = n + ms * ne[0];
            eig[index] = hml[n + n * ne[ms] + ms * ne[0] * ne[0]];
         }
      }

      // Find elower and eupper
      for (int n=0; n<ne[ms]; ++n) 
      {
         double e = eig[n + ms * ne[0]];
         elower = std::min(elower, e);
         eupper = std::max(eupper, e);
      }

      // Initialize Zlower and Zupper
      double Zlower = 0.0, Zupper = 0.0;
      for (int n = 0; n<ne[ms]; ++n) 
      {
         double e = eig[n + ms * ne[0]];
         Zlower += util_occupation_distribution(smeartype, (e-elower)/smearkT);
         Zupper += util_occupation_distribution(smeartype, (e-eupper)/smearkT);
      }

      double flower = Zlower - Z[ms];
      double fupper = Zupper - Z[ms];

      // Check for valid Fermi level
      if (flower * fupper >= 0.0) {
          throw std::runtime_error("Fermi energy not found");
      }

      // Bisection method to find Fermi level
      double emid = 0.0, Zmid = 0.0, fmid = 0.0;
      int it = 0;

      while ((std::abs(fmid) > 1.0e-11 || std::abs(eupper - elower) > 1.0e-11) && it < 50) 
      {
         ++it;
         emid = 0.5 * (elower + eupper);
         Zmid = 0.0;
         for (auto n = 0; n < ne[ms]; ++n) 
         {
            double e = eig[n + ms * ne[0]];
            Zmid += util_occupation_distribution(smeartype, (e-emid)/smearkT);
         }
         fmid = Zmid - Z[ms];

         if (fmid < 0.0) 
         {
            flower = fmid;
            elower = emid;
         } 
         else 
         {
            fupper = fmid;
            eupper = emid;
         }
      }

      smearfermi[ms] = emid;

      // Update occupations and calculate smear corrections
      for (auto n=0; n<ne[ms]; ++n) 
      {
         int index = n + ms*ne[0];
         double e = eig[index];
         double x = (e - smearfermi[ms]) / smearkT;
         double f = util_occupation_distribution(smeartype, x);
         double f0 = occ[index];

         // Update occupations with mixing
         occ[index] = (1.0 - alpha)*f0 + alpha*f;

         // Calculate corrections based on smearing type
         if (smeartype == 1) { // Fermi-Dirac correction
             if (occ[index] > 1.0e-6 && (1.0 - occ[index]) > 1.0e-6) {
                 *smearcorrection += smearkT * (occ[index] * log(occ[index]) +
                                                (1.0 - occ[index]) * log(1.0 - occ[index]));
             }
         } else if (smeartype == 2) { // Gaussian correction
             *smearcorrection -= smearkT * exp(-x * x) / (4.0 * sqrt(M_PI));
         } else if (smeartype == 4) { // Marzari-Vanderbilt correction
             *smearcorrection -= smearkT * exp(-(x + sqrt(0.5)) * (x + sqrt(0.5))) *
                                 (1.0 + sqrt(2.0) * x) / (2.0 * sqrt(M_PI));
         }
      }
   }

   // Adjust smear correction for unpolarized systems
   if (ispin==1) 
      *smearcorrection *= 2.0;
}


} // namespace pwdft
