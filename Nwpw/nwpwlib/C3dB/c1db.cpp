/* c1db.cpp
   Author - Eric Bylaska

        this class is used for defining 3d parallel maps
*/

/**
 * @class c1db
 * @brief Container for operations on distributed 3D blocks.
 *
 * The `c1db` class is designed to handle operations related to distributed 3D blocks.
 * It provides functionality for performing matrix multiplications and other mathematical
 * operations on distributed data structures. This class is often used in conjunction
 * with the `Parallel` and `gdevice2` classes to manage parallelism and device operations.
 */

/*
#include        <cmath>
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>

*/

#include <iostream>
#include "Mapping1.hpp"
#include "Parallel.hpp"
#include "blas.h"
#include "gdevice2.hpp"
#include "c1db.hpp"

#include "iofmt.hpp"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

namespace pwdft {

/***********************************
 *                                 *
 *       CMatrix_start_rot         *
 *                                 *
 ***********************************/
/**
 * @brief Start a rotation operation in a distributed matrix.
 *
 * This function initiates a rotation operation for a distributed matrix. It performs
 * asynchronous communication between neighboring processes to exchange matrix data
 * required for the rotation operation. The rotation operation is typically used in
 * matrix-matrix multiplication routines.
 *
 * @param parall A pointer to the Parallel object managing parallelism.
 * @param j The rotation index or step.
 * @param A Pointer to the source matrix data to be rotated.
 * @param W Pointer to the workspace matrix data required for rotation.
 * @param lda Leading dimension of the matrices A and W.
 * @param na Array of dimensions for distributed blocks along the specified axis.
 * @param rqst_indx An index indicating the request index for asynchronous communication.
 *
 * @note This function is used in the context of distributed matrix operations and
 *       should be called in coordination with similar functions for other matrix
 *       dimensions.
 *
 * @see CMatrix_dgemm1, CMatrix_dgemm2, CMatrix_dgemm3
 */
static void CMatrix_start_rot(Parallel *parall,
                             const int j,
                             double *A, double *W, const int lda, const int *na,
                             const int rqst_indx)
{
   int taskid_j = parall->taskid_j();
   int np_j = parall->np_j();

   int proc_to   = (taskid_j+j) % np_j;
   int proc_from = (taskid_j-j+np_j) % np_j;
   int msgtype   = j;
   int amsglen = lda*na[taskid_j];
   int wmsglen = lda*na[proc_from];

   if (wmsglen>0)
       parall->a2dreceive(rqst_indx,msgtype,proc_from,wmsglen,W);

   if (amsglen>0)
      parall->a2dsend(rqst_indx,msgtype,proc_to,amsglen,A);
}


/***********************************
 *                                 *
 *        CMatrix_end_rot          *
 *                                 *
 ***********************************/
/**
 * @brief Complete a rotation operation in a distributed matrix.
 *
 * This function finalizes a rotation operation for a distributed matrix. It ensures
 * that any asynchronous communication initiated by the corresponding
 * `CMatrix_start_rot` call has completed. The rotation operation is typically used in
 * matrix-matrix multiplication routines.
 *
 * @param parall A pointer to the Parallel object managing parallelism.
 * @param request_indx The index of the request array associated with the rotation operation.
 *
 * @note This function is used in the context of distributed matrix operations and
 *       should be called in coordination with similar functions for other matrix
 *       dimensions.
 *
 * @see CMatrix_start_rot, CMatrix_dgemm1, CMatrix_dgemm2, CMatrix_dgemm3
 */
static void CMatrix_end_rot(Parallel *parall, const int request_indx)
{
   parall->awaitall(request_indx);
}


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/**
 * @brief Constructor for the d1db class.
 *
 * This constructor initializes an instance of the d1db class, which is used for
 * distributed 3D block operations. It sets up the necessary data structures and
 * communication channels for distributed matrix operations.
 *
 * @param inparall A pointer to the Parallel object managing parallelism.
 * @param inmaptype The mapping type for distributed operations.
 * @param inispin The initial spin value.
 * @param inne An array specifying the dimensions of the distributed block.
 *
 * @note This constructor also initiates asynchronous communication channels if
 *       the number of processes in the j-dimension (np_j) is greater than 1.
 *       These channels are used in distributed matrix operations.
 */
c1db::c1db(Parallel *inparall, const int inmaptype, const int inispin, int *inne)
     : Mapping1(inmaptype, inparall->np_j(), inparall->taskid_j(), inispin, inne) 
{
  parall = inparall;

  int np_j = parall->np_j();
  if (np_j>1)
  {
     request1_indx = parall->max_reqstat-2;
     request2_indx = parall->max_reqstat-1;
     parall->astart(request1_indx,2*np_j+2);
     parall->astart(request2_indx,2*np_j+2);
  }
}


/***********************************
 *                                 *
 *      c1db::CMatrix_dgemm1       *
 *                                 *
 ***********************************/

 /*
  Description:
    This function performs a distributed matrix multiplication using the DGEMM algorithm with additional parallelism
    through MPI communication. It calculates the product of matrices A and B and accumulates the result into matrix C,
    considering the specified block size and process grid configuration. The function manages data distribution and
    communication to perform the multiplication efficiently.

  Parameters:
    - parall: Parallel object
    - mygdevice: 
    - m, n, k: Dimensions of the matrices involved in the multiplication.
    - nblock: Block size for matrix multiplication.
    - alpha: Scalar constant for the multiplication operation.
    - A: Pointer to the input matrix A.
    - lda: Leading dimension of matrix A.
    - ma, na: Pointers to arrays specifying row and column dimensions of submatrices for each process.
    - B: Pointer to the input matrix B.
    - ldb: Leading dimension of matrix B.
    - mb, nb: Pointers to arrays specifying row and column dimensions of submatrices for each process.
    - beta: Scalar constant for scaling matrix C.
    - C: Pointer to the output matrix C, where the result is accumulated.
    - ldc: Leading dimension of matrix C.
    - mc, nc: Pointers to arrays specifying row and column dimensions of result submatrices for each process.
    - work1: Temporary workspace array.
    - work2: Temporary workspace array.
*/
void c1db::CMatrix_dgemm1(Parallel *parall, 
                          gdevice2 *mygdevice,
                          int m, int n, int k, int nblock,
                          double alpha,
                          double *A, int lda, int *ma, int *na,
                          double *B, int ldb, int *mb, int *nb,
                          double beta,
                          double *C, int ldc, int *mc, int *nc,
                          double  *work1, double *work2)
{
   auto taskid_i = parall->taskid_i();
   auto taskid_j = parall->taskid_j();
   auto np_i = parall->np_i();
   auto np_j = parall->np_j();

   int one=1;
   int nn=(nc[taskid_j]*mc[taskid_i]);
   DSCAL_PWDFT(nn,beta,C,one);

   int iwrk,shift,ierr;
   int ii=0;
   int jj=0;
   int kk=0;
   int icur=0;
   int jcur=0;
   double rone=1.0;

   /* loop over all row pannels of C */
   while (kk<k)
   {
      iwrk = MIN(nblock, mb[icur]-ii);
      iwrk = MIN(iwrk,   na[jcur]-jj);

      /* pack current iwrk columns of A into work1 */
      if (taskid_j==jcur)
         DLACPY_PWDFT((char*) "G",ma[taskid_i],iwrk,
                      A+jj*lda,lda,
                      work1, ma[taskid_i]);

      /* pack current iwrk rows of B into work2 */
      if (taskid_i==icur)
         DLACPY_PWDFT((char*) "G",iwrk,nb[taskid_j],
                      B+ii,ldb,
                      work2, iwrk);


      /* broadcast work1  within my row */
      //ierr = MPI_Bcast(work1,iwrk*ma[taskid_i],MPI_DOUBLE_PRECISION,jcur,comm_j);
      parall->Brdcst_Values(2,jcur,iwrk*ma[taskid_i],work1);


      /* broadcast work2  within my column */
      //ierr = MPI_Bcast(work2,iwrk*nb[taskid_j],MPI_DOUBLE_PRECISION,icur,comm_i);
      parall->Brdcst_Values(1,icur,iwrk*ma[taskid_j],work2);

      if ((iwrk>0) && (mc[taskid_i]>0) &&  (nc[taskid_j]>0))
         mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                              alpha,
                              work1,ma[taskid_i],
                              work2,iwrk,
                              rone,
                              C,ldc);
         /* DGEMM_PWDFT((char *) "N",(char *) "N",mc[taskid_i],nc[taskid_j],iwrk,
                   alpha,
                   work1,ma[taskid_i],
                   work2,iwrk,
                   rone,
                   C,ldc);
         */

      ii += iwrk;
      jj += iwrk;
      kk += iwrk;

      if (jj>=na[jcur])
      {
        ++jcur;
        jj=0;
      }
      if (ii>=mc[icur])
      {
        ++icur;
        ii=0;
      }

   }
}


/***********************************
 *                                 *
 *         CMatrix_dgemm2          *
 *                                 *
 ***********************************/
 /*
  Description:
    This function performs a distributed matrix multiplication with additional parallelism through MPI communication.
    It calculates the product of matrices transpose(A) and B and accumulates the result into matrix C. The operation is block-wise
    and involves managing data distribution and communication for efficient computation.

  Parameters:
    - parall: Parallel object
    - mygdevice: 
    - m, n, k: Dimensions of the matrices involved in multiplication.
    - nblock: Block size for matrix multiplication.
    - alpha: Scalar constant for the multiplication operation.
    - A: Pointer to the input matrix A.
    - lda: Leading dimension of matrix A.
    - ma, na: Pointers to arrays specifying row and column dimensions of submatrices for each process.
    - B: Pointer to the input matrix B.
    - ldb: Leading dimension of matrix B.
    - mb, nb: Pointers to arrays specifying row and column dimensions of submatrices for each process.
    - beta: Scalar constant for scaling matrix C.
    - C: Pointer to the output matrix C where the result is accumulated.
    - ldc: Leading dimension of matrix C.
    - mc, nc: Pointers to arrays specifying row and column dimensions of result submatrices for each process.
    - work1: Temporary workspace array.
    - work2: Temporary workspace array.
*/
void c1db::CMatrix_dgemm2(Parallel *parall,
                          gdevice2 *mygdevice,
                          int m, int n, int k, int nblock,
                          double alpha,
                          double *A, int lda, int *ma, int *na,
                          double *B, int ldb, int *mb, int *nb,
                          double beta,
                          double *C, int ldc, int *mc, int *nc,
                          double  *work1, double *work2)
{
   auto taskid_i = parall->taskid_i();
   auto taskid_j = parall->taskid_j();
   auto np_i = parall->np_i();
   auto np_j = parall->np_j();

   int one=1;
   int nn=(nc[taskid_j]*mc[taskid_i]);
   DSCAL_PWDFT(nn,beta,C,one);

   int iwrk,shift,ierr;
   int ii=0;
   int jj=0;
   int kk=0;
   int icur=0;
   int jcur=0;
   double rone=1.0;
   double rzero=0.0;

   /* loop over all row pannels of C */
   while (kk<m)
   {
      iwrk = MIN(nblock, mc[icur]-ii);
      iwrk = MIN(iwrk,   na[jcur]-jj);

      /* iwrk*nc(taskid_j) submatrix !=0 */
      if (ma[taskid_i]>0)
      {
         /* pack current iwrk columns of A into work1 */
         if (taskid_j==jcur)
         {
            DLACPY_PWDFT((char*) "G",ma[taskid_i],iwrk,
                         A+jj*lda,lda,
                         work1, ma[taskid_i]);
         }

         /* broadcast work1  within my row */
         //ierr = MPI_Bcast(work1,iwrk*ma[taskid_i],MPI_DOUBLE_PRECISION,jcur,comm_j);
         parall->Brdcst_Values(2,jcur,iwrk*ma[taskid_i],work1);

         if ((iwrk>0) && (nb[taskid_j]>0))
            mygdevice->TN_dgemm2(iwrk,nb[taskid_j],ma[taskid_i],
                                 alpha,
                                 work1,ma[taskid_i],
                                 B,ldb,
                                 rzero,
                                 work2,iwrk);
/*
            DGEMM_PWDFT((char *) "T",(char *) "N",iwrk,nb[taskid_j],ma[taskid_i],
                   alpha,
                   work1, ma[taskid_i],
                   B, ldb,
                   rzero,
                   work2, iwrk);
            */
      }

      /* iwrk*nc(taskid_j+1) submatrix ==0 */
      else
         std::memset(work2,0,nc[taskid_j]*iwrk*sizeof(double));

      /* summ to node that holds current rows of C */
      //ierr = MPI_Reduce(work2,work1,nc[taskid_j]*iwrk,MPI_DOUBLE_PRECISION,MPI_SUM,icur,comm_i);
      parall->Reduce_Values(1,icur,nc[taskid_j]*iwrk,work2,work1);

      /* add to current rows of C */
      if (taskid_i==icur)
      {
         shift = 0;
         for (auto i=ii; i<(ii+iwrk); ++ i)
         {
            DAXPY_PWDFT(nc[taskid_j],rone,work1+shift,iwrk,C+i,mc[taskid_i]);
            ++shift;
         }
      }
      ii += iwrk;
      jj += iwrk;
      kk += iwrk;

      if (jj>=na[jcur])
      {
        ++jcur;
        jj=0;
      }
      if (ii>=mc[icur])
      {
        ++icur;
        ii=0;
      }

   }
}


/***********************************
 *                                 *
 *         CMatrix_dgemm3          *
 *                                 *
 ***********************************/
/*
  Description:
    This function performs a distributed matrix multiplication with additional parallelism through MPI communication.
    It calculates the product of matrices A and transpose(B) and accumulates the result into matrix C. The operation is block-wise
    and involves managing data distribution and communication for efficient computation.

  Parameters:
    - m, n, k: Dimensions of the matrices involved in multiplication.
    - nblock: Block size for matrix multiplication.
    - alpha: Scalar constant for the multiplication operation.
    - A: Pointer to the input matrix A.
    - lda: Leading dimension of matrix A.
    - ma, na: Pointers to arrays specifying row and column dimensions of submatrices for each process.
    - B: Pointer to the input matrix B.
    - ldb: Leading dimension of matrix B.
    - mb, nb: Pointers to arrays specifying row and column dimensions of submatrices for each process.
    - beta: Scalar constant for scaling matrix C.
    - C: Pointer to the output matrix C where the result is accumulated.
    - ldc: Leading dimension of matrix C.
    - mc, nc: Pointers to arrays specifying row and column dimensions of result submatrices for each process.
    - taskid_i, taskid_j: Task IDs within the process grid.
    - np_i, np_j: Number of processes in each dimension.
    - comm_i, comm_j: MPI communicators for process rows and columns.
    - work1: Temporary workspace array.
    - work2: Temporary workspace array.
*/
void c1db::CMatrix_dgemm3(Parallel *parall,
                          gdevice2 *mygdevice,
                          int m, int n, int k, int nblock,
                          double alpha,
                          double *A, int lda, int *ma, int *na,
                          double *B, int ldb, int *mb, int *nb,
                          double beta,
                          double *C, int ldc, int *mc, int *nc,
                          double *work1, double *work2)
{
   auto taskid_i = parall->taskid_i();
   auto taskid_j = parall->taskid_j();
   auto np_i = parall->np_i();
   auto np_j = parall->np_j();

   int one=1;
   int nn=(nc[taskid_j]*mc[taskid_i]);
   DSCAL_PWDFT(nn,beta,C,one);

   int iwrk,shift,ierr;
   int ii=0;
   int jj=0;
   int kk=0;
   int icur=0;
   int jcur=0;
   double rone=1.0;
   double rzero=0.0;

   /* loop over all row pannels of C */
   while (kk<n)
   {
      iwrk = MIN(nblock, mb[icur]-ii);
      iwrk = MIN(iwrk,   nc[jcur]-jj);

      /* pack current iwrk rows of B into work2 */
      if (taskid_i==icur)
         DLACPY_PWDFT((char*) "G",iwrk,nb[taskid_j],
                      B+ii,ldb,
                      work2, iwrk);

      /* broadcast work2  within my column */
      //ierr = MPI_Bcast(work2,iwrk*nb[taskid_j],MPI_DOUBLE_PRECISION,icur,comm_i);
      parall->Brdcst_Values(1,icur,iwrk*nb[taskid_j],work2);

      if ((iwrk>0) && (na[taskid_i]>0) &&  (mc[taskid_j]>0))
         mygdevice->NT_dgemm3(mc[taskid_i],iwrk,na[taskid_j],
                                 alpha,
                                 A,lda,
                                 work2,iwrk,
                                 rzero,
                                 work1,mc[taskid_i]);
         /*
         DGEMM_PWDFT((char *) "N",(char *) "T",mc[taskid_i],iwrk,na[taskid_j],
                   alpha,
                   A,lda,
                   work2,iwrk,
                   rzero,
                   work1,mc[taskid_i]);
         */

      /* reduce work1 within my row */
      //ierr = MPI_Reduce(work1,work2,mc[taskid_i]*iwrk*,MPI_DOUBLE_PRECISION,jcur,comm_j);
      parall->Reduce_Values(2,jcur,mc[taskid_j]*iwrk,work1,work2);

      if (taskid_j==jcur)
      {
         shift = 0;
         for (auto j=jj; j<(jj+iwrk); ++j)
         {
            DAXPY_PWDFT(mc[taskid_i],rone,work2+shift,one,C+j*ldc,one);
            shift += mc[taskid_i];
         }
      }

      ii += iwrk;
      jj += iwrk;
      kk += iwrk;

      if (jj>=nc[jcur])
      {
        ++jcur;
        jj=0;
      }
      if (ii>=mb[icur])
      {
        ++icur;
        ii=0;
      }

   }
}


/***********************************
 *                                 *
 *         CMatrix_dgemm2c         *
 *                                 *
 ***********************************/
/**
 * @brief Perform a distributed matrix-matrix multiplication for Block 2 (Column-Cyclic).
 *
 * This function performs a distributed matrix-matrix multiplication for Block 2
 * using a column-cyclic distribution. It uses a parallel device (mygdevice) to
 * perform the computation efficiently. The input matrices A, B, and C are distributed,
 * and the multiplication is performed as part of a parallel computation.
 *
 * @param parall A pointer to the Parallel object managing parallelism.
 * @param mygdevice A pointer to the parallel device for efficient computation.
 * @param m The number of rows in matrix C.
 * @param n The number of columns in matrix C.
 * @param k The inner dimension for the matrix multiplication.
 * @param nblock The block size for distributed matrices.
 * @param A Pointer to matrix A (distributed).
 * @param B Pointer to matrix B (distributed).
 * @param lda Leading dimension of matrix A.
 * @param ma Array specifying the dimensions of the distributed block in matrix A.
 * @param ma1 Array specifying the dimensions of the distributed block in matrix A.
 * @param na Array specifying the dimensions of the distributed block in matrix A.
 * @param C Pointer to matrix C (distributed).
 * @param ldc Leading dimension of matrix C.
 * @param mc Array specifying the dimensions of the distributed block in matrix C.
 * @param nc Array specifying the dimensions of the distributed block in matrix C.
 * @param work1 Temporary storage for intermediate computation.
 * @param work2 Temporary storage for intermediate computation.
 *
 * @note This function is part of a distributed matrix multiplication algorithm
 *       and involves complex data distribution and computation operations.
 *       It uses parallel device (mygdevice) for optimized computation.
 */
void c1db::CMatrix_dgemm2c(Parallel *parall,
                          gdevice2 *mygdevice,
                          int m, int n, int k, int nblock,
                          double *A, double *B, int lda, int *ma, int *ma1, int *na,
                          double *C, int ldc, int *mc, int *nc,
                          double  *work1, double *work2)
{
   auto taskid_i = parall->taskid_i();
   auto taskid_j = parall->taskid_j();
   auto np_i = parall->np_i();
   auto np_j = parall->np_j();

   int n2     = na[taskid_j];
   int npack2 = ma[taskid_i];
   int nida2  = ma1[taskid_i];

   int one=1;
   int nn=(nc[taskid_j]*mc[taskid_i]);
   std::memset(C,0,nn*sizeof(double));

   int iwrk,shift;
   int ii=0;
   int jj=0;
   int kk=0;
   int icur=0;
   int jcur=0;
   double rone=1.0;
   double rzero=0.0;

   while (kk<m)
   {
      iwrk = MIN(nblock, mc[icur]-ii);
      iwrk = MIN(iwrk,   na[jcur]-jj);

      nn=(nc[taskid_j]*iwrk);
      std::memset(work2,0,nn*sizeof(double));

      /* iwrk*nc(taskid_j) submatrix !=0 */
      if (ma[taskid_i]>0)
      {
         /* pack current iwrk columns of A into work1 */
         if (taskid_j==jcur)
         {
            DLACPY_PWDFT((char*) "G",ma[taskid_i],iwrk,
                         A+jj*lda,lda,
                         work1, ma[taskid_i]);
         }

         /* broadcast work1  within my row */
         //ierr = MPI_Bcast(work1,iwrk*ma[taskid_i],MPI_DOUBLE_PRECISION,jcur,comm_j);
         parall->Brdcst_Values(2,jcur,iwrk*ma[taskid_i],work1);

         if ((iwrk>0) && (na[taskid_j]>0))
            mygdevice->TN_dgemm2c(iwrk,n2,npack2,nida2,work1,B,work2);
      }

      /* summ to node that holds current rows of C */
      //ierr = MPI_Reduce(work2,work1,nc[taskid_j]*iwrk,MPI_DOUBLE_PRECISION,MPI_SUM,icur,comm_i);
      parall->Reduce_Values(1,icur,nc[taskid_j]*iwrk,work2,work1);

      /* add to current rows of C */
      if (taskid_i==icur)
      {
         shift = 0;
         for (auto i=ii; i<(ii+iwrk); ++i)
         {
            DAXPY_PWDFT(nc[taskid_j],rone,work1+shift,iwrk,C+i,mc[taskid_i]);
            ++shift;
         }
      }
      ii += iwrk;
      jj += iwrk;
      kk += iwrk;

      if (jj>=na[jcur])
      {
        ++jcur;
        jj=0;
      }
      if (ii>=mc[icur])
      {
        ++icur;
        ii=0;
      }
   }
}

/***********************************
 *                                 *
 *         CMatrix_dgemm1_rot2     *
 *                                 *
 ***********************************/
/**
 * @brief Perform a distributed matrix-matrix multiplication with rotation for Block 2.
 *
 * This function performs a distributed matrix-matrix multiplication with rotation
 * for Block 2. It uses a parallel device (mygdevice) to perform the computation
 * efficiently. The input matrices A, B, and C are distributed, and the multiplication
 * is performed as part of a parallel computation.
 *
 * @param parall A pointer to the Parallel object managing parallelism.
 * @param mygdevice A pointer to the parallel device for efficient computation.
 * @param m The number of rows in matrix C.
 * @param n The number of columns in matrix C.
 * @param k The inner dimension for the matrix multiplication.
 * @param nblock The block size for distributed matrices.
 * @param alpha The scaling factor for the multiplication.
 * @param A Pointer to matrix A (distributed).
 * @param lda Leading dimension of matrix A.
 * @param ma Array specifying the dimensions of the distributed block in matrix A.
 * @param na Array specifying the dimensions of the distributed block in matrix A.
 * @param B Pointer to matrix B (distributed).
 * @param ldb Leading dimension of matrix B.
 * @param mb Array specifying the dimensions of the distributed block in matrix B.
 * @param nb Array specifying the dimensions of the distributed block in matrix B.
 * @param beta The scaling factor for matrix C.
 * @param C Pointer to matrix C (distributed).
 * @param ldc Leading dimension of matrix C.
 * @param mc Array specifying the dimensions of the distributed block in matrix C.
 * @param nc Array specifying the dimensions of the distributed block in matrix C.
 * @param Bcol Pointer to matrix Bcol (unmodified).
 * @param Bwork Temporary storage for matrix B.
 * @param work1 Temporary storage for intermediate computation.
 * @param work2 Temporary storage for intermediate computation.
 *
 * @note This function is part of a distributed matrix multiplication algorithm
 *       and involves complex data distribution and rotation operations.
 *       It uses parallel device (mygdevice) for optimized computation.
 */
void c1db::CMatrix_dgemm1_rot2(Parallel *parall,
                          gdevice2 *mygdevice,
                          int m, int n, int k, int nblock,
                          double alpha,
                          double *A, int lda, int *ma, int *na,
                          double *B, int ldb, int *mb, int *nb,
                          double beta,
                          double *C, int ldc, int *mc, int *nc,
                          double *Bcol, double *Bwork, double  *work1, double *work2)
{



   int taskid = parall->taskid();
   int taskid_i = parall->taskid_i();
   int taskid_j = parall->taskid_j();
   int np_i = parall->np_i();
   int np_j = parall->np_j();

   int jcur=0;
   double rone=1.0;
   double rzero=0.0;

   int iwrk;
   int ii = 0;
   int bshift = 0;
   int bshift2[np_j+1]; bshift2[0] = 0;

/*
   double *Aall = new (std::nothrow) double[4*lda]();
   double *C2 = new (std::nothrow) double[ldc]();
   std::memset(Aall,0,4*lda*sizeof(double));
   std::memcpy(Aall+lda*taskid_j,A,lda*sizeof(double));
   parall->Vector_SumAll(2,4*lda,Aall);
   std::memcpy(C2,C,ldc*sizeof(double));
   */


   //collect B into columns
   for (auto jj=1; jj<=np_j; ++jj)
   {
      bshift     += na[jj-1]*nb[taskid_j];
      bshift2[jj] = bshift;
   }

   int ne0 = 0;
   for (auto i=0; i<np_i; ++i) 
      ne0 += mb[i];

   int j1 = 0;
   for (auto jj=0; jj<taskid_j; ++jj)
      j1 += nb[jj];

   bshift = 0;
   int iwrk2 = nb[taskid_j];
   ii = 0;
   for (jcur=0; jcur<np_j; ++jcur)
   {
      iwrk = na[jcur];
      for (auto j=0; j<iwrk2; ++j)
      for (auto i=0; i<iwrk;  ++i)
      {
         Bwork[bshift+i+j*iwrk] = B[ii+i+(j+j1)*ne0];
      }
      bshift += iwrk*iwrk2;
      ii     += iwrk;
   }

   //C = beta*C
   double tmphml[16];
   std::memset(tmphml,0,16*sizeof(double));

   int one=1;
   int nn=(nc[taskid_j]*mc[taskid_i]);
   DSCAL_PWDFT(nn,beta,C,one);


/*
              mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],4,
                           alpha,
                           Aall,lda,
                           B+4*taskid_j,4,
                           beta,
                           C2,ldc);
              */

   std::memset(work1,0,lda*na[taskid_j]*sizeof(double));
   CMatrix_start_rot(parall,1,A,work1,lda,na,request1_indx);

   jcur = taskid_j;
   iwrk = na[jcur];

   if ((iwrk>0) && (mc[taskid_i]>0) &&  (nc[taskid_j]>0))
   {
      
      mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                           alpha,
                           A,lda,
                           Bwork+bshift2[jcur],iwrk,
                           rone,
                           C,ldc);
     

     /*
      mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                           alpha,
                           Aall+jcur*lda,lda,
                           Bwork+bshift2[jcur],iwrk,
                           rone,
                           C,ldc);
      */

      tmphml[4*taskid_j] =  Bwork[bshift2[jcur]];
   }

   bool jeven = true;
   for (auto j=2; j<(np_j); ++j)
   {
       if (jeven)
       {
          jeven = false;
          jcur = (jcur-1+np_j) % np_j;
          iwrk = na[jcur];
          CMatrix_end_rot(parall,request1_indx);
          std::memset(work2,0,lda*na[taskid_j]*sizeof(double));
          CMatrix_start_rot(parall,j,A,work2,lda,na,request2_indx);
          if ((iwrk>0) && (mc[taskid_i]>0) &&  (nc[taskid_j]>0))
          {
          
             mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                                  alpha,
                                  work1,lda,
                                  Bwork+bshift2[jcur],iwrk,
                                  rone,
                                  C,ldc);
             /*
              mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                           alpha,
                           Aall+jcur*lda,lda,
                           Bwork+bshift2[jcur],iwrk,
                           rone,
                           C,ldc);
              */

              tmphml[4*taskid_j+1] =  Bwork[bshift2[jcur]];
          }
       }
       else
       {
          jeven = true;
          jcur = (jcur-1+np_j) % np_j;
          iwrk = na[jcur];
          CMatrix_end_rot(parall,request2_indx);
          std::memset(work1,0,lda*na[taskid_j]*sizeof(double));
          CMatrix_start_rot(parall,j,A,work1,lda,na,request1_indx);
          if ((iwrk>0) && (mc[taskid_i]>0) &&  (nc[taskid_j]>0))
          {
           
             mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                                  alpha,
                                  work2,lda,
                                  Bwork+bshift2[jcur],iwrk,
                                  rone,
                                  C,ldc);
            /*
              mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                           alpha,
                           Aall+jcur*lda,lda,
                           Bwork+bshift2[jcur],iwrk,
                           rone,
                           C,ldc);
              */

              tmphml[4*taskid_j+2] =  Bwork[bshift2[jcur]];
          }
       }
   }

   if (jeven)
   {
      jcur = (jcur-1+np_j) % np_j;
      iwrk = na[jcur];
      CMatrix_end_rot(parall,request1_indx);
      if ((iwrk>0) && (mc[taskid_i]>0) &&  (nc[taskid_j]>0))
      {
          
         mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                              alpha,
                              work1,lda,
                              Bwork+bshift2[jcur],iwrk,
                              rone,
                              C,ldc);
         /*
          mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                           alpha,
                           Aall+jcur*lda,lda,
                           Bwork+bshift2[jcur],iwrk,
                           rone,
                           C,ldc);
          */

          tmphml[4*taskid_j+3] =  Bwork[bshift2[jcur]];
      }
   }
   else
   {
      jcur = (jcur-1+np_j) % np_j;
      iwrk = na[jcur];
      CMatrix_end_rot(parall,request1_indx);
      if ((iwrk>0) && (mc[taskid_i]>0) &&  (nc[taskid_j]>0))
      {
    
         mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                              alpha,
                              work2,lda,
                              Bwork+bshift2[jcur],iwrk,
                              rone,
                              C,ldc);
        /*
         mygdevice->NN_dgemm1(mc[taskid_i],nc[taskid_j],iwrk,
                           alpha,
                           Aall+jcur*lda,lda,
                           Bwork+bshift2[jcur],iwrk,
                           rone,
                           C,ldc);
         */

      }
   }
/*
   if (parall->is_master())
      std::cout << std::endl << "B= [" << Efmt(15,10)
                            << B[0] << " " << B[4] << " " << B[8]  << " " << B[12]  << std::endl
                << "      " << B[1] << " " << B[5] << " " << B[9]  << " " << B[13]  << std::endl
                << "      " << B[2] << " " << B[6] << " " << B[10] << " " << B[14]  << std::endl
                << "      " << B[3] << " " << B[7] << " " << B[11] << " " << B[15]  << "]"<< std::endl << std::endl;

   parall->Vector_SumAll(2,16,tmphml);
   if (parall->is_master())
      std::cout << "TMPHML= [" << Efmt(15,10)
                            << tmphml[0] << " " << tmphml[4] << " " << tmphml[8]  << " " << tmphml[12]  << std::endl
                << "      " << tmphml[1] << " " << tmphml[5] << " " << tmphml[9]  << " " << tmphml[13]  << std::endl
                << "      " << tmphml[2] << " " << tmphml[6] << " " << tmphml[10] << " " << tmphml[14]  << std::endl
                << "      " << tmphml[3] << " " << tmphml[7] << " " << tmphml[11] << " " << tmphml[15]  << "]"<< std::endl << std::endl;
                */

   //std::memcpy(C,C2,ldc*sizeof(double));
   //delete [] Aall;
   //delete [] C2;
}


} // namespace pwdft

