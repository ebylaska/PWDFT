/* d1db.C
   Author - Eric Bylaska

        this class is used for defining 3d parallel maps
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
#include "d1db.hpp"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

namespace pwdft {

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

d1db::d1db(Parallel *inparall, const int inmaptype, const int inispin, int *inne)
     : Mapping1(inmaptype, inparall->np_j(), inparall->taskid_j(), inispin, inne) 
{
  parall = inparall;
}


/***********************************
 *                                 *
 *      d1db::DMatrix_dgemm1       *
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

void d1db::DMatrix_dgemm1(Parallel *parall, 
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
 *         DMatrix_dgemm2          *
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

void d1db::DMatrix_dgemm2(Parallel *parall,
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
 *         DMatrix_dgemm3          *
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
void d1db::DMatrix_dgemm3(Parallel *parall,
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
 *         DMatrix_dgemm2c         *
 *                                 *
 ***********************************/

void d1db::DMatrix_dgemm2c(Parallel *parall,
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


} // namespace pwdft

