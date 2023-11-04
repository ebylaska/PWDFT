/* parallel_fortran.cpp
   Author - Eric Bylaska
        this class is used for defining 3d parallel maps
*/

#include "mpi.h"

#define MASTER 0

static MPI_Comm comm_world99;


/**************************************
 *                                    *
 *      parallel_fortran_init         *
 *                                    *
 **************************************/
void parallel_fortran_init(MPI_Comm comm_world0)
{
  comm_world99 = comm_world0;
}


/**************************************
 *                                    *
 *      parallel_fortran_SumAll       *
 *                                    *
 **************************************/
extern "C" void parallel_fortran_SumAll_(double *sum) 
{
   double sumout;
   MPI_Allreduce(sum, &sumout, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world99);
   *sum = sumout;
}


/************************************
 *                                  *
 *    parallel_fortran_SumAll       *
 *                                  *
 ************************************/
extern "C" void parallel_fortran_ISumAll_(int *sum) 
{
   int sumout;
   MPI_Allreduce(sum, &sumout, 1, MPI_INTEGER, MPI_SUM, comm_world99);
   *sum = sumout;
}


/****************************************
 *                                      *
 *    parallel_fortran_Vector_SumAll_   *
 *                                      *
 ****************************************/
extern "C" void fortran_Vector_SumAll_(int *n, double *sum) 
{
   double *sumout = new double[*n];

   MPI_Allreduce(sum, sumout, *n, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world99);
   for (int i = 0; i<(*n); ++i)
      sum[i] = sumout[i];
   delete[] sumout;
}


/******************************************
 *                                        *
 *     parallel_fortran_Vector_ISumAll    *
 *                                        *
 ******************************************/
extern "C" void fortran_Vector_ISumAll_(const int d, const int n, int *sum) 
{
   int *sumout = new int[n];
   MPI_Allreduce(sum, sumout, n, MPI_INTEGER, MPI_SUM, comm_world99);
   for (int i = 0; i < n; ++i)
      sum[i] = sumout[i];
   delete[] sumout;
}

