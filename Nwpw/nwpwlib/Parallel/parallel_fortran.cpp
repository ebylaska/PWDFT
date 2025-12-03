/* parallel_fortran.cpp
   Author - Eric Bylaska
        this class is used for defining 3d parallel maps
*/

//#include <iostream> // DEBUG


#include "mpi.h"

static MPI_Comm comm_world99 = MPI_COMM_WORLD;


/**************************************
 *                                    *
 *      parallel_fortran_init         *
 *                                    *
 **************************************/
extern "C" void parallel_fortran_init(MPI_Comm comm_world0)
{
  comm_world99 = comm_world0;
}


/**************************************
 *                                    *
 *      parallel_fortran_sumall       *
 *                                    *
 **************************************/
extern "C" void parallel_fortran_sumall_(double *sum) 
{
   double sumout;
   MPI_Allreduce(sum, &sumout, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world99);
   *sum = sumout;
}


/************************************
 *                                  *
 *    parallel_fortran_isumall_     *
 *                                  *
 ************************************/
extern "C" void parallel_fortran_isumall_(int *sum) 
{
   int sumout;
   MPI_Allreduce(sum, &sumout, 1, MPI_INTEGER, MPI_SUM, comm_world99);
   *sum = sumout;
}


/****************************************
 *                                      *
 *    parallel_fortran_vector_sumall_   *
 *                                      *
 ****************************************/
extern "C" void parallel_fortran_vector_sumall_(int *n, double *sum) 
{
   double *sumout = new double[*n];

   MPI_Allreduce(sum, sumout, *n, MPI_DOUBLE_PRECISION, MPI_SUM, comm_world99);
   for (int i = 0; i<(*n); ++i)
      sum[i] = sumout[i];
   delete[] sumout;
}


/******************************************
 *                                        *
 *     parallel_fortran_vector_isumall_   *
 *                                        *
 ******************************************/
extern "C" void parallel_fortran_vector_isumall_(int *n, int *sum) 
{
   int *sumout = new int[*n];
   MPI_Allreduce(sum, sumout, *n, MPI_INTEGER, MPI_SUM, comm_world99);
   for (int i = 0; i < (*n); ++i)
      sum[i] = sumout[i];
   delete[] sumout;
}

/******************************************
 *                                        *
 *     parallel_fortran_taskid_           *
 *                                        *
 ******************************************/
extern "C" void parallel_fortran_taskid_(int *myid)
{
    MPI_Comm_rank(comm_world99, myid);
}


/******************************************
 *                                        *
 *     parallel_fortran_numprocs_         *
 *                                        *
 ******************************************/
extern "C" void parallel_fortran_numprocs_(int *nprocs)
{
    MPI_Comm_size(comm_world99, nprocs);
}


