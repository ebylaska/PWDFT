#ifndef _Parallel_HPP_
#define _Parallel_HPP_

/* Parallel.h
   Author - Eric Bylaska
        this class is used defining nd parallel geometries
*/

#include "mpi.h"

namespace pwdft {

#define MASTER 0

/**
 * @class Parallel
 * @brief Provides MPI-based parallelism and communication management.
 *
 * The `Parallel` class encapsulates functionality for managing MPI-based parallel
 * computations. It allows users to distribute tasks across multiple MPI processes,
 * perform collective operations, and handle data communication between processes.
 *
 * This class is designed to work in a multi-dimensional parallel computing environment,
 * supporting both 1D and 2D process grids.
 */

class Parallel {
   int thrid = 0;
   int nthr = 1;
   int max_nthr = 1;
 
   int npi[4], taskidi[4];
   int *reqcnt;
   int *procNd;
   // MPI::Intracomm comm_i[3];
   // MPI::Group     group_i[3];
   // MPI::Request  **request;
   MPI_Comm comm_world;
   MPI_Comm comm_i[4];
   MPI_Group group_i[4];
   MPI_Request **request;
   MPI_Status **statuses;
   bool init2d_called = false;
   bool init3d_called = false;

public:
   int max_reqstat;
   int dim;
   bool base_stdio_print = true;
 
   /* Constructors */
   // Parallel(int, char **);
   Parallel(MPI_Comm);
 
   /* destructor */
   ~Parallel();
 
   /* 2d proc geom constructor */
   void init2d(const int, const int);

   /* 3d proc geom constructor */
   void init3d(const int, const int, const int);
 
   int is_master() { return (taskidi[0] == MASTER); }
   int is_master_d(const int d) { return (taskidi[d] == MASTER); }
 
   int threadid() { return thrid; }
   int nthreads() { return nthr; }
   int maxthreads() { return max_nthr; }
 
   int taskid() { return taskidi[0]; }
   int taskid_i() { return taskidi[1]; }
   int taskid_j() { return taskidi[2]; }
   int taskid_k() { return taskidi[3]; }
   int np() { return npi[0]; }
   int np_i() { return npi[1]; }
   int np_j() { return npi[2]; }
   int np_k() { return npi[3]; }
 
   int convert_taskid_i(const int i) { return procNd[i + taskidi[2] * npi[1]]; }
   int convert_taskid_j(const int j) { return procNd[taskidi[1] + j * npi[1]]; }
   int convert_taskid_ij(const int i, const int j) {
     return procNd[i + j * npi[1]];
   }

   int convert_taskid_ijk(const int i, const int j, const int k) {
     return procNd[i + j * npi[1] + k *npi[1]*npi[2]];
   }
 
   /* Barriers */
   void Barrier() { MPI_Barrier(comm_world); }
   void comm_Barrier(const int i) { MPI_Barrier(comm_i[i]); }
 
   /* SumAll */
   double SumAll(const int, const double);
   void Vector_SumAll(const int, const int, double *);
   void Vector_SumAll_buffer(const int, const int, double *, double *);
   int ISumAll(const int, const int);
   void Vector_ISumAll(const int, const int, int *);
 
   /* MaxAll */
   double MaxAll(const int, const double);
 
   /* Brdcsts */
   void Brdcst_Values(const int, const int, const int, double *);
   void Brdcst_iValues(const int, const int, const int, int *);
   void Brdcst_iValue(const int, const int, int *);
   void Brdcst_cValues(const int, const int, const int, void *);
 
   /* Reduce */
   void Reduce_Values(const int, const int, const int, double *, double *);
 
   /* send/receives */
   void dsend(const int, const int, const int, const int, double *);
   void dreceive(const int, const int, const int, const int, double *);
   void isend(const int, const int, const int, const int, int *);
   void ireceive(const int, const int, const int, const int, int *);
 
   /* asend/areceives */
   void astart(const int, const int);
   void awaitall(const int);
   void aend(const int);
   void adsend(const int, const int, const int, const int, double *);
   void adreceive(const int, const int, const int, const int, double *);

   void a2dsend(const int, const int, const int, const int, double *);
   void a2dreceive(const int, const int, const int, const int, double *);
};

} // namespace pwdft

#endif
