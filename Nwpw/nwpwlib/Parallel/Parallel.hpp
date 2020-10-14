#ifndef _Parallel_H_
#define _Parallel_H_
/* Parallel.h
   Author - Eric Bylaska
	this class is used defining nd parallel geometries
*/

#include	"mpi.h"

#ifdef HAS_SYCL
#include	<cstdio>
#include	<iostream>
#include	<limits>
#include	<CL/sycl.hpp>
#include	"mkl_blas_sycl.hpp"
#endif

#define	MASTER	0

class Parallel {

    int npi[3],taskidi[3];
    int *reqcnt;
    int *procNd;
    //MPI::Intracomm comm_i[3];
    //MPI::Group     group_i[3];
    //MPI::Request  **request;
    MPI_Comm      comm_world;
    MPI_Comm      comm_i[3];
    MPI_Group     group_i[3];
    MPI_Request  **request;
    MPI_Status   **statuses;

#ifdef HAS_SYCL
    auto asyncHandler = [&](cl::sycl::exception_list eL) {
       for (auto& e : eL) {
         try {
           std::rethrow_exception(e);
         } catch (cl::sycl::exception& e) {
           std::cout << e.what() << std::endl;
           std::cout << "fail" << std::endl;
           std::terminate();
         }
       }
    };
    cl::sycl::gpu_selector device_selector;
    cl::sycl::queue device_queue(device_selector,asyncHandler);
#endif

public:
        int dim;

	/* Constructors */
	//Parallel(int, char **);
	Parallel(MPI_Comm);

        /* destructor */
	~Parallel();

        /* 2d proc geom constructor */
        void init2d(const int, const int);

        int is_master() { return (taskidi[0]==MASTER); }

        int taskid()   {return taskidi[0];}
        int taskid_i() {return taskidi[1];}
        int taskid_j() {return taskidi[2];}
        int np()   {return npi[0];}
        int np_i() {return npi[1];}
        int np_j() {return npi[2];}

        int convert_taskid_i(const int i) {return procNd[i+taskidi[2]*npi[1]]; }
        int convert_taskid_j(const int j) {return procNd[taskidi[1]+j*npi[1]]; }
        int convert_taskid_ij(const int i, const int j) {return procNd[i+j*npi[1]]; }

       /* SumAll */
      double SumAll(const int, const double);
      void Vector_SumAll(const int, const int, double *);
      int ISumAll(const int, const int);
      void Vector_ISumAll(const int, const int, int *);

       /* MaxAll */
      double MaxAll(const int, const double);

       /* Brdcsts */
      void Brdcst_Values(const int, const int, const int, double *);
      void Brdcst_iValues(const int, const int, const int, int *);
      void Brdcst_iValue(const int, const int, int *);
      void Brdcst_cValues(const int, const int, const int, void *);

      /* send/receives */
      void    dsend(const int, const int, const int, const int, double *);
      void dreceive(const int, const int, const int, const int, double *);
      void    isend(const int, const int, const int, const int, int    *);
      void ireceive(const int, const int, const int, const int, int    *);

      /* asend/areceives */
      void astart(const int, const int);
      void aend(const int);
      void    adsend(const int, const int, const int, const int, double *);
      void adreceive(const int, const int, const int, const int, double *);

};

#endif
