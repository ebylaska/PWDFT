// NWPW_CUDA Routines


#pragma once

/* can place sycl mkl code here */
#include        "blas.h"

#include	<cuda.h>
#include	<cuda_runtime_api.h>
#include	<cublas_v2.h>
#include	<cufft.h>

/* can place sycl mkl code here */
#include        <cstdio>
#include        <iostream>
#include        <limits>
#include        <CL/sycl.hpp>
#include        <oneapi/mkl.hpp>
#include        <sstream>
#include        <map>
#include        <set>


class Gdevices {
    cufftHandle forward_plan_x, forward_plan_y, forward_plan_z = 0;
    cufftHandle backward_plan_x, backward_plan_y, backward_plan_z = 0;

    cublasOperation_t matT = CUBLAS_OP_T;
    cublasOperation_t matN = CUBLAS_OP_N;


    /* device, host pool memory */
    std::map<size_t, std::set<double*> > free_list_gpu, free_list_host;
    std::map<double *, size_t> live_ptrs_gpu, live_ptrs_host;

    desc_real_t *desc_x;
    desc_cmplx_t *desc_y, *desc_z;

    cl::sycl::event fftevent;

    int nx_fft,ny_fft,nz_fft;

public:
    bool hasgpu = true;

    std::vector<cudaStream_t*> cudaQueues;  // TODO: queues per device
    cudaStream_t* device_queue = nullptr;   // default CUDA queue for now

    /* device memory */
    int    ndev_mem = 0;
    bool   inuse[25];
    size_t ndsize_mem[25];
    double *dev_mem[25];
    int    ia_psi,ia_hpsi;
    int    ib_prj;

    /* tmp memory */
    int    ntmp_mem=0;
    bool   tmpinuse[25];
    size_t tmpndsize_mem[25];
    double *tmp_mem[25];

    /* constructor */
    Gdevices() {
      std::cout << "calling gdevices constructor" << std::endl;
      ndev_mem = 0;

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
      device_queue =  new cl::sycl::queue(cl::sycl::gpu_selector{},
                                      asyncHandler,
                                      cl::sycl::property_list{cl::sycl::property::queue::in_order{}});
    }

    /* deconstructor */
    ~Gdevices() {
       std::cout << "calling gdevices destructor" << std::endl;

       // free fft descriptors
       cufftDestroy(forward_plan_x);
       cufftDestroy(forward_plan_y);
       cufftDestroy(forward_plan_z);
       cufftDestroy(backward_plan_x);
       cufftDestroy(backward_plan_y);
       cufftDestroy(backward_plan_z);


       // free CPU memory
       for (auto i=0; i<ntmp_mem; ++i)
          free(tmp_mem[i]);

       // free GPU memory
       for(std::map<size_t,std::set<double*>>::iterator it=free_list_gpu.begin(); it!=free_list_gpu.end(); ++it) {
          for(std::set<double*>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2) {
	    cudaFree(*it2);
          }
       }
       free_list_gpu.clear();
       live_ptrs_gpu.clear();


       // free host memory
       for(std::map<size_t,std::set<double*>>::iterator it=free_list_host.begin(); it!=free_list_host.end(); ++it) {
          for(std::set<double*>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2) {
             free(*it2);
          }
       }
       free_list_host.clear();
       live_ptrs_host.clear();



       delete device_queue;
    }



    static inline double *resurrect_from_free_list(std::map<size_t, std::set<double*> > &free_map,
                                                   size_t bytes,
                                                   std::map<double*, size_t>& liveset) {
        double* ptr=nullptr;
        assert(free_map.find(bytes) != free_map.end());
        /* assert(free_map.find(bytes)->second.size() > 0); */
        std::set<double*> &st = free_map.find(bytes)->second;
        ptr = *st.begin();
        st.erase(ptr);
        if(st.size()==0)
            free_map.erase(bytes);
        liveset[ptr] = bytes;
        return ptr;
    }

    double* getGpuMem(size_t bytes) {
        double *ptr=nullptr;
        if(free_list_gpu.find(bytes) != free_list_gpu.end()) {
            std::set<double*> &lst = free_list_gpu.find(bytes)->second;
            if(lst.size()!=0) {
                ptr = resurrect_from_free_list(free_list_gpu, bytes, live_ptrs_gpu);
                return ptr;
            }
        }
        else {
            for(std::map<size_t, std::set<double *> >::iterator it=free_list_gpu.begin();
                it != free_list_gpu.end();
                ++it)
            {
                if(it->first >= bytes && it->second.size()>0) {
                    ptr = resurrect_from_free_list(free_list_gpu, it->first, live_ptrs_gpu);
                    return ptr;
                }
            }
        }
    }

    double* getHostMem(size_t bytes) {
        double *ptr=nullptr;
        if(free_list_host.find(bytes)!=free_list_host.end()) {
            std::set<double*> &lst = free_list_host.find(bytes)->second;
            if(lst.size()!=0) {
                ptr = resurrect_from_free_list(free_list_host, bytes, live_ptrs_host);
                return ptr;
            }
        }
        else
        {
            for(std::map<size_t, std::set<double *> >::iterator it=free_list_host.begin();
                it != free_list_host.end();
                ++it)
            {
                if(it->first >= bytes && it->second.size()>0) {
                    ptr = resurrect_from_free_list(free_list_host, it->first, live_ptrs_host);
                    return ptr;
                }
            }
        }

        ptr = (double *)malloc(bytes);
        assert(ptr!=nullptr); /*We hopefully have a pointer*/
        live_ptrs_host[ptr] = bytes;
        return ptr;
    }

    void freeGpuMem(double *p) {
        assert(live_ptrs_gpu.find(p) != live_ptrs_gpu.end());
        size_t bytes = live_ptrs_gpu[p];
        device_queue->memset(p, 0, bytes);
        live_ptrs_gpu.erase(p);
        free_list_gpu[bytes].insert(p);
    }
    void freeHostMem(double *p) {
        assert(live_ptrs_host.find(p) != live_ptrs_host.end());
        size_t bytes = live_ptrs_host[p];
        memset(p, 0, bytes);
        live_ptrs_host.erase(p);
        free_list_host[bytes].insert(p);
    }

    int fetch_dev_mem_indx(const size_t ndsize) {
        int ii = 0;
        while ((((ndsize!=ndsize_mem[ii]) || inuse[ii])) && (ii<ndev_mem))
            ++ii;

        if (ii<ndev_mem) {
            inuse[ii] = true;
        } else {
            ii            = ndev_mem;
            inuse[ii]     = true;
            ndsize_mem[ii] = ndsize;
            cudaMalloc((void**)&(dev_mem[ii]), ndsize);
            ndev_mem += 1;
        }

        return ii;
    }

    int fetch_tmp_mem_indx(const size_t tmpndsize) {
       int ii = 0;
       while (((tmpndsize!=tmpndsize_mem[ii]) || tmpinuse[ii] ) && (ii<ntmp_mem))
         ++ii;

       if (ii<ntmp_mem) {
          tmpinuse[ii] = true;
       } else {
          ii                = ntmp_mem;
          tmpinuse[ii]      = true;
          tmpndsize_mem[ii] = tmpndsize;
          tmp_mem[ii]       = (double *) malloc(tmpndsize*sizeof(double));
          ntmp_mem += 1;
       }
       return ii;
    }



    void TN3_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cbb)
     {
#if 0
        int one = 1;
        int shift1  = 0;
        int mshift1 = 0;

        for (auto k=1; k<=ne; ++k)
        {
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_a,npack,&host_a[shift1],npack,beta,&host_caa[mshift1],k);
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_a,npack,&host_b[shift1],npack,beta,&host_cab[mshift1],k);
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_b,npack,&host_b[shift1],npack,beta,&host_cbb[mshift1],k);
           shift1  += npack;
           mshift1 += ne;
        }
#else
        int ic11 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
        int ic12 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
        int ic22 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

	cudaMemset(dev_mem[ic11],0,ne*ne*sizeof(double));
	cudaMemset(dev_mem[ic12],0,ne*ne*sizeof(double));
	cudaMemset(dev_mem[ic22],0,ne*ne*sizeof(double));

	cudaMemcpy(dev_mem[ia_psi], host_a,npack*ne*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_mem[ia_hpsi],host_b,npack*ne*sizeof(double),cudaMemcpyHostToDevice);

	cublasDgemm(master_handle,
		    matT, matN,
		    ne,ne,npack,&alpha,
		    dev_mem[ia_psi], npack,
		    dev_mem[ia_psi], npack,
		    &beta,dev_mem[ic11],ne);
	cublasDgemm(master_handle,
		    matT, matN,
		    ne,ne,npack,&alpha,
		    dev_mem[ia_psi], npack,
		    dev_mem[ia_hpsi], npack,
		    &beta,dev_mem[ic12],ne);
	cublasDgemm(master_handle,
		    matT, matN,
		    ne,ne,npack,&alpha,
		    dev_mem[ia_hpsi], npack,
		    dev_mem[ia_hpsi], npack,
		    &beta,dev_mem[ic22],ne);

	cudaMemcpy(host_caa,dev_mem[ic11],ne*ne*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(host_cab,dev_mem[ic12],ne*ne*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(host_cbb,dev_mem[ic22],ne*ne*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ic11] = false;
        inuse[ic12] = false;
        inuse[ic22] = false;
#endif
     }

     void NN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
#if 0
        DGEMM_PWDFT((char *) "N",(char *) "N",npack,ne,ne,alpha,host_a,npack,host_b,ne,beta,host_c,npack);
#else
        int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

	cudaMemcpy(dev_mem[ib],host_b,ne*ne*sizeof(double),cudaMemcpyHostToDevice);

	cublasDgemm(master_handle,
		    matN,matN,
		    npack,ne,ne,&alpha,
		    dev_mem[ia_psi],npack,
		    dev_mem[ib],ne,
		    beta,dev_mem[ia_hpsi],npack);

	cudaMemcpy(host_c,dev_mem[ia_hpsi],npack*ne*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ib] = false;

#endif

     }


     void TN_dgemm(int ne, int nprj, int npack, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        //gdevice_TN_dgemm(nn,nprj,ng,rtwo,a,b,rzero,sum);
#if 0
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,nprj,npack,alpha,host_a,npack,host_b,npack,beta,host_c,ne);
        
#else
        //int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        //int ib = fetch_dev_mem_indx(((size_t) npack) * ((size_t) nprj));
        ib_prj = fetch_dev_mem_indx(((size_t) npack) * ((size_t) nprj));
        int ic = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) nprj));

	//cudaMemcpy(dev_mem[ia],host_a, npack*ne*sizeof(double));
	//cudaMemcpy(dev_mem[ib],host_b,npack*nprj*sizeof(double));
	cudaMemcpy(dev_mem[ib_prj],host_b,npack*nprj*sizeof(double),cudaMemcpyHostToDevice);

	cublasDgemm(master_handle,
		    matT,matN,
		    ne,nprj,npack,&alpha,
		    dev_mem[ia_psi],npack,
		    dev_mem[ib_prj],npack,
		    &beta,dev_mem[ic],ne);

	cudaMemcpy(host_c,dev_mem[ic],ne*nprj*sizeof(double),cudaMemcpyDeviceToHost);

	   //inuse[ia] = false;
        //inuse[ib_prj] = false;
        inuse[ic] = false;

#endif
     }

     void NT_dgemm(int npack, int ne, int nprj, double alpha, double *host_a, double *host_b, double beta, double *host_c) {

#if 0
        DGEMM_PWDFT((char *) "N",(char *) "T",npack,ne,nprj,alpha,host_a,npack,host_b,ne,beta,host_c,npack);
#else
        int one = 1;
        int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) nprj));

	cudaMemcpy(dev_mem[ib],host_b,ne*nprj*sizeof(double),cudaMemcpyHostToDevice);

	cublasDgemm(master_handle,
		    matN,matT,
		    npack,ne,nprj,&alpha,
		    dev_mem[ib_prj],npack,
		    dev_mem[ib],ne,
		    &beta,dev_mem[ia_hpsi],npack);

        inuse[ib] = false;
        inuse[ib_prj] = false;
#endif
     }

     /* psi_dev functions*/
     void psi_alloc(int npack, int ne) {
        ia_psi  = fetch_dev_mem_indx(2*((size_t) npack) * ((size_t) ne));
        ia_hpsi = fetch_dev_mem_indx(2*((size_t) npack) * ((size_t) ne));
     }

     void psi_dealloc() {
        inuse[ia_psi]  = false;
        inuse[ia_hpsi] = false;
     }

     void psi_copy_host2gpu(int npack, int ne, double *psi) {
       cudaMemcpy(dev_mem[ia_psi],psi,2*npack*ne*sizeof(double),cudaMemcpyHostToDevice);
     }
     void hpsi_copy_host2gpu(int npack, int ne, double *hpsi) {
       cudaMemcpy(dev_mem[ia_hpsi],hpsi,2*npack*ne*sizeof(double),cudaMemcpyHostToDevice);
     }

     void psi_copy_gpu2host(int npack, int ne, double *psi) {
        cudaMemcpy(psi, dev_mem[ia_psi], 2*ne*npack*sizeof(double),cudaMemcpyDeviceToHost);
     }
     void hpsi_copy_gpu2host(int npack, int ne, double *hpsi) {
        cudaMemcpy(hpsi, dev_mem[ia_hpsi], 2*ne*npack*sizeof(double),cudaMemcpyDeviceToHost);
        device_queue->wait();
     }


     /* fft functions*/
     void batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) {
        cufftPlan1d(&forward_plan_x, nx, CUFFT_D2Z, nq1);
        cufftPlan1d(&forward_plan_y, ny, CUFFT_Z2Z, nq2);
        cufftPlan1d(&forward_plan_z, nz, CUFFT_Z2Z, nq3);

        cufftPlan1d(&backward_plan_x, nx, CUFFT_Z2D, nq1);
        cufftPlan1d(&backward_plan_y, ny, CUFFT_Z2Z, nq2);
        cufftPlan1d(&backward_plan_z, nz, CUFFT_Z2Z, nq3)

     }

     void batch_cfftx(bool forward, int nx, int nq, int n2ft3d, double *a) {

        int indx  = 0;
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));

	cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

	cufftResult result;
	if (forward) {
	  result = cufftExecD2Z(forward_plan_x,
				reinterpret_cast<cufftDoubleReal*>(dev_mem[ia_dev]),
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]));
	}
	else {
	  result = cufftExecZ2D(backward_plan_x,
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				reinterpret_cast<cufftDoubleReal*>(dev_mem[ia_dev]));
	}

	if (result != CUFFT_SUCCESS) {
	  std::stringstream msg;
	  msg << "cuFFT_x Error at " << __FILE__ << " : " << __LINE__ << ", "
	      << cudaPeekAtLastError() << std::endl;
	}

	cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ia_dev] = false;
     }


     void batch_cffty(bool forward, int ny,int nq,int n2ft3d, double *a) {
        int indx  = 0;
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));

	cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

	cufftResult result;
	if (forward) {
	  result = cufftExecZ2Z(forward_plan_y,
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				CUFFT_FORWARD);
	}
	else {
	  result = cufftExecZ2Z(backward_plan_y,
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				CUFFT_INVERSE);
	}

	if (result != CUFFT_SUCCESS) {
	  std::stringstream msg;
	  msg << "cuFFT_y Error at " << __FILE__ << " : " << __LINE__ << ", "
	      << cudaPeekAtLastError() << std::endl;
	}

	cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

	inuse[ia_dev] = false;
     }

     void batch_cfftz(bool forward, int nz,int nq,int n2ft3d, double *a) {
        int indx  = 0;
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));

	cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

	cufftResult result;
	if (forward) {
	  result = cufftExecZ2Z(forward_plan_z,
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				CUFFT_FORWARD);
	}
	else {
	  result = cufftExecZ2Z(backward_plan_z,
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
				CUFFT_INVERSE);
	}

	if (result != CUFFT_SUCCESS) {
	  std::stringstream msg;
	  msg << "cuFFT_z Error at " << __FILE__ << " : " << __LINE__ << ", "
	      << cudaPeekAtLastError() << std::endl;
	}

	cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ia_dev] = false;

     }

};




