// NWPW_CUDA Routines

#pragma once

#include        "blas.h"

#include        "gdevice.hpp"

class Gdevices {

    cufftHandle forward_plan_x = 0 , plan_y = 0, plan_z = 0;
    cufftHandle backward_plan_x = 0;

    cublasHandle_t master_handle = 0;
    cublasOperation_t matT = CUBLAS_OP_T;
    cublasOperation_t matN = CUBLAS_OP_N;

public:
    bool hasgpu = true;

    /* device memory */
    int    ndev_mem = 0;
    bool   inuse[25] = {false};
    size_t ndsize_mem[25];
    double *dev_mem[25];
    int    ia_psi,ia_hpsi;
    int    ib_prj;

    /* constructor */
    Gdevices() {
        ndev_mem = 0;

        NWPW_CUBLAS_ERROR( cublasCreate(&master_handle) );
    }

    /* deconstructor */
    ~Gdevices() {
        cublasDestroy(master_handle);

        // free fft descriptors
        cufftDestroy(forward_plan_x);
        cufftDestroy(plan_y);
        cufftDestroy(plan_z);
        cufftDestroy(backward_plan_x);
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
            NWPW_CUDA_ERROR( cudaMalloc((void**)&(dev_mem[ii]), ndsize*sizeof(double)) );
            ndev_mem += 1;
        }

	NWPW_CUDA_ERROR( cudaMemset(dev_mem[ii], 0, ndsize*sizeof(double)) );
        return ii;
    }

    void TN3_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cbb)
        {
            int ic11 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
            int ic12 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
            int ic22 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

            NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ia_psi], host_a,npack*ne*sizeof(double),cudaMemcpyHostToDevice) );
            NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ia_hpsi],host_b,npack*ne*sizeof(double),cudaMemcpyHostToDevice) );

            NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                           matT, matN,
                                           ne,ne,npack,&alpha,
                                           dev_mem[ia_psi], npack,
                                           dev_mem[ia_psi], npack,
                                           &beta,dev_mem[ic11],ne) );
            NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                           matT, matN,
                                           ne,ne,npack,&alpha,
                                           dev_mem[ia_psi], npack,
                                           dev_mem[ia_hpsi], npack,
                                           &beta,dev_mem[ic12],ne) );
            NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                           matT, matN,
                                           ne,ne,npack,&alpha,
                                           dev_mem[ia_hpsi], npack,
                                           dev_mem[ia_hpsi], npack,
                                           &beta,dev_mem[ic22],ne) );

            NWPW_CUDA_ERROR( cudaMemcpy(host_caa,dev_mem[ic11],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );
            NWPW_CUDA_ERROR( cudaMemcpy(host_cab,dev_mem[ic12],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );
            NWPW_CUDA_ERROR( cudaMemcpy(host_cbb,dev_mem[ic22],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );

            inuse[ic11] = false;
            inuse[ic12] = false;
            inuse[ic22] = false;
        }

    void NN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

        NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ib],host_b,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );

        NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
				       matN,matN,
				       npack,ne,ne,&alpha,
				       dev_mem[ia_psi],npack,
				       dev_mem[ib],ne,
				       &beta,dev_mem[ia_hpsi],npack) );

        NWPW_CUDA_ERROR( cudaMemcpy(host_c,dev_mem[ia_hpsi],npack*ne*sizeof(double),cudaMemcpyDeviceToHost) );

        inuse[ib] = false;
    }


    void TN_dgemm(int ne, int nprj, int npack, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        //gdevice_TN_dgemm(nn,nprj,ng,rtwo,a,b,rzero,sum);

        //int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        //int ib = fetch_dev_mem_indx(((size_t) npack) * ((size_t) nprj));
        ib_prj = fetch_dev_mem_indx(((size_t) npack) * ((size_t) nprj));
        int ic = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) nprj));

        //cudaMemcpy(dev_mem[ia],host_a, npack*ne*sizeof(double));
        //cudaMemcpy(dev_mem[ib],host_b,npack*nprj*sizeof(double));
        cudaMemcpy(dev_mem[ib_prj],host_b,npack*nprj*sizeof(double),cudaMemcpyHostToDevice);

        NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
				       matT,matN,
				       ne,nprj,npack,&alpha,
				       dev_mem[ia_psi],npack,
				       dev_mem[ib_prj],npack,
				       &beta,dev_mem[ic],ne) );

        cudaMemcpy(host_c,dev_mem[ic],ne*nprj*sizeof(double),cudaMemcpyDeviceToHost);

        //inuse[ia] = false;
        //inuse[ib_prj] = false;
        inuse[ic] = false;
    }

    void NT_dgemm(int npack, int ne, int nprj, double alpha, double *host_a, double *host_b, double beta, double *host_c) {

        int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) nprj));
        NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ib],host_b,ne*nprj*sizeof(double),cudaMemcpyHostToDevice) );

        NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
				       matN,matT,
				       npack,ne,nprj,&alpha,
				       dev_mem[ib_prj],npack,
				       dev_mem[ib],ne,
				       &beta,dev_mem[ia_hpsi],npack) );

        inuse[ib] = false;
        inuse[ib_prj] = false;

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
    }


    /* fft functions (uses cuFFT) */
    void batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) {
        NWPW_CUFFT_ERROR( cufftPlan1d(&forward_plan_x, nx, CUFFT_D2Z, nq1) );
        NWPW_CUFFT_ERROR( cufftPlan1d(&backward_plan_x, nx, CUFFT_Z2D, nq1) );

        int y_inembed[] = {ny};
        int y_onembed[] = {ny};
        NWPW_CUFFT_ERROR( cufftPlanMany(&plan_y, 1, &ny, y_inembed, 1, ny, y_onembed, 1, ny, CUFFT_Z2Z, nq2) );

        int z_inembed[] = {nz};
        int z_onembed[] = {nz};
        NWPW_CUFFT_ERROR( cufftPlanMany(&plan_z, 1, &nz, z_inembed, 1, nz, z_onembed, 1, nz, CUFFT_Z2Z, nq3) );
    }

    void batch_cfftx(bool forward, int nx, int nq, int n2ft3d, double *a) {
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
        cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

        if (forward) {
            NWPW_CUFFT_ERROR( cufftExecD2Z(forward_plan_x,
                                           reinterpret_cast<cufftDoubleReal*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev])) );
        }
        else {
            NWPW_CUFFT_ERROR( cufftExecZ2D(backward_plan_x,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleReal*>(dev_mem[ia_dev])) );
        }

        cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ia_dev] = false;
    }

    void batch_cffty(bool forward, int ny,int nq,int n2ft3d, double *a) {
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
        cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

        if (forward) {
            NWPW_CUFFT_ERROR( cufftExecZ2Z(plan_y,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           CUFFT_FORWARD) );
        }
        else {
            NWPW_CUFFT_ERROR( cufftExecZ2Z(plan_y,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           CUFFT_INVERSE) );
        }

        cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ia_dev] = false;
    }

    void batch_cfftz(bool forward, int nz,int nq,int n2ft3d, double *a) {
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
        cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

        if (forward) {
            NWPW_CUFFT_ERROR( cufftExecZ2Z(plan_z,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           CUFFT_FORWARD) );
        }
        else {
            NWPW_CUFFT_ERROR( cufftExecZ2Z(plan_z,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           CUFFT_INVERSE) );
        }

        cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ia_dev] = false;
    }

};
