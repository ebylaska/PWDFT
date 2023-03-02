#pragma once

// NWPW_SYCL Routines

// #include        <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <sycl/sycl.hpp>

#include "blas.h"
#include <oneapi/mkl.hpp>

typedef oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE,
                                     oneapi::mkl::dft::domain::REAL>
    desc_real_t;
typedef oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE,
                                     oneapi::mkl::dft::domain::COMPLEX>
    desc_cmplx_t;

namespace pwdft {

#define NWPW_SYCL_ERROR(EXPR)                                                  \
  try {                                                                        \
    EXPR;                                                                      \
  } catch (sycl::exception const &ex) {                                        \
    std::stringstream msg;                                                     \
    msg << "SYCL Exception at " << __FILE__ << " : " << __LINE__ << std::endl; \
    throw(std::runtime_error(ex.what()));                                      \
  } catch (std::runtime_error const &ex) {                                     \
    std::stringstream msg;                                                     \
    msg << "runtime Exception at " << __FILE__ << " : " << __LINE__            \
        << std::endl;                                                          \
    throw(std::runtime_error(ex.what()));                                      \
  }

#define NWPW_ONEMKL_ERROR(FN)                                                  \
  try {                                                                        \
    FN;                                                                        \
  } catch (oneapi::mkl::exception const &ex) {                                 \
    std::stringstream msg;                                                     \
    msg << "SYCL oneMKL exception at " << __FILE__ << " : " << __LINE__        \
        << std::endl;                                                          \
    throw(std::runtime_error(ex.what()));                                      \
    exit(1);                                                                   \
  }

  namespace detail {
    enum memcpy_direction {
      host_to_host,
      host_to_device,
      device_to_host,
      device_to_device,
      automatic
    };


    // Get actual copy range and make sure it will not exceed range.
    static inline size_t get_copy_range(sycl::range<3> size, size_t slice,
					size_t pitch) {
      return slice * (size.get(2) - 1) + pitch * (size.get(1) - 1) + size.get(0);
    }

    static inline size_t get_offset(sycl::id<3> id, size_t slice,
                                    size_t pitch) {
      return slice * id.get(2) + pitch * id.get(1) + id.get(0);
    }

    static inline std::vector<sycl::event>
    dpct_memcpy(sycl::queue &q, void *to_ptr, const void *from_ptr,
                sycl::range<3> to_range, sycl::range<3> from_range,
                sycl::id<3> to_id, sycl::id<3> from_id,
                sycl::range<3> size, memcpy_direction direction,
                const std::vector<sycl::event> &dep_events = {}) {
      // RAII for host pointer
      class host_buffer {
        void *_buf;
        size_t _size;
        sycl::queue &_q;
        const std::vector<sycl::event> &_deps; // free operation depends

      public:
        host_buffer(size_t size, sycl::queue &q,
                    const std::vector<sycl::event> &deps)
          : _buf(std::malloc(size)), _size(size), _q(q), _deps(deps) {}
        void *get_ptr() const { return _buf; }
        size_t get_size() const { return _size; }
        ~host_buffer() {
          if (_buf) {
            _q.submit([&](sycl::handler &cgh) {
                cgh.depends_on(_deps);
                cgh.host_task([buf = _buf] { std::free(buf); });
              });
          }
        }
      };
      std::vector<sycl::event> event_list;

      size_t to_slice = to_range.get(1) * to_range.get(0),
        from_slice = from_range.get(1) * from_range.get(0);
      unsigned char *to_surface = (unsigned char *)to_ptr + detail::get_offset(to_id, to_slice, to_range.get(0));
      const unsigned char *from_surface = (const unsigned char *)from_ptr +
        detail::get_offset(from_id, from_slice, from_range.get(0));

      if (to_slice == from_slice && to_slice == size.get(1) * size.get(0)) {
	return {q.memcpy(to_surface, from_surface, to_slice * size.get(2), dep_events)};
      }

      //direction = deduce_memcpy_direction(q, to_ptr, from_ptr, direction);
      size_t size_slice = size.get(1) * size.get(0);
      switch (direction) {
      case host_to_device: {
        host_buffer buf(detail::get_copy_range(size, to_slice, to_range.get(0)), q, event_list);
        std::vector<sycl::event> host_events;
        if (to_slice == size_slice) {
          // Copy host data to a temp host buffer with the shape of target.
          host_events = detail::dpct_memcpy(q, buf.get_ptr(), from_surface, to_range, from_range,
                        sycl::id<3>(0, 0, 0), sycl::id<3>(0, 0, 0), size,
                        host_to_host, dep_events);
        } else {
          // Copy host data to a temp host buffer with the shape of target.
          host_events = detail::dpct_memcpy(
                                    q, buf.get_ptr(), from_surface, to_range, from_range,
                                    sycl::id<3>(0, 0, 0), sycl::id<3>(0, 0, 0), size, host_to_host,
                                    // If has padding data, not sure whether it is useless. So fill temp
                                    // buffer with it.
                                    std::vector<sycl::event>{q.memcpy(buf.get_ptr(), to_surface, buf.get_size(),
								      dep_events)});

        }
        // Copy from temp host buffer to device with only one submit.
        event_list.push_back(q.memcpy(to_surface, buf.get_ptr(), buf.get_size(), host_events));
        break;
      }
      case device_to_host: {
        host_buffer buf(detail::get_copy_range(size, from_slice, from_range.get(0)), q,
                        event_list);
        // Copy from host temp buffer to host target with reshaping.
        event_list = detail::dpct_memcpy(q, to_surface, buf.get_ptr(), to_range, from_range, sycl::id<3>(0, 0, 0),
                                 sycl::id<3>(0, 0, 0), size, host_to_host,
                                 // Copy from device to temp host buffer with only one submit.
                                 std::vector<sycl::event>{q.memcpy(buf.get_ptr(), from_surface,
								   buf.get_size(), dep_events)});
        break;
      }
      default:
        throw std::runtime_error("dpct_memcpy: invalid direction value");
      }
      return event_list;
    }

    /// memcpy 2D matrix with pitch.
    static inline std::vector<sycl::event>
    dpct_memcpy(sycl::queue &q, void *to_ptr, const void *from_ptr,
                size_t to_pitch, size_t from_pitch, size_t x, size_t y,
                memcpy_direction direction = automatic) {
      return dpct_memcpy(q, to_ptr, from_ptr, sycl::range<3>(to_pitch, y, 1),
                         sycl::range<3>(from_pitch, y, 1),
                         sycl::id<3>(0, 0, 0), sycl::id<3>(0, 0, 0),
                         sycl::range<3>(x, y, 1), direction);
    }
  }

  inline void syclSetMatrixAsync( int rows, int cols, size_t elem_size,
                                  const void *from_ptr, int from_ld,
                                  void *to_ptr, int to_ld,
                                  sycl::queue* que) {
    if (to_ptr == from_ptr && to_ld == from_ld) {
      return;
    }

    if (to_ld == from_ld) {
      size_t cpoy_size = elem_size * ((cols - 1) * to_ld + rows);
      que->memcpy( to_ptr, from_ptr, cpoy_size );
    } else {
      detail::dpct_memcpy(*que, to_ptr, from_ptr, elem_size * to_ld,
                          elem_size * from_ld, elem_size * rows, cols,
                          detail::host_to_device);
    }
  }

  inline void syclGetMatrixAsync( int rows, int cols, size_t elem_size,
                                  const void *from_ptr, int from_ld,
                                  void *to_ptr, int to_ld,
                                  sycl::queue* que) {
    if (to_ptr == from_ptr && to_ld == from_ld) {
      return;
    }

    if (to_ld == from_ld) {
      size_t cpoy_size = elem_size * ((cols - 1) * to_ld + rows);
      que->memcpy( to_ptr, from_ptr, cpoy_size );
    } else {
      detail::dpct_memcpy(*que, to_ptr, from_ptr, elem_size * to_ld,
                          elem_size * from_ld, elem_size * rows, cols,
                          detail::device_to_host);
    }
  }

sycl::queue *get_syclQue();

class Gdevices {

  oneapi::mkl::transpose matT = oneapi::mkl::transpose::trans;
  oneapi::mkl::transpose matN = oneapi::mkl::transpose::nontrans;

  /* device, host pool memory */
  std::map<size_t, std::set<double *>> free_list_gpu, free_list_host;
  std::map<double *, size_t> live_ptrs_gpu, live_ptrs_host;

  desc_real_t *desc_x;
  desc_cmplx_t *desc_y, *desc_z;

  sycl::event h2d_event, fftevent, d2h_event;

public:
  bool hasgpu = true;

  std::vector<sycl::queue *> stream;

  /* device memory */
  int    ndev_mem = 0;
  bool   inuse[19] = {false};
  size_t ndsize_mem[19];
  double *dev_mem[19];
  int    tile_fac=1;
  int    tile_npack2_max;
  int    tile_npack2[19],tile_start2[19];
  double *a_psi,*a_hpsi,*b_prj;
  int    ia_psi[2],ia_hpsi[2],ib_prj[2];

  /* constructor */
  Gdevices() {
    std::cout << "calling gdevices constructor" << std::endl;
    ndev_mem = 0;

    auto asyncHandler = [&](sycl::exception_list eL) {
      for (auto &e : eL) {
        try {
          std::rethrow_exception(e);
        } catch (sycl::exception &e) {
          std::cout << e.what() << std::endl;
          std::cout << "fail" << std::endl;
          std::terminate();
        }
      }
    };

    // allocate SYCL streams
    for (auto i=0; i<2; ++i) {
      stream.push_back(new sycl::queue(sycl::gpu_selector_v, asyncHandler,
				       sycl::property_list{sycl::property::queue::in_order{}}));
    }
  }

  /* deconstructor */
  ~Gdevices() {

    // free dev_mem
    for (auto i=0; i<ndev_mem; ++i) sycl::free(dev_mem[i], *stream[0]);
    ndev_mem = 0;

    // free cuda streams
    for (auto i=0; i<2; ++i) delete stream[i];
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
      dev_mem[ii] = sycl::malloc_device<double>(ndsize, *stream[0]);
      ndev_mem += 1;
    }

    stream[0]->memset(dev_mem[ii], 0, ndsize*sizeof(double)).wait();
    return ii;
  }

   /**************************************
    *                                    *
    *              TN4_dgemm             *
    *                                    *
    **************************************/
    /* This function computes <host_a|host_a>, <host_a|host_b>, <host_b|host_a>, and <host_b|host_b> overlap matrices.

        host_caa = beta*host_caa + alpha*host_a'*host_a
        host_cab = beta*host_cab + alpha*host_a'*host_b
        host_cba = beta*host_cba + alpha*host_b'*host_a
        host_cbb = beta*host_cbb + alpha*host_b'*host_b

       Entry - npack2,ne: matrix size
               alpha, beta: standard dgemm parameters
               host_a: (npack2xne) matrix
               host_b: (npack2xne) matrix
       Exit - host_caa,host_cab,host_cba, host_cbb: (nexne) matrices
       Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated previously with psi_alloc
            - temporary device memory for (nexne) matrices ic11, ic12, and ic22.
    */
  void TN4_dgemm(int npack2, int ne, double alpha, double *host_a,
                 double *host_b, double beta, double *host_caa,
                 double *host_cab, double *host_cba, double *host_cbb) {
    int ic11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
    int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
    int ic21 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
    int ic22 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));

    if (std::fabs(beta)>0.0) {
      stream[0]->memcpy(dev_mem[ic11], host_caa, ne*ne*sizeof(double));
      stream[0]->memcpy(dev_mem[ic12], host_cab, ne*ne*sizeof(double));
      stream[0]->memcpy(dev_mem[ic21], host_cba, ne*ne*sizeof(double));
      stream[0]->memcpy(dev_mem[ic22], host_cbb, ne*ne*sizeof(double));
    }

    // copy host_a,host_b --> dev_mem
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       &host_a[tile_start2[0]],npack2,
                       dev_mem[ia_psi[0]],tile_npack2[0],stream[0]);
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       &host_b[tile_start2[0]],npack2,
                       dev_mem[ia_hpsi[0]],tile_npack2[0],stream[0]);

    double beta0 = beta;
    for (auto tt=0; tt<tile_fac; ++tt) {
      int ttp1 = tt+1;
      if (ttp1<tile_fac) {
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &host_a[tile_start2[ttp1]],npack2,
                           dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &host_b[tile_start2[ttp1]],npack2,
                           dev_mem[ia_hpsi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
      }
      stream[tt%2]->wait();;

      oneapi::mkl::blas::column_major::gemm(*stream[tt%2],
                  matT, matN,
                  ne,ne,tile_npack2[tt],alpha,
                  dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                  dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                  beta0,dev_mem[ic11],ne);
      oneapi::mkl::blas::column_major::gemm(*stream[tt%2],
                  matT, matN,
                  ne,ne,tile_npack2[tt],alpha,
                  dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                  dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                  beta0,dev_mem[ic12],ne);
      oneapi::mkl::blas::column_major::gemm(*stream[tt%2],
                  matT, matN,
                  ne,ne,tile_npack2[tt],alpha,
                  dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                  dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                  beta0,dev_mem[ic21],ne);
      oneapi::mkl::blas::column_major::gemm(*stream[tt%2],
                  matT, matN,
                  ne,ne,tile_npack2[tt],alpha,
                  dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                  dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                  beta0,dev_mem[ic22],ne);
      beta0 = 1.0;
    }


    stream[0]->memcpy(host_caa, dev_mem[ic11], ne * ne * sizeof(double));
    stream[0]->memcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double));
    stream[0]->memcpy(host_cba, dev_mem[ic21], ne * ne * sizeof(double));
    stream[0]->memcpy(host_cbb, dev_mem[ic22], ne * ne * sizeof(double));

    stream[0]->wait();

    inuse[ic11] = false;
    inuse[ic12] = false;
    inuse[ic21] = false;
    inuse[ic22] = false;
  }


   /**************************************
    *                                    *
    *              TN3_dgemm             *
    *                                    *
    **************************************/
    /* This function computes <host_a|host_a>, <host_a|host_b>, and <host_b|host_b> overlap matrices.

        host_caa = beta*host_caa + alpha*host_a'*host_a
        host_cab = beta*host_cab + alpha*host_a'*host_b
        host_cbb = beta*host_cbb + alpha*host_b'*host_b

       Entry - npack2,ne: matrix size
               alpha, beta: standard dgemm parameters
               host_a: (npack2xne) matrix
               host_b: (npack2xne) matrix
       Exit - host_caa,host_cab,host_cbb: (nexne) matrices
       Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated previously with psi_alloc
            - temporary device memory for (nexne) matrices ic11, ic12, and ic22.
    */
  void TN3_dgemm(int npack2, int ne, double alpha, double *host_a,
                 double *host_b, double beta, double *host_caa,
                 double *host_cab, double *host_cbb) {
    int ic11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
    int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
    int ic22 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));

    if (std::fabs(beta)>0.0) {
      stream[0]->memcpy(dev_mem[ic11], host_caa, ne*ne*sizeof(double));
      stream[0]->memcpy(dev_mem[ic12], host_cab, ne*ne*sizeof(double));
      stream[0]->memcpy(dev_mem[ic22], host_cbb, ne*ne*sizeof(double));
    }

    // copy host_a,host_b --> dev_mem
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       &host_a[tile_start2[0]],npack2,
                       dev_mem[ia_psi[0]],tile_npack2[0],stream[0]);
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       &host_b[tile_start2[0]],npack2,
                       dev_mem[ia_hpsi[0]],tile_npack2[0],stream[0]);

    double beta0 = beta;
    for (auto tt=0; tt<tile_fac; ++tt) {
      int ttp1 = tt+1;
      if (ttp1<tile_fac) {
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &host_a[tile_start2[ttp1]],npack2,
                           dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &host_b[tile_start2[ttp1]],npack2,
                           dev_mem[ia_hpsi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
      }
      stream[tt%2]->wait();;

      oneapi::mkl::blas::column_major::gemm(*stream[tt%2],
                  matT, matN,
                  ne,ne,tile_npack2[tt],alpha,
                  dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                  dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                  beta0,dev_mem[ic11],ne);
      oneapi::mkl::blas::column_major::gemm(*stream[tt%2],
                  matT, matN,
                  ne,ne,tile_npack2[tt],alpha,
                  dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                  dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                  beta0,dev_mem[ic12],ne);
      oneapi::mkl::blas::column_major::gemm(*stream[tt%2],
                  matT, matN,
                  ne,ne,tile_npack2[tt],alpha,
                  dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                  dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                  beta0,dev_mem[ic22],ne);
      beta0 = 1.0;
    }


    stream[0]->memcpy(host_caa, dev_mem[ic11], ne * ne * sizeof(double));
    stream[0]->memcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double));
    stream[0]->memcpy(host_cbb, dev_mem[ic22], ne * ne * sizeof(double));

    stream[0]->wait();

    inuse[ic11] = false;
    inuse[ic12] = false;
    inuse[ic22] = false;
  }

  /**************************************
   *                                    *
   *              TN1_dgemm             *
   *                                    *
   **************************************/
  /* This function computes  <host_a|host_b> overlap matrix.

     host_cab = beta*host_cab + alpha*host_a'*host_b

     Entry - npack2,ne: matrix size
     alpha, beta: standard dgemm parameters
     host_a:  (npack2xne) matrix
     host_b:  (npack2xne) matrix
     Exit - host_cab: (nexne) matrices
     Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated previously with psi_alloc
     - temporary device memory for (nexne) matrix, ic12.
  */
  void TN1_dgemm(int npack2, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_cab) {
    int ic12 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

    if (std::fabs(beta)>0.0) {
      stream[0]->memcpy(dev_mem[ic12], host_cab,ne*ne*sizeof(double)).wait();
    }

    // copy host_a,host_b --> dev_mem
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       &host_a[tile_start2[0]],npack2,
                       dev_mem[ia_psi[0]],tile_npack2[0],stream[0]);
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       &host_b[tile_start2[0]],npack2,
                       dev_mem[ia_hpsi[0]],tile_npack2[0],stream[0]);

    double beta0 = beta;
    for (auto tt=0; tt<tile_fac; ++tt) {
      int ttp1 = tt+1;
      if (ttp1<tile_fac) {
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &host_a[tile_start2[ttp1]],npack2,
                           dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &host_b[tile_start2[ttp1]],npack2,
                           dev_mem[ia_hpsi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
      }
      stream[tt%2]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[0],
                matT, matN,
                ne,ne,tile_npack2[tt],alpha,
                dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                beta0,dev_mem[ic12],ne);
      beta0 = 1.0;
    }

    stream[0]->memcpy(host_cab,dev_mem[ic12],ne*ne*sizeof(double)).wait();

    inuse[ic12] = false;
  }

  /**************************************
   *                                    *
   *              NN_dgemm              *
   *                                    *
   **************************************/
  void NN_dgemm(int npack2, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
    //DGEMM_PWDFT((char *) "N",(char *) "N",npack2,ne,ne,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);

    int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

    syclSetMatrixAsync(ne,ne,sizeof(double),host_b,ne,dev_mem[ib],ne,stream[0]);
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       &host_a[tile_start2[0]],npack2,
                       dev_mem[ia_psi[0]],tile_npack2[0],stream[0]);
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       &host_c[tile_start2[0]],npack2,
                       dev_mem[ia_hpsi[0]],tile_npack2[0],stream[0]);

    //double beta0 = beta;
    for (auto tt=0; tt<tile_fac; ++tt) {
      int ttp1 = tt+1;
      if (ttp1<tile_fac) {
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &host_a[tile_start2[ttp1]],npack2,
                           dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &host_c[tile_start2[ttp1]],npack2,
                           dev_mem[ia_hpsi[ttp1%2]], tile_npack2[ttp1],stream[ttp1%2]);
      }
      stream[tt%2]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[0],
                                            matN,matN,
                                            tile_npack2[tt],ne,ne,alpha,
                                            dev_mem[ia_psi[tt%2]],tile_npack2[tt],
                                            dev_mem[ib],ne,
                                            beta,dev_mem[ia_hpsi[tt%2]],tile_npack2[tt]);
      syclGetMatrixAsync(tile_npack2[tt],ne,sizeof(double),
                         dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                         &host_c[tile_start2[tt]],npack2,stream[tt%2]);

    }

    stream[(tile_fac-1)%2]->wait();

    inuse[ib] = false;
  }


  /**************************************
   *                                    *
   *              TN_dgemm              *
   *                                    *
   **************************************/
  void TN_dgemm(int ne, int nprj, int npack2, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
    //DGEMM_PWDFT((char *) "T",(char *) "N",ne,nprj,npack2,alpha,host_a,npack2,host_b,npack2,beta,host_c,ne);


    //gdevice_TN_dgemm(nn,nprj,ng,rtwo,a,b,rzero,sum);

    //int ia = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) ne));
    //int ib = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) nprj));
    b_prj  = host_b;
    ib_prj[0] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) nprj));
    if (tile_fac>1) ib_prj[1] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) nprj));
    int ic = fetch_dev_mem_indx(((size_t) ne)         * ((size_t) nprj));

    syclSetMatrixAsync(ne,nprj,sizeof(double),host_c,ne,dev_mem[ic],ne,stream[0]);

    if (tile_fac>1)
      syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                         &a_psi[tile_start2[0]],npack2,
                         dev_mem[ia_psi[0]],tile_npack2[0],stream[0]);
    syclSetMatrixAsync(tile_npack2[0],nprj,sizeof(double),
                       &b_prj[tile_start2[0]],npack2,
                       dev_mem[ib_prj[0]],tile_npack2[0],stream[0]);

    double beta0 = beta;
    for (auto tt=0; tt<tile_fac; ++tt) {
      int ttp1 = tt+1;
      if (ttp1<tile_fac) {
        syclSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                           &a_psi[tile_start2[ttp1]],npack2,
                           dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
        syclSetMatrixAsync(tile_npack2[ttp1],nprj,sizeof(double),
                           &b_prj[tile_start2[ttp1]],npack2,
                           dev_mem[ib_prj[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]);
      }

      stream[tt%2]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[0],
                                            matT,matN,
                                            ne,nprj,tile_npack2[tt],alpha,
                                            dev_mem[ia_psi[tt%2]],tile_npack2[tt],
                                            dev_mem[ib_prj[tt%2]],tile_npack2[tt],
                                            beta0,dev_mem[ic],ne);
      beta0 = 1.0;
    }
    stream[0]->memcpy(host_c,dev_mem[ic],ne*nprj*sizeof(double)).wait();

    //inuse[ia] = false;
    //inuse[ib_prj[0]] = false;
    //if (tile_fac>1) inuse[ib_prj[1]] = false;
    inuse[ic] = false;

  }

  void T_free() { inuse[ib_prj[0]] = false; if (tile_fac>1) inuse[ib_prj[1]] = false; }

  void NT_dgemm(int npack2, int ne, int nprj, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
    //DGEMM_PWDFT((char *) "N",(char *) "T",npack2,ne,nprj,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);

    int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) nprj));

    syclSetMatrixAsync(ne,nprj,sizeof(double),
                       host_b,ne,
                       dev_mem[ib],ne,stream[(tile_fac-1)%2]);
    syclSetMatrixAsync(tile_npack2[tile_fac-1],ne,sizeof(double),
                       &host_c[tile_start2[tile_fac-1]],npack2,
                       dev_mem[ia_hpsi[(tile_fac-1)%2]],tile_npack2[tile_fac-1],stream[(tile_fac-1)%2]);
    syclSetMatrixAsync(tile_npack2[tile_fac-1],nprj,sizeof(double),
                       &host_a[tile_start2[tile_fac-1]],npack2,
                       dev_mem[ib_prj[(tile_fac-1)%2]],tile_npack2[tile_fac-1],stream[(tile_fac-1)%2]);
    for (auto tt=tile_fac-1; tt>=0; --tt) {
      int ttm1 = tt-1;
      if (ttm1>=0) {
        syclSetMatrixAsync(tile_npack2[ttm1],ne,sizeof(double),
                           &host_c[tile_start2[ttm1]],npack2,
                           dev_mem[ia_hpsi[ttm1%2]],tile_npack2[ttm1],stream[ttm1%2]);
        syclSetMatrixAsync(tile_npack2[ttm1],nprj,sizeof(double),
                           &host_a[tile_start2[ttm1]],npack2,
                           dev_mem[ib_prj[ttm1%2]], tile_npack2[ttm1],stream[ttm1%2]);
      }
      stream[tt%2]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[0],
                matN,matT,
                tile_npack2[tt],ne,nprj, alpha,
                dev_mem[ib_prj[tt%2]],tile_npack2[tt],
                dev_mem[ib],ne,
                beta,dev_mem[ia_hpsi[tt%2]],tile_npack2[tt]);
      syclGetMatrixAsync(tile_npack2[tt],ne,sizeof(double),
                         dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                         &host_c[tile_start2[tt]],npack2,stream[tt%2]);
    }

    inuse[ib] = false;
    inuse[ib_prj[0]] = false;
    if (tile_fac>1) inuse[ib_prj[1]] = false;
  }


  /**************************************
   *                                    *
   *              MM6_dgemm             *
   *                                    *
   **************************************/
   void MM6_dgemm(int ne,
                   double *host_s21, double *host_s12, double *host_s11,
                   double *host_sa0, double *host_sa1, double *host_st1) {
      double rzero=0.0;
      double rone =1.0;
      int i_s21 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input
      int i_s12 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input
      int i_s11 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input
      int i_sa0 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input
      int i_st1 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //tmp
      int i_sa1 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input-output

      syclSetMatrixAsync(ne,ne,sizeof(double),host_s21,ne,dev_mem[i_s21],ne,stream[0]);
      syclSetMatrixAsync(ne,ne,sizeof(double),host_sa0,ne,dev_mem[i_sa0],ne,stream[0]);
      syclSetMatrixAsync(ne,ne,sizeof(double),host_sa1,ne,dev_mem[i_sa1],ne,stream[0]);

      syclSetMatrixAsync(ne,ne,sizeof(double),host_s12,ne,dev_mem[i_s12],ne,stream[1]);
      syclSetMatrixAsync(ne,ne,sizeof(double),host_s11,ne,dev_mem[i_s11],ne,stream[1]);

      //mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
      stream[0]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[0],
                                            matN,matN,
                                            ne,ne,ne,rone,
                                            dev_mem[i_s21],ne,
                                            dev_mem[i_sa0],ne,
                                            rone,dev_mem[i_sa1],ne);
 
      //mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
      stream[1]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[1],
                                            matN,matN,
                                            ne,ne,ne,rone,
                                            dev_mem[i_sa0],ne,
                                            dev_mem[i_s12],ne,
                                            rone,dev_mem[i_sa1],ne);

      //mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
      oneapi::mkl::blas::column_major::gemm(*stream[1],
                                            matN,matN,
                                            ne,ne,ne,rone,
                                            dev_mem[i_s11],ne,
                                            dev_mem[i_sa0],ne,
                                            rzero,dev_mem[i_st1],ne);

      //mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
      oneapi::mkl::blas::column_major::gemm(*stream[1],
                                            matN,matN,
                                            ne,ne,ne,rone,
                                            dev_mem[i_sa0],ne,
                                            dev_mem[i_st1],ne,
                                            rone,dev_mem[i_sa1],ne);
      syclGetMatrixAsync(ne,ne,sizeof(double),
                         dev_mem[i_sa1],ne,
                         host_sa1,ne,stream[1]);
      stream[1]->wait();

      inuse[i_s21] = false;
      inuse[i_s12] = false;
      inuse[i_s11] = false;
      inuse[i_sa0] = false;
      inuse[i_st1] = false;
      inuse[i_sa1] = false;
   }



  /********************/
  /* psi_dev functions*/
  /********************/
  void psi_alloc(int npack1, int ne, int tile_fac0=1) {
    tile_fac        = tile_fac0;

    tile_npack2_max = (((2*npack1)%tile_fac)==0) ? (2*npack1)/tile_fac : (2*npack1)/tile_fac + 1;
    //for (auto i=0; i<tile_fac; ++i) tile_npack2[i] = (i<((2*npack1)%tile_fac)) ? (2*npack1)/tile_fac+1 : (2*npack1)/tile_fac;
    for (auto i=0; i<tile_fac; ++i) tile_npack2[i] = (2*npack1)/tile_fac;
    for (auto i=0; i<((2*npack1)%tile_fac); ++i) tile_npack2[i] += 1;

    tile_start2[0] = 0;
    for (auto i=1; i<tile_fac; ++i) tile_start2[i] = tile_start2[i-1] + tile_npack2[i-1];

    ia_psi[0]  = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
    ia_hpsi[0] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));

    if (tile_fac>1) {
      ia_psi[1]  = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
      ia_hpsi[1] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
    }
    std::cout << "Into psi_alloc, tile_factor = " << tile_fac << " ndev_mem=" << ndev_mem << std::endl;
  }
  void psi_dealloc() {
    inuse[ia_psi[0]]  = false;
    inuse[ia_hpsi[0]] = false;
    if (tile_fac>1) {
      inuse[ia_psi[1]]  = false;
      inuse[ia_hpsi[1]] = false;
    }
  }
  void psi_copy_host2gpu(int npack1, int ne, double *psi) {
    a_psi = psi;
    syclSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                       psi,2*npack1,
                       dev_mem[ia_psi[0]],tile_npack2[0],stream[0]);
  }
  void hpsi_copy_host2gpu(int npack1, int ne, double *hpsi) {
    int tt = tile_fac-1;
    a_hpsi = hpsi;
    syclSetMatrixAsync(tile_npack2[tt],ne,sizeof(double),
                       &hpsi[tile_start2[tt]],2*npack1,
                       dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],stream[tt%2]);
  }
  void psi_copy_gpu2host(int npack1, int ne, double *psi) {
    if (tile_fac==1) {
      stream[0]->wait();
      syclGetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                         dev_mem[ia_psi[0]],tile_npack2[0],
                         psi,2*npack1,stream[0]);
    }
  }
  void hpsi_copy_gpu2host(int npack1, int ne, double *hpsi) {
    stream[0]->wait();
  }

  /* fft functions*/
  void batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) {
    desc_x = new desc_real_t(nx);
    desc_x->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                      nq1);
    desc_x->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, nx + 2);
    desc_x->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, nx + 2);

    desc_y = new desc_cmplx_t(ny);
    desc_y->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                      nq2);
    desc_y->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, ny);
    desc_y->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, ny);

    desc_z = new desc_cmplx_t(nz);
    desc_z->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                      nq3);
    desc_z->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, nz);
    desc_z->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, nz);

    desc_x->commit(*stream[0]);
    desc_y->commit(*stream[0]);
    desc_z->commit(*stream[0]);
  }

  void batch_fft_end() {
    delete desc_x;
    delete desc_y;
    delete desc_z;

    ndev_mem = 0;
  }

  void batch_cfftx(bool forward, int nx, int nq, int n2ft3d, double *a) {
    int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));

    stream[0]->memcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double));

    if (forward)
      compute_forward(*desc_x, dev_mem[ia_dev]);
    else
      compute_backward(*desc_x, dev_mem[ia_dev]);

    stream[0]->memcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double));
    stream[0]->wait();

    inuse[ia_dev] = false;
  }

  void batch_cffty(bool forward, int ny, int nq, int n2ft3d, double *a) {
    int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));

    stream[0]->memcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double));

    if (forward)
      compute_forward(*desc_y, dev_mem[ia_dev]);
    else
      compute_backward(*desc_y, dev_mem[ia_dev]);

    stream[0]->memcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double));
    stream[0]->wait();

    inuse[ia_dev] = false;
  }

  void batch_cfftz(bool forward, int nz, int nq, int n2ft3d, double *a) {
    int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));

    stream[0]->memcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double));

    if (forward)
      compute_forward(*desc_z, dev_mem[ia_dev]);
    else
      compute_backward(*desc_z, dev_mem[ia_dev]);

    stream[0]->memcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double));
    stream[0]->wait();

    inuse[ia_dev] = false;
  }


  // routines below need to be made into sycl or removed

static void eigsrt_device(double *D, double *V, int n) {
   int i,j,k;
   double p;

   for (i=0; i<(n-1); ++i)
   {
      k = i;
      p = D[i];
      for(j=i+1; j<n; ++j)
         if (D[j]>=p)
         {
            k = j;
            p = D[j];
         }

      if (k!=i)
      {
         D[k] = D[i];
         D[i] = p;
         for (j=0; j<n; ++j)
         {
            p = V[j+i*n];
            V[j+i*n] = V[j+k*n];
            V[j+k*n] = p;
         }
      }
   }
}

  void NN_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) {
       int n,ierr;
       int nn  = ne[0]*ne[0]+14;
       double xmp1[nn];
       //double *xmp1 = new (std::nothrow) double[nn]();

       int shift1 = 0;
       int shift2 = 0;
       for (int ms=0; ms<ispin; ++ms)
       {
          n = ne[ms];

          //eigen_(&n,&n,&hml[shift2],&eig[shift1],xmp1,&ierr);
          // d3db::parall->Barrier();
          EIGEN_PWDFT(n,host_hml+shift2,host_eig+shift1,xmp1,nn,ierr);
          //if (ierr != 0) throw std::runtime_error(std::string("NWPW Error: EIGEN_PWDFT failed!"));

          eigsrt_device(host_eig+shift1,host_hml+shift2,n);
          shift1 += ne[0];
          shift2 += ne[0]*ne[0];
       }
  }


}; // class Gdevices

} // namespace pwdft
