#pragma once

// NWPW_SYCL Routines

#define NDEV_MAX  39
#define DEBUG_IO  false


// #include        <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <sycl/sycl.hpp>
#include <complex>

#include "fft.h"
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

static inline size_t get_offset(sycl::id<3> id, size_t slice, size_t pitch) {
  return slice * id.get(2) + pitch * id.get(1) + id.get(0);
}

static inline std::vector<sycl::event>
dpct_memcpy(sycl::queue &q, void *to_ptr, const void *from_ptr,
            sycl::range<3> to_range, sycl::range<3> from_range,
            sycl::id<3> to_id, sycl::id<3> from_id, sycl::range<3> size,
            memcpy_direction direction,
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
  unsigned char *to_surface =
      (unsigned char *)to_ptr +
      detail::get_offset(to_id, to_slice, to_range.get(0));
  const unsigned char *from_surface =
      (const unsigned char *)from_ptr +
      detail::get_offset(from_id, from_slice, from_range.get(0));

  if (to_slice == from_slice && to_slice == size.get(1) * size.get(0)) {
    return {
        q.memcpy(to_surface, from_surface, to_slice * size.get(2), dep_events)};
  }

  // direction = deduce_memcpy_direction(q, to_ptr, from_ptr, direction);
  size_t size_slice = size.get(1) * size.get(0);
  switch (direction) {
  case host_to_device: {
     host_buffer buf(detail::get_copy_range(size, to_slice, to_range.get(0)), q,
                     event_list);
     std::vector<sycl::event> host_events;
     if (to_slice == size_slice) {
       // Copy host data to a temp host buffer with the shape of target.
       host_events = detail::dpct_memcpy(
           q, buf.get_ptr(), from_surface, to_range, from_range,
           sycl::id<3>(0, 0, 0), sycl::id<3>(0, 0, 0), size, host_to_host,
           dep_events);
     } else {
       // Copy host data to a temp host buffer with the shape of target.
       host_events = detail::dpct_memcpy(
           q, buf.get_ptr(), from_surface, to_range, from_range,
           sycl::id<3>(0, 0, 0), sycl::id<3>(0, 0, 0), size, host_to_host,
           // If has padding data, not sure whether it is useless. So fill temp
           // buffer with it.
           std::vector<sycl::event>{
               q.memcpy(buf.get_ptr(), to_surface, buf.get_size(), dep_events)});
     }
     // Copy from temp host buffer to device with only one submit.
     event_list.push_back(
         q.memcpy(to_surface, buf.get_ptr(), buf.get_size(), host_events));
     break;
  }
  case device_to_host: {
     host_buffer buf(detail::get_copy_range(size, from_slice, from_range.get(0)),
                     q, event_list);
     // Copy from host temp buffer to host target with reshaping.
     event_list = detail::dpct_memcpy(
         q, to_surface, buf.get_ptr(), to_range, from_range,
         sycl::id<3>(0, 0, 0), sycl::id<3>(0, 0, 0), size, host_to_host,
         // Copy from device to temp host buffer with only one submit.
         std::vector<sycl::event>{
             q.memcpy(buf.get_ptr(), from_surface, buf.get_size(), dep_events)});
     break;
  }
  default:
    throw std::runtime_error("dpct_memcpy: invalid direction value");
  }
  return event_list;
}

/// memcpy 2D matrix with pitch.
static inline std::vector<sycl::event>
dpct_memcpy(sycl::queue &q, void *to_ptr, const void *from_ptr, size_t to_pitch,
            size_t from_pitch, size_t x, size_t y,
            memcpy_direction direction = automatic) {
  return dpct_memcpy(q, to_ptr, from_ptr, sycl::range<3>(to_pitch, y, 1),
                     sycl::range<3>(from_pitch, y, 1), sycl::id<3>(0, 0, 0),
                     sycl::id<3>(0, 0, 0), sycl::range<3>(x, y, 1), direction);
}

} // namespace detail

inline void syclSetMatrixAsync(int rows, int cols, size_t elem_size,
                               const void *from_ptr, int from_ld, void *to_ptr,
                               int to_ld, sycl::queue *que) {
   if (to_ptr == from_ptr && to_ld == from_ld) {
     return;
   }
 
   if (to_ld == from_ld) {
     size_t cpoy_size = elem_size * ((cols - 1) * to_ld + rows);
     que->memcpy(to_ptr, from_ptr, cpoy_size);
   } else {
     detail::dpct_memcpy(*que, to_ptr, from_ptr, elem_size * to_ld,
                         elem_size * from_ld, elem_size * rows, cols,
                         detail::host_to_device);
   }
}

inline void syclGetMatrixAsync(int rows, int cols, size_t elem_size,
                               const void *from_ptr, int from_ld, void *to_ptr,
                               int to_ld, sycl::queue *que) {
   if (to_ptr == from_ptr && to_ld == from_ld) {
      return;
   }
 
   if (to_ld == from_ld) {
      size_t cpoy_size = elem_size * ((cols - 1) * to_ld + rows);
      que->memcpy(to_ptr, from_ptr, cpoy_size);
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
   oneapi::mkl::transpose matC = oneapi::mkl::transpose::conjtrans;
 
   /* device, host pool memory */
   std::map<size_t, std::set<double *>> free_list_gpu, free_list_host;
   std::map<double *, size_t> live_ptrs_gpu, live_ptrs_host;
 
   int fftcount = 0;
   desc_real_t *desc_r_x[2];
   desc_cmplx_t *desc_x[2], *desc_y[2], *desc_z[2];
 
   sycl::event h2d_event, fftevent, d2h_event;

public:
   int typegpu = 2;
   bool hasgpu = true;
 
   std::vector<sycl::queue *> stream;
 
   /* device memory */
   int ndev_mem = 0;
   bool inuse[NDEV_MAX] = {false};
   size_t ndsize_mem[NDEV_MAX];
   double *dev_mem[NDEV_MAX];
   int tile_fac = 1;
   int tile_npack2_max;
   int tile_npack2[19], tile_start2[19];
   int tile_npack1[19], tile_start1[19];
   double *a_psi, *a_hpsi, *b_prj;
   int ia_psi[2], ia_hpsi[2], ib_prj[2];
 
   int ifft_dev[15];
   int ifft_n;
 
   /* constructor */
   /**************************************
    *                                    *
    *             Gdevices               *
    *                                    *
    **************************************/
   Gdevices() 
   {
      if (DEBUG_IO) std::cout << "calling gdevices constructor" << std::endl;
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
     
      for (auto i=0; i<12; ++i) {
        stream.push_back(new sycl::queue(
            sycl::gpu_selector_v, asyncHandler,
            sycl::property_list{sycl::property::queue::in_order{}}));
      }
   }
 
   /* deconstructor */
   /**************************************
    *                                    *
    *            ~Gdevices               *
    *                                    *
    **************************************/
   ~Gdevices() 
   {
      // free dev_mem
      for (auto i=0; i<ndev_mem; ++i)
         sycl::free(dev_mem[i], *stream[0]);
      ndev_mem = 0;
    
      // free cuda streams
      for (auto i=0; i<12; ++i)
         delete stream[i];
   }
 
   /**************************************
    *                                    *
    *            fetch_dev_mem_indx      *
    *                                    *
    **************************************/
   int fetch_dev_mem_indx(const size_t ndsize) 
   {
      int ii = 0;
      while ((((ndsize != ndsize_mem[ii]) || inuse[ii])) && (ii < ndev_mem))
         ++ii;
     
      if (ii < ndev_mem) 
      {
         inuse[ii] = true;
      } 
      else 
      {
         ii = ndev_mem;
         inuse[ii] = true;
         ndsize_mem[ii] = ndsize;
         dev_mem[ii] = sycl::malloc_device<double>(ndsize, *stream[0]);
         ndev_mem += 1;
         if (ndev_mem>NDEV_MAX) std::cout << "ERROR: ndev_mem > NDEV_MAX" << std::endl;
      }
     
      stream[0]->memset(dev_mem[ii], 0, ndsize * sizeof(double)).wait();
      return ii;
   }
 
   /**************************************
    *                                    *
    *              TN4_dgemm             *
    *                                    *
    **************************************/
   /* This function computes <host_a|host_a>, <host_a|host_b>, <host_b|host_a>,
      and <host_b|host_b> overlap matrices.
 
       host_caa = beta*host_caa + alpha*host_a'*host_a
       host_cab = beta*host_cab + alpha*host_a'*host_b
       host_cba = beta*host_cba + alpha*host_b'*host_a
       host_cbb = beta*host_cbb + alpha*host_b'*host_b
 
      Entry - npack2,ne: matrix size
              alpha, beta: standard dgemm parameters
              host_a: (npack2xne) matrix
              host_b: (npack2xne) matrix
      Exit - host_caa,host_cab,host_cba, host_cbb: (nexne) matrices
      Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated
      previously with psi_alloc
           - temporary device memory for (nexne) matrices ic11, ic12, and ic22.
   */
   void TN4_dgemm(int npack2, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_caa,
                  double *host_cab, double *host_cba, double *host_cbb) 
   {
      int ic11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
      int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
      int ic21 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
      int ic22 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     
      if (std::fabs(beta) > 0.0) {
        stream[0]->memcpy(dev_mem[ic11], host_caa, ne * ne * sizeof(double));
        stream[0]->memcpy(dev_mem[ic12], host_cab, ne * ne * sizeof(double));
        stream[0]->memcpy(dev_mem[ic21], host_cba, ne * ne * sizeof(double));
        stream[0]->memcpy(dev_mem[ic22], host_cbb, ne * ne * sizeof(double));
      }
     
      // copy host_a,host_b --> dev_mem
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         &host_a[tile_start2[0]], npack2, dev_mem[ia_psi[0]],
                         tile_npack2[0], stream[0]);
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         &host_b[tile_start2[0]], npack2, dev_mem[ia_hpsi[0]],
                         tile_npack2[0], stream[0]);
     
      double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) 
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) 
         {
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &host_a[tile_start2[ttp1]], npack2,
                               dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &host_b[tile_start2[ttp1]], npack2,
                               dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
         }
         stream[tt % 2]->wait();
        
         oneapi::mkl::blas::column_major::gemm(
             *stream[tt % 2], matT, matN, ne, ne, tile_npack2[tt], alpha,
             dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_psi[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic11], ne);
         oneapi::mkl::blas::column_major::gemm(
             *stream[tt % 2], matT, matN, ne, ne, tile_npack2[tt], alpha,
             dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic12], ne);
         oneapi::mkl::blas::column_major::gemm(
             *stream[tt % 2], matT, matN, ne, ne, tile_npack2[tt], alpha,
             dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], dev_mem[ia_psi[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic21], ne);
         oneapi::mkl::blas::column_major::gemm(
             *stream[tt % 2], matT, matN, ne, ne, tile_npack2[tt], alpha,
             dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic22], ne);
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
   /* This function computes <host_a|host_a>, <host_a|host_b>, and
      <host_b|host_b> overlap matrices.
 
       host_caa = beta*host_caa + alpha*host_a'*host_a
       host_cab = beta*host_cab + alpha*host_a'*host_b
       host_cbb = beta*host_cbb + alpha*host_b'*host_b
 
      Entry - npack2,ne: matrix size
              alpha, beta: standard dgemm parameters
              host_a: (npack2xne) matrix
              host_b: (npack2xne) matrix
      Exit - host_caa,host_cab,host_cbb: (nexne) matrices
      Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated
      previously with psi_alloc
           - temporary device memory for (nexne) matrices ic11, ic12, and ic22.
   */
   void TN3_dgemm(int npack2, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_caa,
                  double *host_cab, double *host_cbb) 
   {
      int ic11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
      int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
      int ic22 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     
      if (std::fabs(beta) > 0.0) 
      {
         stream[0]->memcpy(dev_mem[ic11], host_caa, ne * ne * sizeof(double));
         stream[0]->memcpy(dev_mem[ic12], host_cab, ne * ne * sizeof(double));
         stream[0]->memcpy(dev_mem[ic22], host_cbb, ne * ne * sizeof(double));
      }
     
      // copy host_a,host_b --> dev_mem
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         &host_a[tile_start2[0]], npack2, dev_mem[ia_psi[0]],
                         tile_npack2[0], stream[0]);
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         &host_b[tile_start2[0]], npack2, dev_mem[ia_hpsi[0]],
                         tile_npack2[0], stream[0]);
     
      double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) 
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) {
           syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                              &host_a[tile_start2[ttp1]], npack2,
                              dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
                              stream[ttp1 % 2]);
           syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                              &host_b[tile_start2[ttp1]], npack2,
                              dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
                              stream[ttp1 % 2]);
         }
         stream[tt % 2]->wait();
        
         oneapi::mkl::blas::column_major::gemm(
             *stream[tt % 2], matT, matN, ne, ne, tile_npack2[tt], alpha,
             dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_psi[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic11], ne);
         oneapi::mkl::blas::column_major::gemm(
             *stream[tt % 2], matT, matN, ne, ne, tile_npack2[tt], alpha,
             dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic12], ne);
         oneapi::mkl::blas::column_major::gemm(
             *stream[tt % 2], matT, matN, ne, ne, tile_npack2[tt], alpha,
             dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic22], ne);
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
      Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated
      previously with psi_alloc
      - temporary device memory for (nexne) matrix, ic12.
   */
   void TN1_dgemm(int npack2, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_cab) 
   {
      int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     
      if (std::fabs(beta) > 0.0) 
      {
         stream[0]->memcpy(dev_mem[ic12], host_cab, ne*ne*sizeof(double)).wait();
      }
     
      // copy host_a,host_b --> dev_mem
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         &host_a[tile_start2[0]], npack2, dev_mem[ia_psi[0]],
                         tile_npack2[0], stream[0]);
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         &host_b[tile_start2[0]], npack2, dev_mem[ia_hpsi[0]],
                         tile_npack2[0], stream[0]);
     
      double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) 
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) 
         {
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &host_a[tile_start2[ttp1]], npack2,
                               dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &host_b[tile_start2[ttp1]], npack2,
                               dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
         }
         stream[tt % 2]->wait();
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matT, matN, ne, ne, tile_npack2[tt], alpha,
             dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic12], ne);
         beta0 = 1.0;
      }
     
      stream[0]->memcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double)).wait();
     
      inuse[ic12] = false;
   }
 
   /**************************************
    *                                    *
    *              NN_dgemm              *
    *                                    *
    **************************************/
   void NN_dgemm(int npack2, int ne, double alpha, double *host_a,
                 double *host_b, double beta, double *host_c) 
   {
      // DGEMM_PWDFT((char *) "N",(char *)
      // "N",npack2,ne,ne,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);
     
      int ib = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     
      syclSetMatrixAsync(ne, ne, sizeof(double), host_b, ne, dev_mem[ib], ne, stream[0]);
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         &host_a[tile_start2[0]], npack2, dev_mem[ia_psi[0]],
                         tile_npack2[0], stream[0]);
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         &host_c[tile_start2[0]], npack2, dev_mem[ia_hpsi[0]],
                         tile_npack2[0], stream[0]);
     
      // double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) 
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) 
         {
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &host_a[tile_start2[ttp1]], npack2,
                               dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &host_c[tile_start2[ttp1]], npack2,
                               dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
         }
         stream[tt % 2]->wait();
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matN, matN, tile_npack2[tt], ne, ne, alpha,
             dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ib], ne, beta,
             dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt]);
         syclGetMatrixAsync(tile_npack2[tt], ne, sizeof(double),
                            dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt],
                            &host_c[tile_start2[tt]], npack2, stream[tt % 2]);
      }
     
      stream[(tile_fac - 1) % 2]->wait();
     
      inuse[ib] = false;
   }
 
   /**************************************
    *                                    *
    *              NN_dgemm1             *
    *                                    *
    **************************************/
   void NN_dgemm1(int m, int n, int k,
                  double alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double beta,
                  double *host_c,int ldc) 
   {
      int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)k));
      int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)n));
      int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
 
      syclSetMatrixAsync(lda,k,sizeof(double),host_a,lda,dev_mem[ia],lda,stream[0]);
      syclSetMatrixAsync(ldb,n,sizeof(double),host_b,ldb,dev_mem[ib],ldb,stream[0]);
 
      stream[0]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[0], 
         matN,matN,m,n,k, 
         alpha,
         dev_mem[ia],lda, 
         dev_mem[ib],ldb, 
         beta,
         dev_mem[ic],ldc);
 
      syclGetMatrixAsync(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]);
      stream[0]->wait();
 
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }
 
   /**************************************
    *                                    *
    *              TN_dgemm2             *
    *                                    *
    **************************************/
   void TN_dgemm2(int m, int n, int k,
                  double alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double beta,
                  double *host_c,int ldc) 
   {
      int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)m));
      int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)n));
      int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
 
      syclSetMatrixAsync(lda,m,sizeof(double),host_a,lda,dev_mem[ia],lda,stream[0]);
      syclSetMatrixAsync(ldb,n,sizeof(double),host_b,ldb,dev_mem[ib],ldb,stream[0]);
 
      stream[0]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[0], 
         matT,matN,m,n,k, 
         alpha,
         dev_mem[ia],lda, 
         dev_mem[ib],ldb, 
         beta,
         dev_mem[ic],ldc);
 
      syclGetMatrixAsync(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]);
      stream[0]->wait();
 
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }
 
   /**************************************
    *                                    *
    *              TN_dgemm2c            *
    *                                    *
    **************************************/
/**
  * @brief Performs a specialized matrix multiplication and addition operation.
  *
  * This function computes the matrix product of host_a (transposed) and host_b, 
  * multiplies the result by a scalar, and adds the result to host_c. The computation 
  * is distributed between the CPU and a GPU device for efficiency. The CPU computation 
  * is done using the DGEMM_PWDFT function, while the GPU computation uses the MKL's gemm function.
  *
  * @param n       Number of rows in the resulting matrix.
  * @param m       Number of columns in the resulting matrix.
  * @param npack2  Leading dimension of host_a and host_b.
  * @param nida2   Parameter for the DGEMM_PWDFT function to determine the extent of computation.
  * @param host_a  Pointer to the first input matrix.
  * @param host_b  Pointer to the second input matrix.
  * @param host_c  Pointer to the output matrix. This matrix is also used as an input for the addition operation.
  *
  * @note The matrices host_a and host_b are transferred to the device asynchronously. 
  *       The CPU starts its computation immediately after initiating the data transfer to the GPU. 
  *       Once the GPU computation is complete, the results are added to host_c, which already contains 
  *       the results of the CPU computation.
  */
    void TN_dgemm2c(int n, int m, int npack2, int nida2, double *host_a, double *host_b, double *host_c) 
    {
       constexpr double rtwo  = 2.0;
       constexpr double rone  = 1.0;
       constexpr double rmone = -1.0;
       constexpr double rzero = 0.0;
      
       int ia = fetch_dev_mem_indx(static_cast<size_t>(npack2) * n);
       int ib = fetch_dev_mem_indx(static_cast<size_t>(npack2) * m);
       int ic = fetch_dev_mem_indx(static_cast<size_t>(n) * m);
 
       syclSetMatrixAsync(npack2,n,sizeof(double),host_a,npack2,dev_mem[ia],npack2,stream[0]);
       syclSetMatrixAsync(npack2,m,sizeof(double),host_b,npack2,dev_mem[ib],npack2,stream[0]);
 
       // Start the DGEMM_PWDFT operation on the CPU
      
       stream[0]->wait();
       oneapi::mkl::blas::column_major::gemm(*stream[0], 
          matT,matN,n,m,npack2, 
          rtwo,
          dev_mem[ia],npack2, 
          dev_mem[ib],npack2, 
          rzero,
          dev_mem[ic],n);
      
       syclGetMatrixAsync(n,m,sizeof(double),dev_mem[ic],n,host_c,n,stream[0]);
       stream[0]->wait();
 
       if (nida2 > 0) 
       {
          DGEMM_PWDFT((char *) "T", (char *) "N", n, m, nida2, rmone, host_a, npack2, host_b, npack2, rzero, host_c, n);
       }
 
       inuse[ia] = false;
       inuse[ib] = false;
       inuse[ic] = false;
    }
    
 
   /**************************************
    *                                    *
    *              NT_dgemm3             *
    *                                    *
    **************************************/
   void NT_dgemm3(int m, int n, int k,
                  double alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double beta,
                  double *host_c,int ldc) 
   {
      int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)k));
      int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)k));
      int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
 
      syclSetMatrixAsync(lda,k,sizeof(double),host_a,lda,dev_mem[ia],lda,stream[0]);
      syclSetMatrixAsync(ldb,k,sizeof(double),host_b,ldb,dev_mem[ib],ldb,stream[0]);
 
      stream[0]->wait();
      oneapi::mkl::blas::column_major::gemm(*stream[0], 
         matN,matT,m,n,k, 
         alpha,
         dev_mem[ia],lda, 
         dev_mem[ib],ldb, 
         beta,
         dev_mem[ic],ldc);
 
      syclGetMatrixAsync(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]);
      stream[0]->wait();
 
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }
 
 
   /**************************************
    *                                    *
    *              TN_dgemm              *
    *                                    *
    **************************************/
   void TN_dgemm(int ne, int nprj, int npack2, double alpha, double *host_a,
                 double *host_b, double beta, double *host_c) 
   {
      // DGEMM_PWDFT((char *) "T",(char *)
      // "N",ne,nprj,npack2,alpha,host_a,npack2,host_b,npack2,beta,host_c,ne);
     
      // gdevice_TN_dgemm(nn,nprj,ng,rtwo,a,b,rzero,sum);
     
      // int ia = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) ne));
      // int ib = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) nprj));
      b_prj = host_b;
      ib_prj[0] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)nprj));
      if (tile_fac > 1)
         ib_prj[1] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)nprj));
      int ic = fetch_dev_mem_indx(((size_t)ne) * ((size_t)nprj));
     
      syclSetMatrixAsync(ne, nprj, sizeof(double), host_c, ne, dev_mem[ic], ne,
                         stream[0]);
     
      if (tile_fac > 1)
         syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                            &a_psi[tile_start2[0]], npack2, dev_mem[ia_psi[0]],
                            tile_npack2[0], stream[0]);
      syclSetMatrixAsync(tile_npack2[0], nprj, sizeof(double),
                         &b_prj[tile_start2[0]], npack2, dev_mem[ib_prj[0]],
                         tile_npack2[0], stream[0]);
     
      double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) 
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) 
         {
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &a_psi[tile_start2[ttp1]], npack2,
                               dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
            syclSetMatrixAsync(tile_npack2[ttp1], nprj, sizeof(double),
                               &b_prj[tile_start2[ttp1]], npack2,
                               dev_mem[ib_prj[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
         }
        
         stream[tt % 2]->wait();
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matT, matN, ne, nprj, tile_npack2[tt], alpha,
             dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ib_prj[tt % 2]],
             tile_npack2[tt], beta0, dev_mem[ic], ne);
         beta0 = 1.0;
      }
      stream[0]->memcpy(host_c, dev_mem[ic], ne * nprj * sizeof(double)).wait();
     
      // inuse[ia] = false;
      // inuse[ib_prj[0]] = false;
      // if (tile_fac>1) inuse[ib_prj[1]] = false;
      inuse[ic] = false;
   }
 
   /**************************************
    *                                    *
    *              T_free                *
    *                                    *
    **************************************/
   void T_free() 
   {
      inuse[ib_prj[0]] = false;
      if (tile_fac > 1)
         inuse[ib_prj[1]] = false;
   }
 
   /**************************************
    *                                    *
    *              NT_dgemm              *
    *                                    *
    **************************************/
   void NT_dgemm(int npack2, int ne, int nprj, double alpha, double *host_a,
                 double *host_b, double beta, double *host_c) 
   {
      // DGEMM_PWDFT((char *) "N",(char *)
      // "T",npack2,ne,nprj,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);
     
      int ib = fetch_dev_mem_indx(((size_t)ne) * ((size_t)nprj));
     
      syclSetMatrixAsync(ne, nprj, sizeof(double), host_b, ne, dev_mem[ib], ne,
                         stream[(tile_fac - 1) % 2]);
      syclSetMatrixAsync(tile_npack2[tile_fac - 1], ne, sizeof(double),
                         &host_c[tile_start2[tile_fac - 1]], npack2,
                         dev_mem[ia_hpsi[(tile_fac - 1) % 2]],
                         tile_npack2[tile_fac - 1], stream[(tile_fac - 1) % 2]);
      syclSetMatrixAsync(tile_npack2[tile_fac - 1], nprj, sizeof(double),
                         &host_a[tile_start2[tile_fac - 1]], npack2,
                         dev_mem[ib_prj[(tile_fac - 1) % 2]],
                         tile_npack2[tile_fac - 1], stream[(tile_fac - 1) % 2]);
      for (auto tt = tile_fac - 1; tt >= 0; --tt) 
      {
         int ttm1 = tt - 1;
         if (ttm1 >= 0) 
         {
            syclSetMatrixAsync(tile_npack2[ttm1], ne, sizeof(double),
                               &host_c[tile_start2[ttm1]], npack2,
                               dev_mem[ia_hpsi[ttm1 % 2]], tile_npack2[ttm1],
                               stream[ttm1 % 2]);
            syclSetMatrixAsync(tile_npack2[ttm1], nprj, sizeof(double),
                               &host_a[tile_start2[ttm1]], npack2,
                               dev_mem[ib_prj[ttm1 % 2]], tile_npack2[ttm1],
                               stream[ttm1 % 2]);
         }
         stream[tt % 2]->wait();
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matN, matT, tile_npack2[tt], ne, nprj, alpha,
             dev_mem[ib_prj[tt % 2]], tile_npack2[tt], dev_mem[ib], ne, beta,
             dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt]);
         syclGetMatrixAsync(tile_npack2[tt], ne, sizeof(double),
                            dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt],
                            &host_c[tile_start2[tt]], npack2, stream[tt % 2]);
      }
     
      inuse[ib] = false;
      inuse[ib_prj[0]] = false;
      if (tile_fac > 1)
         inuse[ib_prj[1]] = false;
   }
 
   /**************************************
    *                                    *
    *              MM6_dgemm             *
    *                                    *
    **************************************/
   void MM6_dgemm(int ne, double *host_s21, double *host_s12, double *host_s11,
                  double *host_sa0, double *host_sa1, double *host_st1) 
   {
      double rzero = 0.0;
      double rone = 1.0;
      int i_s21 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
      int i_s12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
      int i_s11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
      int i_sa0 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
      int i_st1 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // tmp
      int i_sa1 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input-output
     
      syclSetMatrixAsync(ne, ne, sizeof(double), host_s21, ne, dev_mem[i_s21], ne, stream[0]);
      syclSetMatrixAsync(ne, ne, sizeof(double), host_sa0, ne, dev_mem[i_sa0], ne, stream[0]);
      syclSetMatrixAsync(ne, ne, sizeof(double), host_sa1, ne, dev_mem[i_sa1], ne, stream[0]);
     
      syclSetMatrixAsync(ne, ne, sizeof(double), host_s12, ne, dev_mem[i_s12], ne, stream[1]);
      syclSetMatrixAsync(ne, ne, sizeof(double), host_s11, ne, dev_mem[i_s11], ne, stream[1]);
     
      // mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
      stream[0]->wait();
      oneapi::mkl::blas::column_major::gemm(
          *stream[0], matN, matN, ne, ne, ne, rone, dev_mem[i_s21], ne,
          dev_mem[i_sa0], ne, rone, dev_mem[i_sa1], ne);
     
      // mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
      stream[1]->wait();
      oneapi::mkl::blas::column_major::gemm(
          *stream[1], matN, matN, ne, ne, ne, rone, dev_mem[i_sa0], ne,
          dev_mem[i_s12], ne, rone, dev_mem[i_sa1], ne);
     
      // mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
      oneapi::mkl::blas::column_major::gemm(
          *stream[1], matN, matN, ne, ne, ne, rone, dev_mem[i_s11], ne,
          dev_mem[i_sa0], ne, rzero, dev_mem[i_st1], ne);
     
      // mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
      oneapi::mkl::blas::column_major::gemm(
          *stream[1], matN, matN, ne, ne, ne, rone, dev_mem[i_sa0], ne,
          dev_mem[i_st1], ne, rone, dev_mem[i_sa1], ne);
      syclGetMatrixAsync(ne, ne, sizeof(double), dev_mem[i_sa1], ne, host_sa1, ne,
                         stream[1]);
      stream[1]->wait();
     
      inuse[i_s21] = false;
      inuse[i_s12] = false;
      inuse[i_s11] = false;
      inuse[i_sa0] = false;
      inuse[i_st1] = false;
      inuse[i_sa1] = false;
   }




   /**************************************
    *                                    *
    *             isCmplxZero            *
    *                                    *
    **************************************/
   bool isCmplxZero(const double *beta)
   {
      const double EPSILON = 1e-9; // Threshold for comparison
      return std::abs(beta[0]) < EPSILON && std::abs(beta[1]) < EPSILON;
   }
      
   /**************************************
    *                                    *
    *              NN1_zgemm             *
    *                                    *
    **************************************/
   void NN1_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a, double *host_b,
                 double *beta, double *host_c) 
   {
      //ZGEMM_PWDFT((char *)"N", (char *)"N", npack, ne, ne, alpha, host_a, npack1_max,
      //            host_b, ne, beta, host_c, npack1_max);
      //

      int ib = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t)ne));

      syclSetMatrixAsync(ne, ne, sizeof(double), host_b, 2*ne, dev_mem[ib], 2*ne, stream[0]);

      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         host_a + tile_start2[0], 2*npack1_max, dev_mem[ia_psi[0]], tile_npack2[0], stream[0]);

      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double),
                         host_c + tile_start2[0], 2*npack1_max, dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]);

       // double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt)
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac)
         {
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &host_a[tile_start2[ttp1]], 2*npack1_max,
                               dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
            syclSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double),
                               &host_c[tile_start2[ttp1]], 2*npack1_max,
                               dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
                               stream[ttp1 % 2]);
         }
         stream[tt % 2]->wait();
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matN, matN, tile_npack1[tt], ne, ne, 
	     reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia_psi[tt % 2]]), tile_npack1[tt], 
	     reinterpret_cast<const std::complex<double> *>(dev_mem[ib]), ne, 
	     reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[ia_hpsi[tt % 2]]), tile_npack1[tt]);
         syclGetMatrixAsync(tile_npack2[tt], ne, sizeof(double),
                            dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt],
                            &host_c[tile_start2[tt]], 2*npack1_max, stream[tt % 2]);
      }

      stream[(tile_fac - 1) % 2]->wait();

      inuse[ib] = false;
   }

        
   /**************************************
    *                                    *
    *              CN1_zgemm             *
    *                                    *
    **************************************/
   void CN1_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) 
   {
      //ZGEMM_PWDFT((char *)"C", (char *)"N", ne, ne, npack, alpha, host_a, npack1_max,
      //            host_b, npack1_max, beta, host_c, ne);

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ic = fetch_dev_mem_indx(((size_t) 2*ne)         * ((size_t) ne));

      syclSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]);
      syclSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_b, 2*npack1_max, dev_mem[ib], 2*npack1_max, stream[0]);
      if (!isCmplxZero(beta))
         syclSetMatrixAsync(2*ne, ne, sizeof(double), host_c, 2*ne, dev_mem[ic], 2*ne, stream[0]);

      stream[0]->wait();
      oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, ne, ne, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[ic]), ne);

      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[ic], 2*ne, host_c, 2*ne, stream[0]);
      stream[0]->wait();

      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }                
        
   /**************************************
    *                                    *
    *              CN2_zgemm             *
    *                                    *
    **************************************/
   void CN2_zgemm(int ne, int nprj, int npack, int npack1_max, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) 
   {
      //ZGEMM_PWDFT((char *)"C", (char *)"N", ne, nprj, npack, alpha, host_a, npack1_max,
      //            host_b, npack1_max, beta, host_c, ne);

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) nprj));
      int ic = fetch_dev_mem_indx(((size_t) 2*ne)     * ((size_t) nprj));

      syclSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]);
      syclSetMatrixAsync(2*npack1_max, nprj, sizeof(double), host_b, 2*npack1_max, dev_mem[ib], 2*npack1_max, stream[0]);
      if (!isCmplxZero(beta))
         syclSetMatrixAsync(2*ne, nprj, sizeof(double), host_c, 2*ne, dev_mem[ic], 2*ne, stream[0]);

      stream[0]->wait();
      oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, ne, nprj, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[ic]), ne);

      syclGetMatrixAsync(2*ne, nprj, sizeof(double), dev_mem[ic], 2*ne, host_c, 2*ne, stream[0]);
      stream[0]->wait();

      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }

   /**************************************
    *                                    *
    *           CN2_stride_zgemm         *
    *                                    *
    **************************************/
   void CN2_stride_zgemm(int ne, int nprj, int npack, int npack1_max, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) 
   {
      ZGEMM_PWDFT((char *)"C", (char *)"N", ne, nprj, npack, alpha, host_a, npack1_max,
                  host_b, npack1_max, beta, host_c, ne);
      
   }

   /**************************************
    *                                    *
    *              NC2_zgemm             *
    *                                    *
    **************************************/
   void NC2_zgemm(int npack1_max, int npack, int ne, int nprj,  double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) 
   {
      //ZGEMM_PWDFT((char *)"N", (char *)"C", npack,ne,nprj, alpha, host_a, npack1_max,
      //            host_b, ne, beta, host_c, npack1_max);

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) nprj));
      int ib = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) nprj));
      int ic = fetch_dev_mem_indx(((size_t) 2*npack1_max)     * ((size_t) ne));

      syclSetMatrixAsync(2*npack1_max, nprj, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]);
      syclSetMatrixAsync(2*ne, nprj, sizeof(double), host_b, 2*ne, dev_mem[ib], 2*ne, stream[0]);
      if (!isCmplxZero(beta))
         syclSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_c, 2*npack1_max, dev_mem[ic], 2*npack1_max, stream[0]);

      stream[0]->wait();
      oneapi::mkl::blas::column_major::gemm(
             *stream[0], matN, matC, npack, ne, nprj, 
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]), ne,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[ic]), npack1_max);

      syclGetMatrixAsync(2*npack1_max, ne, sizeof(double), dev_mem[ic], 2*npack1_max, host_c, 2*npack1_max, stream[0]);
      stream[0]->wait();

      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }


   /**************************************
    *                                    *
    *           NC2_stride_zgemm         *
    *                                    *
    **************************************/
   void NC2_stride_zgemm(int npack1_max, int npack, int ne, int nprj,  double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) 
   {
      ZGEMM_PWDFT((char *)"N", (char *)"C", npack,ne,nprj, alpha, host_a, npack1_max,
                  host_b, ne, beta, host_c, npack1_max);
   }

   /**************************************
    *                                    *
    *            CN4_zgemm               *
    *                                    *
    **************************************/
   void CN4_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                 double *host_b, double *beta, double *host_caa,
                 double *host_cab, double *host_cba, double *host_cbb) 
   {
      int one = 1;
      int shift1 = 0;
      int mshift1 = 0;

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int icaa = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icab = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icba = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icbb = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));

      syclSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]);
      syclSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_b, 2*npack1_max, dev_mem[ib], 2*npack1_max, stream[1]);

      if (!isCmplxZero(beta))
      {
         syclSetMatrixAsync(2*ne, ne, sizeof(double), host_caa, 2*ne, dev_mem[icaa], 2*ne, stream[0]);
         syclSetMatrixAsync(2*ne, ne, sizeof(double), host_cab, 2*ne, dev_mem[icab], 2*ne, stream[1]);
         syclSetMatrixAsync(2*ne, ne, sizeof(double), host_cba, 2*ne, dev_mem[icba], 2*ne, stream[1]);
         syclSetMatrixAsync(2*ne, ne, sizeof(double), host_cbb, 2*ne, dev_mem[icbb], 2*ne, stream[1]);
      }

      stream[0]->wait();
      for (auto k = 1; k <= ne; ++k)
      {
	 /*
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_a, npack1_max,
                     host_a + shift1, npack1_max,
                     beta,
                     host_caa + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_a, npack1_max,
                     host_b + shift1, npack1_max,
                     beta,
                     host_cab + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_b, npack1_max,
                     host_a + shift1, npack1_max,
                     beta,
                     host_cba + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_b, npack1_max,
                     host_b + shift1, npack1_max,
                     beta,
                     host_cbb + mshift1, k);
	 */
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, k, one, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]+shift1), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[icaa]+mshift1), k);

         shift1  += 2*npack1_max;
         mshift1 += 2*ne;
      }
      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[icaa], 2*ne, host_caa, 2*ne, stream[0]);

      shift1 = 0;
      mshift1 = 0;
      stream[1]->wait();
      for (auto k = 1; k <= ne; ++k)
      {
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, k, one, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]+shift1), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[icab]+mshift1), k);
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, k, one, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]+shift1), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[icba]+mshift1), k);
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, k, one, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]+shift1), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[icbb]+mshift1), k);

         shift1  += 2*npack1_max;
         mshift1 += 2*ne;
      }
      stream[0]->wait();
      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[icab], 2*ne, host_cab, 2*ne, stream[1]);
      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[icba], 2*ne, host_cba, 2*ne, stream[1]);
      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[icbb], 2*ne, host_cbb, 2*ne, stream[1]);
      stream[1]->wait();

      inuse[ia] = false;
      inuse[ib] = false;
      inuse[icaa] = false;
      inuse[icab] = false;
      inuse[icba] = false;
      inuse[icbb] = false;
   }

   /**************************************
    *                                    *
    *            CN3__zgemm              *
    *                                    *
    **************************************/
   void CN3_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_caa,
                  double *host_cab, double *host_cbb) 
   {
      int one = 1;   
      int shift1 = 0;
      int mshift1 = 0;

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int icaa = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icab = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icbb = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));

      syclSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]);
      syclSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_b, 2*npack1_max, dev_mem[ib], 2*npack1_max, stream[1]);

      if (!isCmplxZero(beta))
      {
         syclSetMatrixAsync(2*ne, ne, sizeof(double), host_caa, 2*ne, dev_mem[icaa], 2*ne, stream[0]);
         syclSetMatrixAsync(2*ne, ne, sizeof(double), host_cab, 2*ne, dev_mem[icab], 2*ne, stream[1]);
         syclSetMatrixAsync(2*ne, ne, sizeof(double), host_cbb, 2*ne, dev_mem[icbb], 2*ne, stream[1]);
      }

                     
      stream[0]->wait();
      for (auto k = 1; k <= ne; ++k)  
      {              
         /*
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha, 
                     host_a, npack1_max,
                     host_a + shift1, npack1_max,
                     beta,  
                     host_caa + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_a, npack1_max,
                     host_b + shift1, npack1_max,
                     beta,
                     host_cab + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_b, npack1_max,
                     host_b + shift1, npack1_max,
                     beta,
                     host_cbb + mshift1, k);
         */
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, k, one, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]+shift1), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[icaa]+mshift1), k);

         shift1 += 2*npack1_max;
         mshift1 += 2*ne;
      }
      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[icaa], 2*ne, host_caa, 2*ne, stream[0]);

      shift1 = 0;
      mshift1 = 0;
      stream[1]->wait();
      for (auto k = 1; k <= ne; ++k)  
      {              
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, k, one, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ia]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]+shift1), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[icab]+mshift1), k);
         oneapi::mkl::blas::column_major::gemm(
             *stream[0], matC, matN, k, one, npack,
             reinterpret_cast<const std::complex<double> *>(alpha),
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]), npack1_max,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ib]+shift1), npack1_max,
             reinterpret_cast<const std::complex<double> *>(beta),
             reinterpret_cast<      std::complex<double> *>(dev_mem[icbb]+mshift1), k);

         shift1 += 2*npack1_max;
         mshift1 += 2*ne;
      }
      stream[0]->wait();
      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[icab], 2*ne, host_cab, 2*ne, stream[1]);
      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[icbb], 2*ne, host_cbb, 2*ne, stream[1]);
      stream[1]->wait();

      inuse[ia] = false;
      inuse[ib] = false;
      inuse[icaa] = false;
      inuse[icab] = false;
      inuse[icbb] = false;
   }



   /**************************************
    *                                    *
    *             WW6_zgemm              *
    *                                    *
    **************************************/
   void WW6_zgemm(int ne, double *host_s21, double *host_s12, double *host_s11,
                  double *host_sa0, double *host_sa1, double *host_st1) 
   {
      double rone[2]  = {1.0,0.0};
      double rzero[2] = {0.0,0.0};

      int is11 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int is12 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int is21 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int isa0 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int isa1 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int ist1 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));

      syclSetMatrixAsync(2*ne, ne, sizeof(double), host_s21, 2*ne, dev_mem[is21], 2*ne, stream[0]);
      syclSetMatrixAsync(2*ne, ne, sizeof(double), host_sa0, 2*ne, dev_mem[isa0], 2*ne, stream[0]);
      syclSetMatrixAsync(2*ne, ne, sizeof(double), host_sa1, 2*ne, dev_mem[isa1], 2*ne, stream[0]);

      syclSetMatrixAsync(2*ne, ne, sizeof(double), host_s12, 2*ne, dev_mem[is12], 2*ne, stream[1]);
      syclSetMatrixAsync(2*ne, ne, sizeof(double), host_s11, 2*ne, dev_mem[is11], 2*ne, stream[2]);

      stream[0]->wait();
      
      // www_Multiply1(ms, s21, sa0, 1.0, sa1, 1.0);
      //ZGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_s21, ne,
      //            host_sa0, ne, rone, host_sa1, ne);
      oneapi::mkl::blas::column_major::gemm(
             *stream[0], matN, matN, ne, ne, ne,
             reinterpret_cast<const std::complex<double> *>(rone),
             reinterpret_cast<const std::complex<double> *>(dev_mem[is21]), ne,
             reinterpret_cast<const std::complex<double> *>(dev_mem[isa0]), ne,
             reinterpret_cast<const std::complex<double> *>(rone),
             reinterpret_cast<      std::complex<double> *>(dev_mem[isa1]), ne);
     
      stream[1]->wait();

      // www_Multiply2(ms, sa0, s12, 1.0, sa1, 1.0);
      //ZGEMM_PWDFT((char *)"C", (char *)"N", ne, ne, ne, rone, host_sa0, ne,
      //            host_s12, ne, rone, host_sa1, ne);
      oneapi::mkl::blas::column_major::gemm(
             *stream[0], matN, matN, ne, ne, ne,
             reinterpret_cast<const std::complex<double> *>(rone),
             reinterpret_cast<const std::complex<double> *>(dev_mem[isa0]), ne,
             reinterpret_cast<const std::complex<double> *>(dev_mem[is12]), ne,
             reinterpret_cast<const std::complex<double> *>(rone),
             reinterpret_cast<      std::complex<double> *>(dev_mem[isa1]), ne);

      stream[2]->wait();

      // www_Multiply3(ms, s11, sa0, 1.0, st1, 0.0);
      //ZGEMM_PWDFT((char *)"N", (char *)"C", ne, ne, ne, rone, host_s11, ne,
      //             host_sa0, ne, rzero, host_st1, ne);
      oneapi::mkl::blas::column_major::gemm(
             *stream[0], matN, matN, ne, ne, ne,
             reinterpret_cast<const std::complex<double> *>(rone),
             reinterpret_cast<const std::complex<double> *>(dev_mem[is11]), ne,
             reinterpret_cast<const std::complex<double> *>(dev_mem[isa0]), ne,
             reinterpret_cast<const std::complex<double> *>(rzero),
             reinterpret_cast<      std::complex<double> *>(dev_mem[ist1]), ne);
                    
      // www_Multiply1(ms, sa0, st1, 1.0, sa1, 1.0);
      //ZGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_sa0, ne,
      //             host_st1, ne, rone, host_sa1, ne);
      oneapi::mkl::blas::column_major::gemm(
             *stream[0], matN, matN, ne, ne, ne,
             reinterpret_cast<const std::complex<double> *>(rone),
             reinterpret_cast<const std::complex<double> *>(dev_mem[isa0]), ne,
             reinterpret_cast<const std::complex<double> *>(dev_mem[ist1]), ne,
             reinterpret_cast<const std::complex<double> *>(rone),
             reinterpret_cast<      std::complex<double> *>(dev_mem[isa1]), ne);

      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[ist1], 2*ne, host_st1, 2*ne, stream[0]);
      syclGetMatrixAsync(2*ne, ne, sizeof(double), dev_mem[isa1], 2*ne, host_sa1, 2*ne, stream[0]);
      stream[0]->wait();

      inuse[is11] = false;
      inuse[is12] = false;
      inuse[is21] = false;
      inuse[isa0] = false;
      inuse[isa1] = false;
      inuse[ist1] = false;
   } 



/*   void NN_zgemm(int m, int n, int k,
                 std::complex<double> *alpha,
                 std::complex<double> *host_a, int lda,
                 std::complex<double> *host_b, int ldb,
                 std::complex<double> *beta,
                 std::complex<double> *host_c, int ldc)
   {
      // Calculate indices for device memory
      int ia = fetch_dev_mem_indx(static_cast<size_t>(2*lda) * static_cast<size_t>(k));
      int ib = fetch_dev_mem_indx(static_cast<size_t>(2*ldb) * static_cast<size_t>(n));
      int ic = fetch_dev_mem_indx(static_cast<size_t>(2*ldc) * static_cast<size_t>(n));
 
      // Cast complex pointers to oneMKL complex data type
      auto complex_A = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_a);
      auto complex_B = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_b);
      auto complex_C = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_c);
 
      // Perform asynchronous memory transfer to device
      syclSetMatrixAsync(lda, k, sizeof(std::complex<double>), host_a, lda, dev_mem[ia], lda, stream[0]);
      syclSetMatrixAsync(ldb, n, sizeof(std::complex<double>), host_b, ldb, dev_mem[ib], ldb, stream[0]);
 
      // Wait for the memory transfer to complete
      stream[0]->wait();
 
      // Perform complex matrix multiplication
      oneapi::mkl::blas::column_major::gemm(*stream[0],
                                            matN, matN,
                                            m, n, k,
                                            alpha,
                                            complex_A, lda,
                                            complex_B, ldb,
                                            beta,
                                            complex_C, ldc);
 
      // Perform asynchronous memory transfer back to host
      syclGetMatrixAsync(ldc, n, sizeof(std::complex<double>), dev_mem[ic], ldc, host_c, ldc, stream[0]);
 
      // Wait for the memory transfer to complete
      stream[0]->wait();
 
      // Mark memory as unused
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   } */

   /**************************************
    *                                    *
    *              NN_zgemm              *
    *                                    *
    **************************************/
   void NN_zgemm(int m, int n, int k,
                 double *alpha, 
                 double *host_a, int lda,
                 double *host_b, int ldb,
                 double *beta,
                 double *host_c,int ldc) 
   {
      ZGEMM_PWDFT((char *)"N", (char *)"N", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   }   
   


   /*
   void CN_zgemm(int m, int n, int k,
                 std::complex<double> *alpha,
                 std::complex<double> *host_a, int lda,
                 std::complex<double> *host_b, int ldb,
                 std::complex<double> *beta,
                 std::complex<double> *host_c, int ldc)
   {
      // Calculate indices for device memory
      int ia = fetch_dev_mem_indx(static_cast<size_t>(2 * lda) * static_cast<size_t>(k));
      int ib = fetch_dev_mem_indx(static_cast<size_t>(2 * ldb) * static_cast<size_t>(n));
      int ic = fetch_dev_mem_indx(static_cast<size_t>(2 * ldc) * static_cast<size_t>(n));
 
      // Cast complex pointers to oneMKL complex data type
      auto complex_A = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_a);
      auto complex_B = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_b);
      auto complex_C = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_c);
 
      // Asynchronously set matrices A and B on the device
      syclSetMatrixAsync(lda, k, sizeof(std::complex<double>), host_a, lda, dev_mem[ia], lda, stream[0]);
      syclSetMatrixAsync(ldb, n, sizeof(std::complex<double>), host_b, ldb, dev_mem[ib], ldb, stream[0]);
 
      // Wait for the asynchronous operations to complete
      stream[0]->wait();
 
      // Perform complex matrix multiplication on the device with conjugate transpose of A
      oneapi::mkl::blas::column_major::gemm(*stream[0],
                                            matC, matN,
                                            m, n, k,
                                            alpha,
                                            complex_A, lda,
                                            complex_B, ldb,
                                            beta,
                                            complex_C, ldc);
     
      // Asynchronously retrieve matrix C from the device
      syclGetMatrixAsync(ldc, n, sizeof(std::complex<double>), dev_mem[ic], ldc, host_c, ldc, stream[0]);
    
      // Wait for the asynchronous operation to complete
      stream[0]->wait();
                  
      // Mark device memory as unused
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   } */

   /**************************************
    *                                    *
    *              CN_zgemm              *
    *                                    *
    **************************************/
   void CN_zgemm(int m, int n, int k,
                  double *alpha, 
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double *beta,
                  double *host_c,int ldc) {
      ZGEMM_PWDFT((char *)"C", (char *)"N", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   }     
      
   /*
   void NC_zgemm(int m, int n, int k,
                 std::complex<double> *alpha,
                 std::complex<double> *host_a, int lda,
                 std::complex<double> *host_b, int ldb,
                 std::complex<double> *beta,
                 std::complex<double> *host_c, int ldc)
   {
       // Calculate indices for device memory
       int ia = fetch_dev_mem_indx(static_cast<size_t>(2*lda) * static_cast<size_t>(k));
       int ib = fetch_dev_mem_indx(static_cast<size_t>(2*ldb) * static_cast<size_t>(n));
       int ic = fetch_dev_mem_indx(static_cast<size_t>(2*ldc) * static_cast<size_t>(n));
 
       // Cast complex pointers to oneMKL complex data type
       auto complex_A = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_a);
       auto complex_B = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_b);
       auto complex_C = reinterpret_cast<oneapi::mkl::blas::complex_double *>(host_c);
 
       // Asynchronously set matrices A and B on the device
       syclSetMatrixAsync(lda, k, sizeof(std::complex<double>), host_a, lda, dev_mem[ia], lda, stream[0]);
       syclSetMatrixAsync(ldb, n, sizeof(std::complex<double>), host_b, ldb, dev_mem[ib], ldb, stream[0]);
 
       // Wait for the asynchronous operations to complete
       stream[0]->wait();
 
       // Perform complex matrix multiplication on the device
       oneapi::mkl::blas::column_major::gemm(*stream[0],
                                             matN, matC,
                                             m, n, k,
                                             alpha,
                                             complex_A, lda,
                                             complex_B, ldb,
                                             beta,
                                             complex_C, ldc);

       // Asynchronously retrieve matrix C from the device
       syclGetMatrixAsync(ldc, n, sizeof(std::complex<double>), dev_mem[ic], ldc, host_c, ldc, stream[0]);
 
       // Wait for the asynchronous operation to complete
       stream[0]->wait();
                   
       // Mark device memory as unused
       inuse[ia] = false;
       inuse[ib] = false;
       inuse[ic] = false;
   }
*/
   /**************************************
    *                                    *
    *              NC_zgemm              *
    *                                    *
    **************************************/
   void NC_zgemm(int m, int n, int k,
                  double *alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double *beta,
                  double *host_c,int ldc) 
   {
      ZGEMM_PWDFT((char *)"N", (char *)"C", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   } 

 
   /********************/
   /* psi_dev functions*/
   /********************/

   /**************************************
    *                                    *
    *             psi_alloc              *
    *                                    *
    **************************************/
   void psi_alloc(int npack1, int ne, int tile_fac0 = 1) 
   {
      tile_fac = tile_fac0;
     
      tile_npack2_max = (((2 * npack1) % tile_fac) == 0)
                            ? (2 * npack1) / tile_fac
                            : (2 * npack1) / tile_fac + 1;
      // for (auto i=0; i<tile_fac; ++i) tile_npack2[i] =
      // (i<((2*npack1)%tile_fac)) ? (2*npack1)/tile_fac+1 : (2*npack1)/tile_fac;
      //
      for (auto i = 0; i < tile_fac; ++i)
         tile_npack1[i] = (npack1) / tile_fac;
      for (auto i = 0; i < ((npack1) % tile_fac); ++i)
         tile_npack1[i] += 1;

      for (auto i = 0; i < tile_fac; ++i)
         tile_npack2[i] = (2 * npack1) / tile_fac;
      for (auto i = 0; i < ((2 * npack1) % tile_fac); ++i)
         tile_npack2[i] += 1;
     
      tile_start1[0] = 0;
      tile_start2[0] = 0;
      for (auto i = 1; i < tile_fac; ++i)
      {
         tile_start1[i] = tile_start1[i - 1] + tile_npack1[i - 1];
         tile_start2[i] = tile_start2[i - 1] + tile_npack2[i - 1];
      }
     
      ia_psi[0] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)ne));
      ia_hpsi[0] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)ne));
     
      if (tile_fac > 1) 
      {
         ia_psi[1] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)ne));
         ia_hpsi[1] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)ne));
      }
      if (DEBUG_IO) std::cout << "Into psi_alloc, tile_factor = " << tile_fac << " ndev_mem=" << ndev_mem << std::endl;
   }

   /**************************************
    *                                    *
    *             psi_dealloc            *
    *                                    *
    **************************************/
   void psi_dealloc() 
   {
      inuse[ia_psi[0]] = false;
      inuse[ia_hpsi[0]] = false;
      if (tile_fac > 1) 
      {
         inuse[ia_psi[1]] = false;
         inuse[ia_hpsi[1]] = false;
      }
   }

   /**************************************
    *                                    *
    *          psi_copy_host2gpu         *
    *                                    *
    **************************************/
   void psi_copy_host2gpu(int npack1, int ne, double *psi) 
   {
      a_psi = psi;
      syclSetMatrixAsync(tile_npack2[0], ne, sizeof(double), psi, 2 * npack1,
                         dev_mem[ia_psi[0]], tile_npack2[0], stream[0]);
   }

   /**************************************
    *                                    *
    *          hpsi_copy_host2gpu        *
    *                                    *
    **************************************/
   void hpsi_copy_host2gpu(int npack1, int ne, double *hpsi) 
   {
      int tt = tile_fac - 1;
      a_hpsi = hpsi;
      syclSetMatrixAsync(
          tile_npack2[tt], ne, sizeof(double), &hpsi[tile_start2[tt]], 2 * npack1,
          dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], stream[tt % 2]);
   }

   /**************************************
    *                                    *
    *          psi_copy_gpu2host         *
    *                                    *
    **************************************/
   void psi_copy_gpu2host(int npack1, int ne, double *psi) 
   {
      if (tile_fac == 1) 
      {
         stream[0]->wait();
         syclGetMatrixAsync(tile_npack2[0], ne, sizeof(double), dev_mem[ia_psi[0]],
                            tile_npack2[0], psi, 2*npack1, stream[0]);
      }
   }

   /**************************************
    *                                    *
    *         hpsi_copy_gpu2host         *
    *                                    *
    **************************************/
   void hpsi_copy_gpu2host(int npack1, int ne, double *hpsi) 
   {
      stream[0]->wait();
   }

 
   /* fft functions*/

   /**************************************
    *                                    *
    *         batch_fft_init             *
    *                                    *
    **************************************/
   int batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) 
   {
      
      desc_r_x[fftcount] = new desc_real_t(nx);
      desc_r_x[fftcount]->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                        nq1);
      desc_r_x[fftcount]->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, nx + 2);
      desc_r_x[fftcount]->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, nx/2 + 1);

      desc_x[fftcount] = new desc_cmplx_t(nx);
      desc_x[fftcount]->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, nq2);
      desc_x[fftcount]->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, nx);
      desc_x[fftcount]->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, nx);

      desc_y[fftcount] = new desc_cmplx_t(ny);
      desc_y[fftcount]->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, nq2);
      desc_y[fftcount]->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, ny);
      desc_y[fftcount]->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, ny);
     
      desc_z[fftcount] = new desc_cmplx_t(nz);
      desc_z[fftcount]->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                        nq3);
      desc_z[fftcount]->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, nz);
      desc_z[fftcount]->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, nz);
     
      desc_x[fftcount]->commit(*stream[0]);
      desc_y[fftcount]->commit(*stream[0]);
      desc_z[fftcount]->commit(*stream[0]);
     
      int tag = fftcount;
      ++fftcount;
     
      return tag;
   }
 
   /**************************************
    *                                    *
    *       batch_fft_pipeline_init      *
    *                                    *
    **************************************/
   void batch_fft_pipeline_mem_init(const int nstages, const int n2ft3d) {
      ifft_n = nstages;
 
      // allocate memory and cuda streams
      for (auto i=0; i<ifft_n; ++i)
         ifft_dev[i] = fetch_dev_mem_indx(((size_t)n2ft3d));
   }
 
   /**************************************
    *                                    *
    *            batch_fft_end           *
    *                                    *
    **************************************/
   void batch_fft_end(const int tag) 
   {
     delete desc_r_x[tag];
     delete desc_x[tag];
     delete desc_y[tag];
     delete desc_z[tag];
     --fftcount;
 
     ndev_mem = 0;
   }


   /**************************************
    *                                    *
    *            batch_rfft_tmpx         *
    *                                    *
    **************************************/
   void batch_rfftx_tmpx(bool forward, int nx, int nq, int n2ft3d, double *a, double *tmpx)
   {                                 
      int nxh2 = nx + 2;
      if (forward)
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)                                         
         {                                                                     
            drfftf_(&nx, a + indx, tmpx);                                      
            indx += nxh2;
         }
         indx = 1;
         for (auto j = 0; j < (nq); ++j)
         {
            for (auto i = nx; i >= 2; --i)                                     
            {
               a[indx + i - 1] = a[indx + i - 2];                              
            }
            a[indx + 1 - 1] = 0.0;
            a[indx + nx + 1 - 1] = 0.0;
            indx += nxh2;
         }                                                                     
      }                                                                        
      else
      {
         int indx = 1;
         for (auto j = 0; j < nq; ++j)                                         
         {
            for (auto i = 2; i <= nx; ++i)
               a[indx + i - 2] = a[indx + i - 1];
            indx += nxh2;
         }
         indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            drfftb_(&nx, a + indx, tmpx);
            indx += nxh2;
         }
      }
   }


   /**************************************
    *                                    *
    *            batch_cfftx_tmpx        *
    *                                    *
    **************************************/
   void batch_cfftx_tmpx(bool forward, int nx, int nq, int n2ft3d, double *a, double *tmpx)
   {
      if (forward)
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            dcfftf_(&nx, a + indx, tmpx);
            indx += (2*nx);
         }
      }
      else
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            dcfftb_(&nx, a + indx, tmpx);
            indx += (2*nx);
         }
      }
   }

   /**************************************
    *                                    *
    *            batch_cffty_tmpy        *
    *                                    *
    **************************************/
   void batch_cffty_tmpy(bool forward, int ny, int nq, int n2ft3d, double *a, double *tmpy)
   {
      if (forward)
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            dcfftf_(&ny, a + indx, tmpy);
            indx += (2 * ny);
         }
      }
      else
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            dcfftb_(&ny, a + indx, tmpy);
            indx += (2 * ny);
         }
      }
   }
 
 
   /**************************************
    *                                    *
    *       batch_cffty_tmpy_zero        *
    *                                    *
    **************************************/
   void batch_cffty_tmpy_zero(bool forward, int ny, int nq, int nffts, int n2ft3d, double *a, double *tmpy, bool *zero)
   {
      if (forward)
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            if (!zero[q])
               dcfftf_(&ny, a + indx, tmpy);
            indx += (2 * ny);
         }
      }
      else
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            if (!zero[q])
               dcfftb_(&ny, a + indx, tmpy);
            indx += (2 * ny);
         }
      }
   }
 
 
 
   /**************************************
    *                                    *
    *          batch_cfftz_tmpz          *
    *                                    *
    **************************************/
   void batch_cfftz_tmpz(bool forward, int nz, int nq, int n2ft3d, double *a, double *tmpz) 
   {
      if (forward) 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            dcfftf_(&nz, a + indx, tmpz);
            indx += (2 * nz);
         }
      } 
      else 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            dcfftb_(&nz, a + indx, tmpz);
            indx += (2 * nz);
         }
      }
   }
 
 
   /**************************************
    *                                    *
    *        batch_cfftz_tmpz_zero       *
    *                                    *
    **************************************/
   void batch_cfftz_tmpz_zero(bool forward, int nz, int nq, int nffts, int n2ft3d, double *a, double *tmpz, bool *zero)
   {
      if (forward)
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            if (!zero[q])
               dcfftf_(&nz, a + indx, tmpz);
            indx += (2 * nz);
         }
      }
      else
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q)
         {
            if (!zero[q])
               dcfftb_(&nz, a + indx, tmpz);
            indx += (2 * nz);
         }
      }
   }
 
 
   
 
 
 
   // routines below need to be made into sycl or removed
 
   /**************************************
    *                                    *
    *            eigsrt_device           *
    *                                    *
    **************************************/
   static void eigsrt_device(double *D, double *V, int n) 
   {
      int i, j, k;
      double p;
     
      for (i = 0; i < (n - 1); ++i) 
      {
         k = i;
         p = D[i];
         for (j = i + 1; j < n; ++j)
            if (D[j] >= p) 
            {
               k = j;
               p = D[j];
            }
        
         if (k != i) 
         {
            D[k] = D[i];
            D[i] = p;
            for (j = 0; j < n; ++j) 
            {
               p = V[j + i * n];
               V[j + i * n] = V[j + k * n];
               V[j + k * n] = p;
            }
         }
      }
   }

 
   /**************************************
    *                                    *
    *           NN_eigensolver           *
    *                                    *
    **************************************/
   void NN_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) 
   {
      int n, ierr;
      int nn = ne[0] * ne[0] + 14;
      double xmp1[nn];
      // double *xmp1 = new (std::nothrow) double[nn]();
     
      int shift1 = 0;
      int shift2 = 0;
      for (int ms = 0; ms < ispin; ++ms) 
      {
         n = ne[ms];
        
         // eigen_(&n,&n,&hml[shift2],&eig[shift1],xmp1,&ierr);
         //  d3db::parall->Barrier();
         EIGEN_PWDFT(n, host_hml + shift2, host_eig + shift1, xmp1, nn, ierr);
         // if (ierr != 0) throw std::runtime_error(std::string("NWPW Error:
         // EIGEN_PWDFT failed!"));
        
         eigsrt_device(host_eig + shift1, host_hml + shift2, n);
         shift1 += ne[0];
         shift2 += ne[0] * ne[0];
      }
   }

   /**************************************
    *                                    *
    *           WW_eigensolver           *
    *                                    *
    **************************************/
   void WW_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) 
   {
      int n, ierr;
      int nn = ne[0] * ne[0] + 14;
      double xmp1[nn];
      double rmp1[nn];
      // double *xmp1 = new (std::nothrow) double[nn]();
     
      int shift1 = 0;
      int shift2 = 0;
      for (int ms=0; ms<ispin; ++ms) 
      {
         n = ne[ms];
        
         // eigen_(&n,&n,&hml[shift2],&eig[shift1],xmp1,&ierr);
         //  d3db::parall->Barrier();
         MKL_Complex16* hml = reinterpret_cast<MKL_Complex16*>(host_hml+shift2);
         ZEIGEN_PWDFT(n, hml, host_eig + shift1, xmp1, nn, rmp1, ierr);

         //ZEIGEN_PWDFT(n, host_hml + shift2, host_eig + shift1, xmp1, nn, rmp1, ierr);
         // if (ierr != 0) throw std::runtime_error(std::string("NWPW Error:
         // EIGEN_PWDFT failed!"));
       
         //eigsrt_device(host_eig + shift1, host_hml + shift2, n);
         shift1 += 2*ne[0];
         shift2 += 4*ne[0]*ne[0];
      } 
   }


   ////////////////////////// special complex-complex fft ////////////////////////////

typedef std::complex<double> complex_t;
         
   // Define a constants for the radix values
static constexpr int radix_values[] = {17, 16, 11, 9, 8, 7, 6, 5, 4, 3, 2};
//static constexpr int radix_values[] = {17, 11, 9, 8, 7, 6, 5, 4, 3, 2};
//static constexpr int radix_values[] = {17, 11, 9, 8, 7, 6, 5, 4, 3, 2};
            
   /**************************************
    *                                    *
    *             fft_radix              *
    *                                    *
    **************************************/
   inline static void fft_radix(const int n, const int s, bool eo, const int radix, const complex_t* twiddle, complex_t* x, complex_t *y)
   {
      if (n == 1) 
      { 
         // Use std::copy to efficiently copy elements from x to y
         if (eo) std::copy(x, x + s, y); 
         //if (eo) std::memcpy(y, x, s * sizeof(double));
         return;
      }

      const int m = n/radix;
      complex_t Atwiddle[radix*radix];
 
      // Precompute twiddle factors for matrix multiplication
      for (int r2=0; r2<radix; ++r2)
      for (int r1=0; r1<radix; ++r1)
         Atwiddle[r1+r2*radix] = twiddle[(r1*r2)%radix];
 
      for (int p=0; p<m; ++p) 
      {
         for (int q=0; q<s; ++q) 
         {
            complex_t* x1 = x+q+s*p;
            complex_t* y1 = y+q+s*radix*p;
 
            // Initialize the output y1 vector
            for (int r1=0; r1<radix; ++r1) 
               y1[s * r1] = 0.0;
 
            // Matrix-vector multiplication
            for (int r2 = 0; r2 < radix; ++r2)
            for (int r1 = 0; r1 < radix; ++r1)
               y1[s*r1] += x1[s*m*r2] * Atwiddle[r1 + r2*radix];
               //y1[s*r1] += x1[s*m*r2] * twiddle[(r1*r2)%radix];
 
            // Apply phase factor to each result
            for (int r1 = 0; r1 < radix; ++r1)
            {
               //complex_t fac = std::pow(twiddle[radix + r1], p);
               complex_t fac = complex_t(1.0, 0.0);
               for (int ps = 0; ps < p; ++ps)
                  fac *= twiddle[radix + r1];
 
               y1[s * r1] *= fac;
            }
         }
      }
   }

   /**************************************
    *                                    *
    *            fft_twiddle             *
    *                                    *
    **************************************/
   inline static void fft_twiddle(const int n, const complex_t* twiddle, complex_t* x) // Fourier transform
   {
      //complex_t* y = new complex_t[n];
      complex_t y[n];
      int eo = 0;
      int s  = 1;
      int nn = n;
      int nsize = 0;
      while (s<=n) 
      {
         // Identify the largest radix applicable for current nn
         int radix = 2;  // Default to radix-2
         for (int r : radix_values) {
            if (nn % r == 0) {
               radix = r;
               break;
            }
         }
 
         // Perform FFT with the determined radix
         if (eo)
            fft_radix(nn, s, eo, radix, twiddle + nsize, y, x);
         else
            fft_radix(nn, s, eo, radix, twiddle + nsize, x, y);

         nsize += 2*radix;
         nn /= radix;
         s *= radix;
         eo = !eo;  // Toggle the 'even-odd' flag
      }
   }


   /**************************************
    *                                    *
    *         set_sub_fft_twiddle        *
    *                                    *
    **************************************/
   static void set_sub_fft_twiddle(const int isgn, const int n, const int radix, complex_t *twiddle)
   {
      const double theta0 = 2*M_PI/((double) n);
      const double theta_radix = 2*M_PI/((double) radix);
 
      // Calculate radix-specific twiddle factors
      for (int r=0; r<radix; ++r)
         twiddle[r] = complex_t(cos(r*theta_radix), isgn*sin(r*theta_radix));
 
      // Calculate the main twiddle factors for the FFT
      for (int r=0; r<radix; ++r)
         twiddle[radix+r] = complex_t(cos(r*theta0), isgn*sin(r*theta0));
   }



   /**************************************
    *                                    *
    *            set_fft_twiddle         *
    *                                    *
    **************************************/
   void set_fft_twiddle(const int isgn, const int n, double *twiddle) 
   {
      complex_t* complex_twiddle = reinterpret_cast<complex_t*>(twiddle);
      int nsize = 0;
      int s = 1;
      int nn = n;
 
      while (s <= n) 
      {
         bool found = false;
         // Loop through possible radix values to find the largest factor of nn
         for (int radix : radix_values) 
         {
            if (nn % radix == 0) 
            {
               set_sub_fft_twiddle(isgn, nn, radix, complex_twiddle + nsize);
               nsize += 2*radix;
               nn /= radix;
               s *= radix;
               found = true;
               break;
            }
         }
         if (!found) break;
      }
   }

   /**************************************
    *                                    *
    *            size_fft_twiddle        *
    *                                    *
    **************************************/
   int size_fft_twiddle(const int n) 
   {
      int nsize = 0;
      int s = 1;
      int nn = n;
 
      while (s <= n) 
      {
         bool found = false;
         // Loop through possible radix values to find the largest factor of nn
         for (int radix : radix_values) 
         {
            if (nn % radix == 0) 
            {
               nsize += 2 * radix;
               nn /= radix;
               s *= radix;
               found = true;
               break;
            }
         }
         if (!found) break;
      }
 
      return nsize;
   }

   /**************************************
    *                                    *
    *              batch_cfft            *
    *                                    *
    **************************************/
   void batch_cfft(const bool forward, const int nz, const int nq, const int nffts, const int nfft3d, double *a, const double *twiddle, const double *tmpz)
   {
      // Ensure the function processes the right type of data
      // If twiddle is indeed of type complex_t, this cast is necessary
      complex_t*       complex_a       = reinterpret_cast<complex_t*>(a);
      const complex_t* complex_twiddle = reinterpret_cast<const complex_t*>(twiddle);

      int shift = nfft3d;
      //int shift = 2*nfft3d;
      int indx = 0;

      // Process each FFT batch
      for (auto q=0; q<nq; ++q)
      {
         fft_twiddle(nz, complex_twiddle, complex_a+indx);
         //(forward ? dcfftf_(&nz, a + indx, tmpz) : dcfftb_(&nz, a + indx, tmpz));

         indx += (shift);
      }
   }

   /**************************************
    *                                    *
    *          batch_cfft_zero           *
    *                                    *
    **************************************/
   void batch_cfft_zero(const bool forward, int nz, int nq, int nffts, int nfft3d, double *a, const double *twiddle, const double *tmpz, const bool *zero) 
   {
      // Ensure the function processes the right type of data
      // If twiddle is indeed of type complex_t, this cast is necessary
      complex_t*       complex_a       = reinterpret_cast<complex_t*>(a);
      const complex_t* complex_twiddle = reinterpret_cast<const complex_t*>(twiddle);

      int shift = nfft3d;
      //int shift = 2*nfft3d;
      int indx = 0;

      // Process each FFT batch
      for (auto q=0; q<nq; ++q)
      {
         if (!zero[q])
            fft_twiddle(nz, complex_twiddle, complex_a+indx);
            //(forward ? dcfftf_(&nz, a + indx, tmpz) : dcfftb_(&nz, a + indx, tmpz));

         indx += (shift);
      }
   }


   ////////////////////////// special complex-complex fft ////////////////////////////


}; // class Gdevices

} // namespace pwdft

