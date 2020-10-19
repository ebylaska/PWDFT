#include	"gdevices.hpp"

static Gdevices mygdevice;

Gdevices::Gdevices() : ndev_mem(0) {

  cl::sycl::async_handler asyncHandler = [&](cl::sycl::exception_list eL) {
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

void gdevice_TN3_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *caa, double *cab, double *cbb)
{
  mygdevice.TN3_dgemm(npack,ne,alpha,a,b,beta,caa,cab,cbb);
}

void gdevice_TN_dgemm(int npack, int ne, int nprj, double alpha, double *a, double *b, double beta, double *c)
{
  mygdevice.TN_dgemm(npack,ne,nprj,alpha,a,b,beta,c);
}

void gdevice_NN_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *c)
{
  mygdevice.NN_dgemm(npack,ne,alpha,a,b,beta,c);
}

void gdevice_NT_dgemm(int npack, int ne, int nprj, double alpha, double *a, double *b, double beta, double *c)
{
  mygdevice.NT_dgemm(npack,ne,nprj,alpha,a,b,beta,c);
}
