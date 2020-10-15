#include	"gdevices.hpp"

static Gdevices mygdevice;


void gdevice_ffm_dgemm(const int npack, const int ne, const double alpha, const double *a, const double *b, const double beta, double *c)
{ 
   mygdevice.ffm_dgemm(npack,ne,alpha,a,b,beta,c);
}

void gdevice_fmf_dgemm(const int npack, const int ne, const double alpha, const double *a, const double *b, const double beta, double *c)
{ 
   mygdevice.fmf_dgemm(npack,ne,alpha,a,b,beta,c);
}
