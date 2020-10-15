#include	"gdevices.hpp"

static Gdevices mygdevice;


void gdevice_ffm_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *c)
{ 
   mygdevice.ffm_dgemm(npack,ne,alpha,a,b,beta,c);
}

void gdevice_fmf_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *c)
{ 
   mygdevice.fmf_dgemm(npack,ne,alpha,a,b,beta,c);
}
