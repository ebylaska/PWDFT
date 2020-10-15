#include	"gdevices.hpp"

static Gdevices mygdevice;

void gdevice_TN3_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *caa, double *cab, double *cbb)
{ 
   mygdevice.TN3_dgemm(npack,ne,alpha,a,b,beta,caa,cab,cbb);
}

void gdevice_TN_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *c)
{ 
   mygdevice.TN_dgemm(npack,ne,alpha,a,b,beta,c);
}

void gdevice_NN_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *c)
{ 
   mygdevice.NN_dgemm(npack,ne,alpha,a,b,beta,c);
}
