#ifndef _BLAS_H_
#define _BLAS_H_


#if defined(NWPW_INTEL_MKL)
#include "mkl.h"

#define	DSCAL_PWDFT(n,alpha,a,ida)		cblas_sscal(n,alpha,a,ida);
#define DCOPY_PWDFT(n,a,ida,b,idb)              cblas_scopy(n,a,ida,b,idb)
#define DAXPY_PWDFT(n,alpha,a,ida,b,idb)        cblas_saxpy(n,alpha,a,ida,b,idb)

#define TRANSCONV(a)    ( a=="N" ?  CblasNoTrans : CblasTrans )
#define DGEMM_PWDFT(s1,s2,n,m,k,alpha,a,ida,b,idb,beta,c,idc) cblas_sgemm(CblasColMajor,TRANSCONV(s1),TRANSCONV(s2),n,m,k,alpha,a,ida,b,idb,beta,c,idc)

#define IDAMAX_PWDFT(nn,hml,one)	cblas_isamax(nn,hml,one)

#define EIGEN_PWDFT(n,hml,eig,xtmp,nn,ierr)	ierr=LAPACKE_ssyev(LAPACK_COL_MAJOR,'V','U',n,hml,n,eig)

#define	DDOT_PWDFT(n,a,ida,b,idb)	cblas_sdot(n,a,ida,b,idb);

#else


extern "C" void dcopy_(int *, float *, int *, float *, int *);
extern "C" float ddot_(int *, float *, int *, float *, int *);
extern "C" void daxpy_(int *, float *, float *, int *, float *, int *);
extern "C" void dscal_(int *, float *, float *, int *);
extern "C" void dgemm_(char *, char *, int *, int *, int *,
                       float *, 
                       float *, int *,
                       float *, int *,
                       float *,
                       float *, int *);

//extern "C" void eigen_(int *, int *, float *, float *, float *, int *);

extern "C" void dsyev_(char *, char *, int *,
                       float *, int *,
                       float *,
                       float *, int *, int*);
extern "C" int  idamax_(int *, float *, int *);


#define	DSCAL_PWDFT(n,alpha,a,ida)		sscal_(&(n),&(alpha),a,&(ida));
#define DCOPY_PWDFT(n,a,ida,b,idb)              scopy_(&(n),a,&(ida),b,&(idb))
#define DAXPY_PWDFT(n,alpha,a,ida,b,idb)        saxpy_(&(n),&(alpha),a,&(ida),b,&(idb))
#define DGEMM_PWDFT(s1,s2,n,m,k,alpha,a,ida,b,idb,beta,c,idc) sgemm_(s1,s2,&(n),&(m),&(k),&(alpha),a,&(ida),b,&(idb),&(beta),c,&(idc))

#define IDAMAX_PWDFT(nn,hml,one)	idamax_(&(nn),hml,&(one))

#define EIGEN_PWDFT(n,hml,eig,xtmp,nn,ierr) ssyev_((char *) "V",(char *) "U", &(n),hml,&(n),eig,xtmp,&(nn),&ierr)

#define	DDOT_PWDFT(n,a,ida,b,idb)	sdot_(&(n),a,&ida,b,&(idb));

#endif


#endif

