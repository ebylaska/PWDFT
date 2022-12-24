#pragma once

namespace pwdft {

void gdevice_TN3_dgemm(int, int, double, double *, double *, double, double *, double *, double *);
void gdevice_TN1_dgemm(int, int, double, double *, double *, double, double *);
void gdevice_TN_dgemm( int, int, int, double, double *, double *, double, double *);
void gdevice_T_free();

void gdevice_NN_dgemm(int, int, double, double *, double *, double, double *);
void gdevice_NT_dgemm(int, int, int, double, double *, double *, double, double *);

void gdevice_psi_alloc(int, int, int);
void gdevice_psi_dealloc();

void gdevice_psi_copy_host2gpu(int, int, double *);
void gdevice_hpsi_copy_host2gpu(int, int, double *);

void gdevice_psi_copy_gpu2host(int, int, double *);
void gdevice_hpsi_copy_gpu2host(int, int, double *);

void gdevice_batch_fft_init(int, int, int, int, int, int);
void gdevice_batch_fft_end();
void gdevice_batch_cfftx(bool,int,int,int,double *);
void gdevice_batch_cffty(bool,int,int,int,double *);
void gdevice_batch_cfftz(bool,int,int,int,double *);
void gdevice_batch_cfftx_tmpx(bool,int,int,int,double *,double *);
void gdevice_batch_cffty_tmpy(bool,int,int,int,double *,double *);
void gdevice_batch_cfftz_tmpz(bool,int,int,int,double *,double *);

void gdevice_batch_cffty_tmpy_zero(bool,int,int,int,double *,double *,bool *);
void gdevice_batch_cfftz_tmpz_zero(bool,int,int,int,double *,double *,bool *);

} // namespace pwdft
