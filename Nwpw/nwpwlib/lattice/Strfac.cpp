/* Strfac.C -
   Author - Eric Bylaska
*/

#include        <cstdlib>
#include        <iostream>
#include        <cstdio>
#include	<cstring>
#include        <cmath>

#include	"Strfac.hpp"


/* Constructors */


/*********************************
 *                               *
 *         Strfac::Strfac        *
 *                               *
 *********************************/

Strfac::Strfac(Ion *inion, PGrid *ingrid)
{
   int i,j,k,l,k1,k2,k3;
   int nx,ny,nz,nxh,nb,p,indx;
   int tnp,tid;
   int *ii_indx,*jj_indx,*kk_indx;

   mygrid = ingrid;
   myion  = inion;
   tnp = mygrid->np;
   tid = mygrid->taskid;
   Lattice *lattice = mygrid->lattice;

   nx = mygrid->nx;
   ny = mygrid->ny;
   nz = mygrid->nz;
   nxh = nx/2;

   for (j=0; j<3; ++j)
   for (i=0; i<3; ++i)
   {
      unitg[i+j*3] = lattice->unitg(i,j);
      unita[i+j*3] = lattice->unita(i,j);
   }

   /* allocate memory */
   wx1 = new double [2*(myion->nion)*(mygrid->nx)];
   wy1 = new double [2*(myion->nion)*(mygrid->ny)];
   wz1 = new double [2*(myion->nion)*(mygrid->nz)];

#ifdef NWPW_SYCL
   wx1_sycl = cl::sycl::malloc_device<double>(2*(myion->nion)*(mygrid->nx), *get_syclQue());
   wy1_sycl = cl::sycl::malloc_device<double>(2*(myion->nion)*(mygrid->ny), *get_syclQue());
   wz1_sycl = cl::sycl::malloc_device<double>(2*(myion->nion)*(mygrid->nz), *get_syclQue());

   i1_indx_sycl = cl::sycl::malloc_device<int>(mygrid->npack(1), *get_syclQue());
   j1_indx_sycl = cl::sycl::malloc_device<int>(mygrid->npack(1), *get_syclQue());
   k1_indx_sycl = cl::sycl::malloc_device<int>(mygrid->npack(1), *get_syclQue());
#endif

   i_indx[0] = new int[mygrid->npack(0)];
   j_indx[0] = new int[mygrid->npack(0)];
   k_indx[0] = new int[mygrid->npack(0)];
   i_indx[1] = new int[mygrid->npack(1)];
   j_indx[1] = new int[mygrid->npack(1)];
   k_indx[1] = new int[mygrid->npack(1)];

   ii_indx = new int[mygrid->nfft3d];
   jj_indx = new int[mygrid->nfft3d];
   kk_indx = new int[mygrid->nfft3d];
   for (nb=0; nb<2; ++nb)
   {
      for (k=0; k<nz;   ++k)
      for (j=0; j<ny;   ++j)
      for (i=0; i<=nxh; ++i)
      {
         p = mygrid->ijktop(i,j,k);
         if (p==tid)
         {
            indx = mygrid->ijktoindex(i,j,k);
            ii_indx[indx] = i;
            jj_indx[indx] = j;
            kk_indx[indx] = k;
         }
      }
      mygrid->i_pack(nb,ii_indx);
      mygrid->i_pack(nb,jj_indx);
      mygrid->i_pack(nb,kk_indx);
      mygrid->ii_pack_copy(nb,ii_indx,i_indx[nb]);
      mygrid->ii_pack_copy(nb,jj_indx,j_indx[nb]);
      mygrid->ii_pack_copy(nb,kk_indx,k_indx[nb]);
   }

   delete [] kk_indx;
   delete [] jj_indx;
   delete [] ii_indx;

#ifdef NWPW_SYCL
   get_syclQue()->memcpy(i1_indx_sycl, i_indx[1], mygrid->npack(1)*sizeof(int));
   get_syclQue()->memcpy(j1_indx_sycl, j_indx[1], mygrid->npack(1)*sizeof(int));
   get_syclQue()->memcpy(k1_indx_sycl, k_indx[1], mygrid->npack(1)*sizeof(int));
#endif
}


/*********************************
 *                               *
 *          Strfac::phafac       *
 *                               *
 *********************************/
void Strfac::phafac()
{
   int i,k,nxh,nyh,nzh,nx,ny,nz;
   double a,b,sw1,sw2,sw3,pi;
   double cw1x,cw2x,cw3x;
   double cw1y,cw2y,cw3y;

   pi  = 4.00*atan(1.0);

   nx = (mygrid->nx);
   ny = (mygrid->ny);
   nz = (mygrid->nz);
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;

   for (i=0; i<(myion->nion); ++i)
   {
      sw1 = unitg[0]*myion->rion1[0+3*i]
          + unitg[1]*myion->rion1[1+3*i]
          + unitg[2]*myion->rion1[2+3*i]+pi;
      sw2 = unitg[3]*myion->rion1[0+3*i]
          + unitg[4]*myion->rion1[1+3*i]
          + unitg[5]*myion->rion1[2+3*i]+pi;
      sw3 = unitg[6]*myion->rion1[0+3*i]
          + unitg[7]*myion->rion1[1+3*i]
          + unitg[8]*myion->rion1[2+3*i]+pi;

      cw1x=cos(sw1); cw1y=-sin(sw1);
      cw2x=cos(sw2); cw2y=-sin(sw2);
      cw3x=cos(sw3); cw3y=-sin(sw3);

      wx1[2*i*nx] = 1.0; wx1[2*i*nx+1] = 0.0;
      wy1[2*i*ny] = 1.0; wy1[2*i*ny+1] = 0.0;
      wz1[2*i*nz] = 1.0; wz1[2*i*nz+1] = 0.0;
      for (k=1; k<=nxh; ++k)
      {
         a = wx1[2*(k-1 + i*nx)];
         b = wx1[2*(k-1 + i*nx)+1];
         wx1[2*(k + i*nx)]   = a*cw1x - b*cw1y;
         wx1[2*(k + i*nx)+1] = a*cw1y + b*cw1x;
         wx1[2*(nx-k + i*nx)]   =  wx1[2*(k + i*nx)];
         wx1[2*(nx-k + i*nx)+1] = -wx1[2*(k + i*nx)+1];
      }
      for (k=1; k<=nyh; ++k)
      {
         a = wy1[2*(k-1 + i*ny)];
         b = wy1[2*(k-1 + i*ny)+1];
         wy1[2*(k + i*ny)]   = a*cw2x - b*cw2y;
         wy1[2*(k + i*ny)+1] = a*cw2y + b*cw2x;
         wy1[2*(ny-k + i*ny)]   =  wy1[2*(k + i*ny)];
         wy1[2*(ny-k + i*ny)+1] = -wy1[2*(k + i*ny)+1];
      }
      for (k=1; k<=nzh; ++k)
      {
         a = wz1[2*(k-1 + i*nz)];
         b = wz1[2*(k-1 + i*nz)+1];
         wz1[2*(k + i*nz)]   = a*cw3x - b*cw3y;
         wz1[2*(k + i*nz)+1] = a*cw3y + b*cw3x;
         wz1[2*(nz-k + i*nz)]   =  wz1[2*(k + i*nz)];
         wz1[2*(nz-k + i*nz)+1] = -wz1[2*(k + i*nz)+1];
      }

      wx1[2*(nxh+i*nx)] = 0.0; wx1[2*(nxh+i*nx)+1] = 0.0;
      wy1[2*(nyh+i*ny)] = 0.0; wy1[2*(nyh+i*ny)+1] = 0.0;
      wz1[2*(nzh+i*nz)] = 0.0; wz1[2*(nzh+i*nz)+1] = 0.0;

   }
#ifdef NWPW_SYCL
   get_syclQue()->memcpy(wx1_sycl, wx1,  2*(myion->nion)*(mygrid->nx)*sizeof(double));
   get_syclQue()->memcpy(wy1_sycl, wy1,  2*(myion->nion)*(mygrid->ny)*sizeof(double));
   get_syclQue()->memcpy(wz1_sycl, wz1,  2*(myion->nion)*(mygrid->nz)*sizeof(double));

   //get_syclQue()->wait();
#endif
}

/*********************************
 *                               *
 *       Strfac::strfac_pack     *
 *                               *
 *********************************/
void Strfac::strfac_pack(const int nb, const int ii, double *strx)
{
   int npack,nx,ny,nz;
   npack = mygrid->npack(nb);
   nx = mygrid->nx;
   ny = mygrid->ny;
   nz = mygrid->nz;

   const int *indxi = i_indx[nb];
   const int *indxj = j_indx[nb];
   const int *indxk = k_indx[nb];

   const double *exi = &wx1[2*ii * nx];
   const double *exj = &wy1[2*ii * ny];
   const double *exk = &wz1[2*ii * nz];

   double ai, aj, ak, bi, bj, bk;
   double c, d;
   for (int i=0; i<npack; ++i)
   {
      ai = exi[2*indxi[i]]; bi = exi[2*indxi[i]+1];
      aj = exj[2*indxj[i]]; bj = exj[2*indxj[i]+1];
      ak = exk[2*indxk[i]]; bk = exk[2*indxk[i]+1];
      c  = aj*ak - bj*bk;
      d  = aj*bk + ak*bj;
      strx[2*i]   = (ai*c - bi*d);
      strx[2*i+1] = (ai*d + bi*c);
   }
}

#ifdef NWPW_SYCL
void Strfac::strfac_pack_sycl(const int nb, const int ii, double *strx)
{
    assert(nb==1);    // following only works with nb=1

    int npack,nx,ny,nz;
    npack = mygrid->npack(nb);
    nx = mygrid->nx;
    ny = mygrid->ny;
    nz = mygrid->nz;

    const double *exi = &wx1_sycl[2*ii * nx];
    const double *exj = &wy1_sycl[2*ii * ny];
    const double *exk = &wz1_sycl[2*ii * nz];

    const int *i_idx_sycl = i1_indx_sycl;
    const int *j_idx_sycl = j1_indx_sycl;
    const int *k_idx_sycl = k1_indx_sycl;

    cl::sycl::range<1> threads(32);
    cl::sycl::range<1> blocks((npack + threads[0] - 1)/threads[0]);

    auto event = get_syclQue()->submit([&](cl::sycl::handler &cgh) {
	    auto global_range = blocks * threads;
            cgh.parallel_for<class strfac_pack>(cl::sycl::nd_range<1>(global_range, threads),
						[=](cl::sycl::nd_item<1> item)
    {
	size_t i = item.get_global_id(0);

	if ( i < npack ) {
	    int i_id = i_idx_sycl[i];
	    int j_id = j_idx_sycl[i];
	    int k_id = k_idx_sycl[i];

	    double ai = exi[2*i_id];
	    double bi = exi[2*i_id+1];
	    double aj = exj[2*j_id];
	    double bj = exj[2*j_id+1];
	    double ak = exk[2*k_id];
	    double bk = exk[2*k_id+1];

	    double c  = aj*ak - bj*bk;
	    double d  = aj*bk + ak*bj;

	    strx[2*i]   = (ai*c - bi*d);
	    strx[2*i+1] = (ai*d + bi*c);
	}

    });
        });
}

void Strfac::generate_projectors_sycl(const int nb,
				      const int ii,
				      const int nprj,
				      const int* sd_func,
				      const double* vnl,
				      double *prj)
{
    assert(nb==1); // following only works with nb=1

    int npack,nx,ny,nz;
    npack = mygrid->npack(nb);
    nx = mygrid->nx;
    ny = mygrid->ny;
    nz = mygrid->nz;

    const double *exi = &wx1_sycl[2*ii * nx];
    const double *exj = &wy1_sycl[2*ii * ny];
    const double *exk = &wz1_sycl[2*ii * nz];

    const int *i_idx_sycl = i1_indx_sycl;
    const int *j_idx_sycl = j1_indx_sycl;
    const int *k_idx_sycl = k1_indx_sycl;

    cl::sycl::range<1> threads(32);
    cl::sycl::range<1> blocks((npack + threads[0] - 1)/threads[0]);

// [[intel::kernel_args_restrict]]
    auto event = get_syclQue()->submit([&](cl::sycl::handler &cgh) {
	    auto global_range = blocks * threads;
            cgh.parallel_for<class gen_proj>(cl::sycl::nd_range<1>(global_range, threads),
	    				     [=](cl::sycl::nd_item<1> item)
    {
	size_t i = item.get_global_id(0);
	if ( i < npack ) {
	    int i_id = i_idx_sycl[i];
	    int j_id = j_idx_sycl[i];
	    int k_id = k_idx_sycl[i];

	    int exi_id = 2 * i_id;
	    int exj_id = 2 * j_id;
	    int exk_id = 2 * k_id;

	    double ai = exi[exi_id];
	    double bi = exi[exi_id+1];
	    double aj = exj[exj_id];
	    double bj = exj[exj_id+1];
	    double ak = exk[exk_id];
	    double bk = exk[exk_id+1];

	    double c  = aj*ak - bj*bk;
	    double d  = aj*bk + ak*bj;

	    double r_strfac = (ai*c - bi*d);
	    double i_strfac = (ai*d + bi*c);

	    size_t i2 = i*2;
	    for( int l=0; l<nprj; l++) {
		double vnl_val = vnl[i + l*npack];
		size_t lnpack2 = l*2*npack;
		if( sd_func[l] ) {
		    // mypneb->tcc_Mul_sycl()

		    prj[i2   + lnpack2] = vnl_val * r_strfac;
		    prj[i2+1 + lnpack2] = vnl_val * i_strfac;
		} else {
		    // mypneb->tcc_iMul_sycl()

		    prj[i2   + lnpack2] = vnl_val * -i_strfac;
		    prj[i2+1 + lnpack2] = vnl_val *  r_strfac;
		}

		//cc_pack_indot()

	    }
	}
    });
        });

    // mypneb->cc_pack_indot(1,nn,psi,prj,&sw1[l*nn]);
    // void PGrid::cc_pack_indot(const int nb, const int nn, double *a, double *b, double *sum)
    // {
    // 	int ng  = 2*(nida[nb]+nidb[nb]);
    // 	int ng0 = 2*nida[nb];
    // 	double* sw1_local = &(sw1[l * nn]);
    // 	for (int i=0; i<nn; ++i)
    // 	{
    // 	    sum[i] = 2.0 * DDOT_PWDFT(ng,&(a[i*ng]),1,b,1);
    // 	    sum[i] -= DDOT_PWDFT(ng0,&(a[i*ng]),1,b,1);

    // 	    oneapi::mkl::blas::dot(ng, &(psi[i*ng]), 1, prj, 1, &sw1_local[i]);

    // 	}

    // 	// 2 * (mypneb->neq[0]+mypneb->neq[1]) * mypneb->npack(1) // size of psi
    // 	// nn = mypneb->neq[0]+mypneb->neq[1];
    // 	// ng = 2*mypneb->neq[0]+mypneb->neq[1];
    // 	oneapi::mkl::blas::dot(nn*ng, psi, ng, prj, 1, &sw1_local[i]);
    // }

#ifdef NWPW_SYCL_ENABLE_PROFILE
    event.wait();
    auto submit_time = event.get_profiling_info<cl::sycl::info::event_profiling::command_submit>();
    auto start_time = event.get_profiling_info<cl::sycl::info::event_profiling::command_start>();
    auto end_time = event.get_profiling_info<cl::sycl::info::event_profiling::command_end>();
    auto submission_time = (start_time - submit_time) / 1000000.0f;
    auto execution_time = (end_time - start_time) / 1000000.0f;
    std::cout << "gen_projectors: " << execution_time << ", " << submission_time << std::endl;
#endif
}

/*
void Strfac::generate_projectors_all(const int nb,
				     const int nions,
				     const int* katm,
				     const int nprjMax,
				     const int* nprj,
				     const int* sd_func,
				     double** vnl,
				     double *prjAll)
{
    assert(nb==1); // following only works with nb=1

    int npack,nx,ny,nz;
    npack = mygrid->npack(nb);
    nx = mygrid->nx;
    ny = mygrid->ny;
    nz = mygrid->nz;

    const int *i_idx_sycl = i1_indx_sycl;
    const int *j_idx_sycl = j1_indx_sycl;
    const int *k_idx_sycl = k1_indx_sycl;

    const double* phafacx = wx1_sycl;
    const double* phafacy = wy1_sycl;
    const double* phafacz = wz1_sycl;

    auto event = get_syclQue()->submit([&](cl::sycl::handler &cgh) {
            cgh.parallel_for_work_group<class gen_proj>(cl::sycl::range<1>(nions), cl::sycl::range<1>(npack),
                                                        [=](cl::sycl::group<1> group)
    {
        size_t ii = group.get_id(0);
        size_t ii2 = 2*ii;

        const double* exi = &phafacx[ii2 * nx];
        const double* exj = &phafacy[ii2 * ny];
        const double* exk = &phafacz[ii2 * nz];
	const size_t ia = katm[ii];
	const double* vnl_ion = vnl[ia];
        const int n_proj = nprj[ia];
	const int* sdFunc = &(sd_func[ii * nprjMax]);
	double* prj = &(prjAll[ii * nprjMax*2*npack]);

        group.parallel_for_work_item([&](cl::sycl::h_item<1> item) {
                size_t i = item.get_global_id();
		size_t i2 = i * 2;

                int i_id = i_idx_sycl[i];
                int j_id = j_idx_sycl[i];
                int k_id = k_idx_sycl[i];

                int exi_id = 2 * i_id;
                int exj_id = 2 * j_id;
                int exk_id = 2 * k_id;

                double ai = exi[exi_id];
                double bi = exi[exi_id+1];
                double aj = exj[exj_id];
                double bj = exj[exj_id+1];
                double ak = exk[exk_id];
                double bk = exk[exk_id+1];

                double c  = aj*ak - bj*bk;
                double d  = aj*bk + ak*bj;

                double r_strfac = ai*c - bi*d;
                double i_strfac = ai*d + bi*c;

                for( int l=0; l<n_proj; l++) {
		    size_t lnpack2 = 2*l*npack;
                    double vnl_val = vnl_ion[i + l*npack];
                    if(sdFunc[l]) {
                        // mypneb->tcc_Mul_sycl()
                        prj[i2   + lnpack2] = vnl_val * r_strfac;
                        prj[i2+1 + lnpack2] = vnl_val * i_strfac;
                    } else {
                        // mypneb->tcc_iMul_sycl()
                        prj[i2   + lnpack2] = vnl_val * -i_strfac;
                        prj[i2+1 + lnpack2] = vnl_val *  r_strfac;
                    }
                }
            });
    });
        });

    // event.wait();
    // auto submit_time = event.get_profiling_info<cl::sycl::info::event_profiling::command_submit>();
    // auto start_time = event.get_profiling_info<cl::sycl::info::event_profiling::command_start>();
    // auto end_time = event.get_profiling_info<cl::sycl::info::event_profiling::command_end>();
    // auto submission_time = (start_time - submit_time) / 1000000.0f;
    // auto execution_time = (end_time - start_time) / 1000000.0f;
    // std::cout << "gen_projectors: " << execution_time << ", " << submission_time << std::endl;
}
*/
#endif
