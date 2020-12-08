#ifndef _PGRID_H_
#define _PGRID_H_
/* d3db.h
   Author - Eric Bylaska

*/

#pragma once

#include	<cmath>
#include	"Parallel.hpp"
#include	"d3db.hpp"
#include	"Lattice.hpp"
#include	"Balance.hpp"

class PGrid : public d3db {

   Balance *mybalance;
   int balanced;
   float *Garray,*Gpack[2];
   float Gmax,Gmin;
   int     *masker[2],*packarray[2];
   int     nwave[2],nwave_entire[2],nwave_all[2],nida[2],nidb[2],nidb2[2];


   /* pfft data */
   int *zero_row2[2], *zero_row3[2], *zero_slab23[2];

   /* pfft_queue data */
   int aqmax,aqsize,alast_index;
   int    *aqindx,*aqstatus;
   float *atmp;


public:
        Lattice *lattice;

        /* constructor */
	PGrid(Parallel *, Lattice *, int, int, int, int, int);
	PGrid(Parallel *, Lattice *, Control2&);

        /* destructor */
        ~PGrid() {
            delete [] Garray;
            delete [] Gpack[0];
            //delete [] Gpack[1];
            delete [] masker[0];
            //delete [] masker[1];
            delete [] packarray[0];
            //delete [] packarray[1];
            if (balanced) delete mybalance;
            delete [] zero_row3[0];
            delete [] zero_row3[1];
            delete [] zero_row2[0];
            delete [] zero_row2[1];
            delete [] zero_slab23[0];
            delete [] zero_slab23[1];
        }

        float *Gxyz(const int i) { return &Garray[i*nfft3d]; }
        float *Gpackxyz(const int nb,const int i) { return &(Gpack[nb][i*(nida[nb]+nidb[nb])]); }
        float Gmax_ray() { return Gmax;}
        float Gmin_ray() { return Gmin;}
        float dGmin_ray() { return 0.01*Gmin;}
        int n_ray() {
           int nray0 = (int) ceil(100*(Gmax/Gmin) + 1.0);
           if (nray0<10) nray0 = 10;
           return nray0;
        }
        float *generate_G_ray() {
           int nray0 = (int) ceil(100*(Gmax/Gmin) + 1.0);
           if (nray0<10) nray0 = 10;
           float *g_ray = new float[nray0];
           float dGmin = 0.01*Gmin;
           for (auto i=0; i<nray0; ++i)
              g_ray[i] = dGmin*i;
           return g_ray;
        }

        int nzero(const int nb) {return nida[nb];}
        int npack(const int nb)     {return (nida[nb] + nidb[nb]);}
        int npack_all(const int nb) {return nwave_all[nb];}
        int isbalanced() { return balanced;  }

        float *c_pack_allocate(const int nb) {
           float *ptr;
           ptr = new float [2*npack(nb)];
           return ptr;
        }
        void c_pack_deallocate(float *ptr) { delete [] ptr;}


        void c_unpack(const int, float *);
        void c_pack(const int, float *);
        void cc_pack_copy(const int, float *, float *);
        float cc_pack_dot(const int, float *, float *);
        float cc_pack_idot(const int, float *, float *);
        void cc_pack_indot(const int, const int, float *, float *, float *);
        float tt_pack_dot(const int, float *, float *);
        float tt_pack_idot(const int, float *, float *);

        void cc_pack_inprjdot(const int, int, int, float *, float *, float *);

        void t_unpack(const int, float *);
        void t_pack(const int, float *);
        void tt_pack_copy(const int, float *, float *);
        void t_pack_nzero(const int, const int, float *);

        void tc_pack_copy( const int, float *, float *);

        void tcc_Mul( const int, float *, float *, float *);
        void tcc_iMul(const int, float *, float *, float *);
        void tcc_MulSum2( const int, float *, float *, float *);
        void cc_Sum2(const int, float *, float *);
        void c_zero(const int, float *);
        void c_SMul(const int, float, float *);
        void cc_SMul(const int, float, float *, float *);
        void cc_daxpy(const int, float, float *, float *);
        void cct_iconjgMul(const int, const float *, const float *, float *);
        void cct_iconjgMulb(const int, const float *, const float *, float *);

        void i_pack(const int, int *);
        void ii_pack_copy(const int, int *, int *);

        void cr_pfft3b_queuein(const int, float *);
        void cr_pfft3b_queueout(const int, float *);
        int  cr_pfft3b_queuefilled();
        //int  cr_pfft3b(const int, float *);
        //int  rc_pfft3f(const int, float *);

#ifdef NWPW_SYCL
        void cc_pack_inprjdot_sycl(const int, int, int, float *, float *, float *);
        void tcc_Mul_sycl( const int, const float *, const float *, float *);
        void tcc_iMul_sycl(const int, const float *, const float *, float *);
#endif

};

#endif
