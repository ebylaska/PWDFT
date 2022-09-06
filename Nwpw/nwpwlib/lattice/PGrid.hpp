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

namespace pwdft {


class PGrid : public d3db {

   Balance *mybalance;
   int balanced;
   double *Garray,*Gpack[2];
   double Gmax,Gmin;
   int     *masker[2],*packarray[2];
   int     nwave[2],nwave_entire[2],nwave_all[2],nida[2],nidb[2],nidb2[2];


   /* pfft data */
   int *zero_row2[2], *zero_row3[2], *zero_slab23[2];

   /* pfft_queue data */
   int aqmax,aqsize,alast_index;
   int    *aqindx,*aqstatus;
   double *atmp;

   /* zplane data */
   double *zplane_tmp1,*zplane_tmp2;


public:
        /* lattice pointer */
        Lattice *lattice;

        /* r_grid data */
        bool has_r_grid = false;
        double *r_grid;


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
            delete [] zplane_tmp1;
            delete [] zplane_tmp2;
            if (has_r_grid) delete [] r_grid;
        }

        double *Gxyz(const int i) { return &Garray[i*nfft3d]; }
        double *Gpackxyz(const int nb,const int i) { return &(Gpack[nb][i*(nida[nb]+nidb[nb])]); }
        double Gmax_ray() { return Gmax;}
        double Gmin_ray() { return Gmin;}
        double dGmin_ray() { return 0.01*Gmin;}
        int n_ray() {
           int nray0 = (int) ceil(100*(Gmax/Gmin) + 1.0);
           if (nray0<10) nray0 = 10;
           return nray0;
        }
        double *generate_G_ray() {
           int nray0 = (int) ceil(100*(Gmax/Gmin) + 1.0);
           if (nray0<10) nray0 = 10;
           double *g_ray = new (std::nothrow) double[nray0]();
           double dGmin = 0.01*Gmin;
           for (auto i=0; i<nray0; ++i)
              g_ray[i] = dGmin*i;
           return g_ray;
        }

        int nzero(const int nb) {return nida[nb];}
        int npack(const int nb)     {return (nida[nb] + nidb[nb]);}
        int npack_all(const int nb) {return nwave_all[nb];}
        int isbalanced() { return balanced;  }

        double *c_pack_allocate(const int nb) {
           double *ptr;
           ptr = new (std::nothrow) double [2*npack(nb)]();
           return ptr;
        }
        void c_pack_deallocate(double *ptr) { delete [] ptr;}


        void c_unpack(const int, double *);
        void c_pack(const int, double *);
        void cc_pack_copy(const int, double *, double *);
        double cc_pack_dot(const int, double *, double *);
        double cc_pack_idot(const int, double *, double *);
        void cc_pack_indot(const int, const int, double *, double *, double *);
        double tt_pack_dot(const int, double *, double *);
        double tt_pack_idot(const int, double *, double *);

        void cc_pack_inprjdot(const int, int, int, double *, double *, double *);

        void t_unpack(const int, double *);
        void t_pack(const int, double *);
        void tt_pack_copy(const int, double *, double *);
        void t_pack_nzero(const int, const int, double *);

        void tc_pack_copy( const int, double *, double *);

        void tcc_pack_Mul( const int, const double *, const double *, double *);
        void tc_pack_Mul(const int, const double *, double *);
        void tcc_pack_iMul(const int, const double *, const double *, double *);
        void tc_pack_iMul(const int, const double *, double *);

        void tcc_pack_MulSum2(const int, const double *, const double *, double *);
        void cc_pack_Sum2(const int, const double *, double *);
        void cccc_pack_Sum(const int, const double *, const double *, const double *, double *);
        void tcc_pack_aMulAdd(const int, const double, const double *, const double *, double *);



        void c_pack_addzero(const int, const double, double *);

        void c_pack_zero(const int, double *);
        void c_pack_SMul(const int, const double, double *);
        void cc_pack_SMul(const int, const double, const double *, double *);
        void cc_pack_daxpy(const int, const double, const double *, double *);
        void cct_pack_iconjgMul(const int, const double *, const double *, double *);
        void cct_pack_iconjgMulb(const int, const double *, const double *, double *);

        void i_pack(const int, int *);
        void ii_pack_copy(const int, int *, int *);

        void cr_pfft3b_queuein(const int, double *);
        void cr_pfft3b_queueout(const int, double *);
        int  cr_pfft3b_queuefilled();
        //int  cr_pfft3b(const int, double *);
        //int  rc_pfft3f(const int, double *);

};

}

#endif
