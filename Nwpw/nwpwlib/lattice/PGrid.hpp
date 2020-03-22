#ifndef _PGRID_H_
#define _PGRID_H_
/* d3db.h
   Author - Eric Bylaska

*/

#include	"Parallel.hpp"
#include	"d3db.hpp"
#include	"lattice.hpp"
#include	"Balance.hpp"

class PGrid : public d3db {

   Balance *mybalance;
   int balanced;
   double *Garray,*Gpack[2];
   int     *masker[2],*packarray[2];
   int     nwave[2],nwave_entire[2],nwave_all[2],nida[2],nidb[2],nidb2[2];

   /* pfft data */
   int *zero_row2[2], *zero_row3[2], *zero_slab23[2];

   /* pfft_queue data */
   int aqmax,aqsize,alast_index;
   int    *aqindx,*aqstatus;
   double *atmp;

public:

        /* constructor */
	PGrid(Parallel *);

        /* destructor */
        ~PGrid() { 
            delete [] Garray; 
            delete [] Gpack[0]; 
            delete [] masker[0];
            delete [] packarray[0];
            if (balanced) delete mybalance;
            delete [] zero_row3[0];
            delete [] zero_row3[1];
            delete [] zero_row2[0];
            delete [] zero_row2[1];
            delete [] zero_slab23[0];
            delete [] zero_slab23[1];
        }

        double *Gxyz(const int i) { return &Garray[i*nfft3d]; }
        double *Gpackxyz(const int nb,const int i) { return &(Gpack[nb][i*(nida[nb]+nidb[nb])]); }

        int nzero(const int nb) {return nida[nb];}
        int npack(const int nb)     {return (nida[nb] + nidb[nb]);}
        int npack_all(const int nb) {return nwave_all[nb];}
        int isbalanced() { return balanced;  }

        double *c_pack_allocate(const int nb) {
           double *ptr;
           ptr = new double [2*npack(nb)];
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

        void cr_pfft3b_queuein(const int, double *);
        void cr_pfft3b_queueout(const int, double *);
        int  cr_pfft3b_queuefilled();

        void t_pack(const int, double *);
        void tt_pack_copy(const int, double *, double *);
        void tcc_Mul( const int, double *, double *, double *);
        void tcc_iMul(const int, double *, double *, double *);
        void tcc_MulSum2( const int, double *, double *, double *);
        void cc_Sum2(const int, double *, double *);
        void c_zero(const int, double *);
        void c_SMul(const int, double, double *);
        void cc_SMul(const int, double, double *, double *);
        void cc_daxpy(const int, double, double *, double *);
        void cct_iconjgMul(const int, const double *, const double *, double *);
        void cct_iconjgMulb(const int, const double *, const double *, double *);

        void i_pack(const int, int *);
        void ii_pack_copy(const int, int *, int *);

};

#endif
