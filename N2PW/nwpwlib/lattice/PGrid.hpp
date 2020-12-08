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
   float *Garray;
   int     *masker[2],*packarray[2];
   int     nwave[2],nwave_entire[2],nwave_all[2],nida[2],nidb[2],nidb2[2];

   /* pfft data */
   int *zero_row2[2], *zero_row3[2], *zero_slab23[2];

   /* pfft_queue data */
   int aqmax,aqsize,alast_index;
   int    *aqindx,*aqstatus;
   float *atmp;

public:

        /* constructor */
	PGrid(Parallel *);

        /* destructor */
        ~PGrid() { 
            delete [] Garray; 
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

        float *Gxyz(const int i) { return &Garray[i*nfft3d]; }

        int nzero(const int i) {return nida[i];}
        int npack(const int i)     {return (nida[i] + nidb[i]);}
        int npack_all(const int i) {return nwave_all[i];}
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

        void cr_pfft3b_queuein(const int, float *);
        void cr_pfft3b_queueout(const int, float *);
        int  cr_pfft3b_queuefilled();

        void t_pack(const int, float *);
        void tt_pack_copy(const int, float *, float *);
        void tcc_Mul( const int, float *, float *, float *);
        void tcc_iMul(const int, float *, float *, float *);
        void tcc_MulSum2( const int, float *, float *, float *);
        void cc_Sum2(const int, float *, float *);
        void c_zero(const int, float *);
        void c_SMul(const int, float, float *);
        void cc_SMul(const int, float, float *, float *);
        void cc_daxpy(const int, float, float *, float *);

        void i_pack(const int, int *);
        void ii_pack_copy(const int, int *, int *);

};

#endif
