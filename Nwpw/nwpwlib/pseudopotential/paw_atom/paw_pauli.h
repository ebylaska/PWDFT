#ifndef _PAW_PAULI_H_
#define _PAW_PAULI_H_

/*
   $Id$
*/


extern int   paw_R_Pauli(int n,
                                   int l,
                                   float  Z,
                                   float* v,
                                   float *Eig,
                                   float* u,
                                   float* uprime);

extern  int paw_R_Pauli_Fixed_E(
        int n,
        int l,
        float Z,
        float *v,
        int match,
        float E,
        float *u,
        float *uprime
    );



#endif



