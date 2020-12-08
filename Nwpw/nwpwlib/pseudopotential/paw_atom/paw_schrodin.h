#ifndef _PAW_SCHRODINGER_H_
#define _PAW_SCHRODINGER_H_

/*
   $Id$
*/

/* Schrodinger.h - 6/9/95
   author     - Eric Bylaska

   This file contains routines for integrating the radial
   Schodinger equation.

*/

extern int   paw_R_Schrodinger(int n,
                                   int l,
                                   float* v,
                                   float *Eig,
                                   float* u,
                                   float* uprime);

extern  int paw_R_Schrodinger_Fixed_E(
        int l,
        float *v,
        int match,
        float E,
        float *u,
        float *uprime
    );
extern int paw_R_Schrodinger_Fixed_Logderiv(
        int n,
        int l,
        float *v,
        int match,
        float u_logderiv,
        float *Eig,
        float *u,
        float *uprime
    );


extern  void paw_R_Schrodinger_Fixed_E1(
        int l,
        float *v,
        float *f,
        int match,
        float E,
        float *u,
        float *uprime
    );

extern int paw_R_Schrodinger_Fixed_Logderiv1(
        int n,
        int l,
        float *v,
        int match,
        float u_logderiv,
        float *Eig,
        float *u,
        float *uprime
    );


#endif



