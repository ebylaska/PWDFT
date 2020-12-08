/*
 $Id$
*/
#ifndef _PRED_CORR_H_
#define _PRED_CORR_H_

/* Pred_Corr.h - 6/9/95
   author      - Eric Bylaska

   This file contains the 5th order predictor-corrector
formulas for integrating inward and outward.
This file also contains 5th order derivatives.

*/
extern	float  Predictor_In();
extern	float	Predictor_Out();

extern	float	Corrector_In();
extern  float Corrector_In_F();
extern	float	Corrector_Out();

extern	float	Derivative5_1();
extern	float	Derivative5_2();
extern	float	Derivative5_3();
extern	float	Derivative5_4();
extern	float	Derivative5_5();

extern	float	Derivative7_1();
extern	float	Derivative7_2();
extern	float	Derivative7_3();
extern	float	Derivative7_4();
extern	float	Derivative7_5();
extern	float	Derivative7_6();
extern	float	Derivative7_7();
extern	float	Laplacian7_4();

#endif
