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
extern double Predictor_In(const int, const double [], const double []);
extern double Predictor_Out(const int, const double [], const double []);

extern double Corrector_In(const int, const double [], const double []);
extern double Corrector_In_F(const int, const double []);
extern double Corrector_Out(const int, const double [], const double []);

extern double Derivative5_1(const int, const double []);
extern double Derivative5_2(const int, const double []);
extern double Derivative5_3(const int, const double []);
extern double Derivative5_4(const int, const double []);
extern double Derivative5_5(const int, const double []);

extern double Derivative7_1(const int, const double []);
extern double Derivative7_2(const int, const double []);
extern double Derivative7_3(const int, const double []);
extern double Derivative7_4(const int, const double []);
extern double Derivative7_5(const int, const double []);
extern double Derivative7_6(const int, const double []);
extern double Derivative7_7(const int, const double []);
extern double Laplacian7_4(const int, const double []);

#endif
