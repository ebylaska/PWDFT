#ifndef _PAW_MY_MEMORY_H_
#define _PAW_MY_MEMORY_H_
/*
   $Id$
*/


extern float   *paw_alloc_1d_array(int);
extern float   **paw_alloc_2d_array(int,int);
extern void     paw_dealloc_1d_array(float*);
extern void     paw_dealloc_2d_array(int,int,float**);


#endif

