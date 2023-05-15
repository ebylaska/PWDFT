/*
 $Id$
*/
#ifndef _GRID_H_
#define _GRID_H_
/* grid.h
   author - Eric Bylaska
*/

extern void init_Grids(int);
extern void end_Grids(void);
extern double *alloc_Grid(void);
extern void dealloc_Grid(double *);

#endif
