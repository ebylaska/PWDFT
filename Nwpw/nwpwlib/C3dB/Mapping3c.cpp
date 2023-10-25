/* Mapping3c.cpp
   Author - Eric Bylaska

        this class is used for defining 3d parallel maps
*/

/**
 * @class Mapping3c
 * @brief Class for defining 3D parallel maps.
 */

/*
#include        <cmath>
#include        <cstdlib>
#include        <iostream>

*/

#include "Mapping3c.hpp"
#include "hilbert.hpp"
#include <iostream>

namespace pwdft {

/***********************************
 *                                 *
 *       c_generate_map_indexes    *
 *                                 *
 ***********************************/
/**
 * @brief Generate mapping indexes for parallel tasks.
 *
 * This function generates mapping indexes for parallel tasks based on the specified
 * parameters. It assigns processor IDs and task IDs (P-map and Q-map) to grid
 * points in a distributed manner, allowing for load balancing among processes.
 *
 * @param taskid The ID of the current task.
 * @param np The total number of parallel processes.
 * @param ny The size of the Y dimension of the grid.
 * @param nz The size of the Z dimension of the grid.
 * @param pmap An array to store the processor IDs (P-map) for each grid point.
 * @param qmap An array to store the task IDs (Q-map) for each grid point.
 *
 * @return The number of grid points assigned to the current task (taskid).
 */
static int c_generate_map_indexes(int taskid, int np, int ny, int nz, int *pmap, int *qmap) 
{
   int nq, q, p;
   int *indx_proc;
   int *indx_q;
   int nq1 = (ny * nz) / np;
   int rmdr1 = (ny * nz) % np;
   int nq2 = nq1;
 
   indx_proc = new (std::nothrow) int[ny * nz]();
   indx_q = new (std::nothrow) int[ny * nz]();
 
   if (rmdr1 > 0)
     ++nq2;
   nq = 0;
   p = 0;
   q = 0;
   for (int i = 0; i < (ny * nz); ++i) {
     indx_proc[i] = p;
     indx_q[i] = q;
     if (taskid == p)
       ++nq;
     ++q;
     if (q >= nq2) {
       q = 0;
       ++p;
       p = p % np;
       if (p >= rmdr1)
         nq2 = nq1;
     }
   }
 
   for (int k = 0; k < nz; ++k)
     for (int j = 0; j < ny; ++j) {
       p = indx_proc[pmap[j + k * ny]];
       q = indx_q[pmap[j + k * ny]];
       pmap[j + k * ny] = p;
       qmap[j + k * ny] = q;
     }
   delete[] indx_proc;
   delete[] indx_q;
 
   return nq;
}


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/**
 * @brief Destructor for the Mapping3 class.
 *
 * This destructor releases the dynamically allocated memory for various member variables.
 */
Mapping3c::Mapping3c() 
{
   maptype = 0;
   cn2ft3d = 0;
   cnfft3d = 0;
   nx = 0;
   ny = 0;
   nz = 0;
   qmap[0] = NULL;
   qmap[1] = NULL;
   qmap[2] = NULL;
   qmap[3] = NULL;
   qmap[4] = NULL;
   qmap[5] = NULL;
   pmap[0] = NULL;
   pmap[1] = NULL;
   pmap[2] = NULL;
   pmap[3] = NULL;
   pmap[4] = NULL;
   pmap[5] = NULL;
   kmap = NULL;
}

/**
 * @brief Constructor for the Mapping3 class with parameters.
 *
 * This constructor initializes the `Mapping3` object with the provided parameters.
 *
 * @param mapin The mapping type.
 * @param npin The number of parallel processes.
 * @param taskidin The ID of the current task.
 * @param nxin The size of the X dimension.
 * @param nyin The size of the Y dimension.
 * @param nzin The size of the Z dimension.
 */
Mapping3c::Mapping3c(const int mapin, const int npin, const int taskidin,
                     const int nxin,  const int nyin, const int nzin) 
{
   int p, q;
 
   qmap[0] = NULL;
   qmap[1] = NULL;
   qmap[2] = NULL;
   qmap[3] = NULL;
   qmap[4] = NULL;
   qmap[5] = NULL;
   pmap[0] = NULL;
   pmap[1] = NULL;
   pmap[2] = NULL;
   pmap[3] = NULL;
   pmap[4] = NULL;
   pmap[5] = NULL;
   kmap = NULL;
 
   np = npin;
   taskid = taskidin;
   nx = nxin;
   ny = nyin;
   nz = nzin;
   maptype = mapin;
 
   /* slab mapping */
   if ((maptype == 1) || (maptype == -1)) {
     maptype = 1;
     qmap[0] = new (std::nothrow) int[nz]();
     pmap[0] = new (std::nothrow) int[nz]();
     kmap = new (std::nothrow) int[nz]();
 
     /* cyclic mapping */
     p = 0;
     q = 0;
     for (auto k = 0; k < nz; ++k) 
     {
        qmap[0][k] = q;
        pmap[0][k] = p;
        if (p == taskid)
           nq = q + 1;
        ++p;
        if (p >= np) {
           p = 0;
           ++q;
        }
     }
     cnfft3d = (nx) * ny * nq;
     cn2ft3d = cnfft3d;
 
     cnfft3d_map = cnfft3d;
     cn2ft3d_map = cn2ft3d;
 
   }
 
   /* hilbert or hcurve  mapping */
   else 
   {
      qmap[0] = new (std::nothrow) int[ny * nz]();
      pmap[0] = new (std::nothrow) int[ny * nz]();
      qmap[1] = new (std::nothrow) int[nz * (nx)]();
      pmap[1] = new (std::nothrow) int[nz * (nx)]();
      qmap[2] = new (std::nothrow) int[(nx) * ny]();
      pmap[2] = new (std::nothrow) int[(nx) * ny]();
     
      qmap[3] = new (std::nothrow) int[ny * nz]();
      pmap[3] = new (std::nothrow) int[ny * nz]();
      qmap[4] = new (std::nothrow) int[nz * (nx)]();
      pmap[4] = new (std::nothrow) int[nz * (nx)]();
      qmap[5] = new (std::nothrow) int[(nx) * ny]();
      pmap[5] = new (std::nothrow) int[(nx) * ny]();
     
      if (maptype > 0) {
        if (maptype == 2) {
          hilbert2d_map(ny, nz, pmap[0]);
          hilbert2d_map(nz, nx, pmap[1]);
          hilbert2d_map(nx, ny, pmap[2]);
     
          hilbert2d_map(ny, nz, pmap[3]);
          hilbert2d_map(nz, nx, pmap[4]);
          hilbert2d_map(nx, ny, pmap[5]);
        }
        if (maptype == 3) {
          hcurve2d_map(ny, nz, pmap[0]);
          hcurve2d_map(nz, nx, pmap[1]);
          hcurve2d_map(nx, ny, pmap[2]);
     
          hcurve2d_map(ny, nz, pmap[3]);
          hcurve2d_map(nz, nx, pmap[4]);
          hcurve2d_map(nx, ny, pmap[5]);
        }
        nq1 = c_generate_map_indexes(taskid, np, ny, nz, pmap[0], qmap[0]);
        nq2 = c_generate_map_indexes(taskid, np, nz, nx, pmap[1], qmap[1]);
        nq3 = c_generate_map_indexes(taskid, np, nx, ny, pmap[2], qmap[2]);
     
        nqr1 = c_generate_map_indexes(taskid, np, ny, nz, pmap[3], qmap[3]);
        nqr2 = c_generate_map_indexes(taskid, np, nz, nx, pmap[4], qmap[4]);
        nqr3 = c_generate_map_indexes(taskid, np, nx, ny, pmap[5], qmap[5]);
     
        /* double grid map1 defined wrt to single grid         */
        /* makes expand and contract routines trivial parallel */
      } 
      
      cnfft3d = nx * nq1;
      if ((ny * nq2) > cnfft3d)
         cnfft3d = ny * nq2;
      if ((nz * nq3) > cnfft3d)
         cnfft3d = nz * nq3;
      cn2ft3d = cnfft3d;
      cnfft3d_map = nz * nq3;
      cn2ft3d_map = nx * nq1;
     
      nrft3d = nx * nqr1;
      if ((ny * nqr2) > cnfft3d)
         nrft3d = ny * nqr2;
      if ((nz * nqr3) > cnfft3d)
         nrft3d = nz * nqr3;
      nrft3d_map = nx * nqr1;
   }
}


/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
/**
 * @brief Destructor for the Mapping3 class.
 *
 * This destructor releases the dynamically allocated memory for various member variables.
 */
Mapping3c::~Mapping3c() 
{
   for (int i = 0; i < 6; ++i) 
   {
      if (qmap[i])
         delete[] pmap[i];
      if (pmap[i])
         delete[] qmap[i];
   }
   if (kmap)
      delete[] kmap;
 
   if ((maptype == 2) || (maptype == 3)) {
     // int *h_iq_to_i1[6],*h_iq_to_i2[6];
     // int *h_i1_start[6],*h_i2_start[6];
   }
}


} // namespace pwdft
