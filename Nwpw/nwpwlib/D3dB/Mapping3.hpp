#ifndef _MAPPING3_HPP_
#define _MAPPING3_HPP_

#pragma once

/* Mapping3.hpp
   Author - Eric Bylaska

        this class is used defining 3d parallel maps
*/

/**
 * @file Mapping3.hpp
 * @brief Defines the Mapping3 class for 3D parallel maps.
 * @author Eric Bylaska
 */


namespace pwdft {

/**
 * @class Mapping3
 * @brief Class for defining 3D parallel maps.
 */
class Mapping3 {

   // int *qmap[3],*pmap[3],*kmap;
   int *qmap[6], *pmap[6], *kmap;

public:
   int np, taskid;
   int maptype, n2ft3d, nfft3d, n2ft3d_map, nfft3d_map, nrft3d, nrft3d_map;
   int nx, ny, nz;
   int nq, nq1, nq2, nq3, nqr, nqr1, nqr2, nqr3;
 
   /* Constructors */
   /**
    * @brief Default constructor for Mapping3.
    */
   Mapping3();
 
   /**
    * @brief Constructor for Mapping3.
    *
    * @param np        Total number of processors.
    * @param taskid    Processor ID.
    * @param inmaptype Mapping type.
    * @param inne      Array containing additional parameters.
    * @param innx      Size of the x-dimension.
    * @param inny      Size of the y-dimension.
    * @param innz      Size of the z-dimension.
    */
   Mapping3(const int, const int, const int, const int, const int, const int);
 
   /* copy constructor */
   // Mapping3(const Mapping3&);
 
   /* destructor */
   /**
    * @brief Destructor for Mapping3.
    */
   ~Mapping3();
 
   //complex mapping
   /**
    * @brief Calculate the index for a given (i, j, k) position in the grid.
    *
    * Calculates the index in the grid for a specified (i, j, k) position taking into
    * account the mapping type. If the mapping type is 1, it uses qmap[0] for calculations;
    * otherwise, it uses qmap[2]. This function is used for complex mapping.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return The calculated index in the grid.
    */
   int ijktoindex(const int i, const int j, const int k) 
   {
      if (maptype == 1)
         return (i + j * (nx/2+1) + qmap[0][k]*(nx/2+1)*ny);
      else
         return (k + qmap[2][i + j * (nx / 2 + 1)] * nz);
   }
 
   /**
    * @brief Calculate the index for a given (i, j, k) position in the grid.
    *
    * Calculates the index in the grid for a specified (i, j, k) position taking into
    * account the mapping type. If the mapping type is 1, it uses qmap[0] and qmap[1] for calculations;
    * otherwise, it uses qmap[1]. This function is used for complex mapping.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return The calculated index in the grid.
    */
   int ijktoindex1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx / 2 + 1) + qmap[0][j] * (nx / 2 + 1) * nz);
     else
       return (j + qmap[1][k + i * nz] * ny);
   }

   /**
    * @brief Calculate the index for a given (i, j, k) position in the grid.
    *
    * Calculates the index in the grid for a specified (i, j, k) position taking into
    * account the mapping type. If the mapping type is 1, it uses qmap[0] and qmap[1] for calculations;
    * otherwise, it uses qmap[0]. This function is used for complex mapping.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return The calculated index in the grid.
    */
   int ijktoindex2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx + 2) + qmap[0][j] * (nx + 2) * ny);
     else
       return (i + qmap[0][j + k * ny] * (nx + 2));
   }

   /**
    * @brief Calculate the index for a given (i, j, k) position in the grid (transposed version).
    *
    * Calculates the index in the grid for a specified (i, j, k) position taking into account
    * the mapping type. If the mapping type is 1, it uses qmap[0] and qmap[1] for calculations;
    * otherwise, it uses qmap[0]. This function is a transposed version of ijktoindex2.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return The calculated index in the grid.
    */
   int ijktoindex2t(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx + 2) + qmap[0][j] * (nx / 2 + 1) * ny);
     else
       return (i + qmap[0][j + k * ny] * (nx / 2 + 1));
   }
 
   /**
    * @brief Get the pmap value for a given (i, j, k) position in the grid.
    *
    * Retrieves the pmap value for a specified (i, j, k) position in the grid based on
    * the mapping type. If the mapping type is 1, it uses pmap[0] for calculations;
    * otherwise, it uses pmap[2]. This function is used to determine the pmap value
    * for a given grid position.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return The pmap value corresponding to the given position.
    */
   int ijktop(const int i, const int j, const int k) 
   {
      if (maptype == 1)
         return (pmap[0][k]);
      else
         return (pmap[2][i + j * (nx/2+1)]);
   }

   /**
    * @brief Get the pmap value for a given (i, j, k) position in the grid.
    *
    * Retrieves the pmap value for a specified (i, j, k) position in the grid based on
    * the mapping type. If the mapping type is 1, it uses pmap[0] for calculations;
    * otherwise, it uses pmap[1]. This function is used to determine the pmap value
    * for a given grid position.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return The pmap value corresponding to the given position.
    */
   int ijktop1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (pmap[0][j]);
     else
       return (pmap[1][k + i * nz]);
   }

   /**
    * @brief Get the pmap value for a given (i, j, k) position in the grid.
    *
    * Retrieves the pmap value for a specified (i, j, k) position in the grid based on
    * the mapping type. If the mapping type is 1, it uses pmap[0] for calculations;
    * otherwise, it uses pmap[0]. This function is used to determine the pmap value
    * for a given grid position.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return The pmap value corresponding to the given position.
    */
   int ijktop2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (pmap[0][j]);
     else
       return (pmap[0][j + k * ny]);
   }
 
   /**
    * @brief Get the qmap value for a given (i, j, k) position in the grid.
    *
    * Retrieves the qmap value for a specified (i, j, k) position in the grid based on
    * the mapping type. If the mapping type is 1, it uses qmap[0] for calculations;
    * otherwise, it uses qmap[2]. This function is used to determine the qmap value
    * for a given grid position.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return The qmap value corresponding to the given position.
    */
   int ijktoq(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][k]);
     else
       return (qmap[2][i + j * (nx / 2 + 1)]);
   }

   /**
    * @brief Get the qmap value (alternative) for a given (i, j, k) position in the grid.
    *
    * Retrieves an alternative qmap value for a specified (i, j, k) position in the grid
    * based on the mapping type. If the mapping type is 1, it uses qmap[0][j] for calculations;
    * otherwise, it uses qmap[1][k + i * nz]. This function is an alternative way to determine
    * the qmap value for a given grid position.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return An alternative qmap value corresponding to the given position.
    */
   int ijktoq1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][j]);
     else
       return (qmap[1][k + i * nz]);
   }

   /**
    * @brief Get the qmap value (alternative) for a given (i, j, k) position in the grid.
    *
    * Retrieves an alternative qmap value for a specified (i, j, k) position in the grid
    * based on the mapping type. If the mapping type is 1, it uses qmap[0][j] for calculations;
    * otherwise, it uses qmap[0][j + k * ny]. This function is an alternative way to determine
    * the qmap value for a given grid position.
    *
    * @param i The x-coordinate of the position.
    * @param j The y-coordinate of the position.
    * @param k The z-coordinate of the position.
    * @return An alternative qmap value corresponding to the given position.
    */
   int ijktoq2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][j]);
     else
       return (qmap[0][j + k * ny]);
   }
 
   // real mapping
   int ijkrtoindex(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + j * (nx + 2) + qmap[0][k] * (nx + 2) * ny);
     else
       return (k + qmap[5][i + j * (nx + 2)] * nz);
   }
   int ijkrtoindex1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx + 2) + qmap[0][j] * (nx + 2) * nz);
     else
       return (j + qmap[4][k + i * nz] * ny);
   }
   int ijkrtoindex2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx + 2) + qmap[0][j] * (nx + 2) * ny);
     else
       return (i + qmap[3][j + k * ny] * (nx + 2));
   }
   int ijkrtoindex2t(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx + 2) + qmap[0][j] * (nx + 2) * ny);
     else
       return (i + qmap[3][j + k * ny] * (nx + 2));
   }
 
   int ijkrtop(const int i, const int j, const int k) {
     if (maptype == 1)
       return (pmap[0][k]);
     else
       return (pmap[5][i + j * (nx + 2)]);
   }
   int ijkrtop1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (pmap[0][j]);
     else
       return (pmap[4][k + i * nz]);
   }
   int ijkrtop2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (pmap[0][j]);
     else
       return (pmap[3][j + k * ny]);
   }
 
   int ijkrtoq(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][k]);
     else
       return (qmap[5][i + j * (nx + 2)]);
   }
   int ijkrtoq1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][j]);
     else
       return (qmap[4][k + i * nz]);
   }
   int ijkrtoq2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][j]);
     else
       return (qmap[3][j + k * ny]);
   }
};

} // namespace pwdft

#endif
