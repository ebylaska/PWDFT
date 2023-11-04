#ifndef _MAPPING3c_HPP_
#define _MAPPING3c_HPP_

#pragma once

/* Mapping3c.hpp
   Author - Eric Bylaska

        this class is used defining 3d parallel maps for complex FFTs.
*/

/**
 * @file Mapping3c.hpp
 * @brief Defines the Mapping3c class for 3D parallel maps.
 * @author Eric Bylaska
 */


namespace pwdft {

/**
 * @class Mapping3c
 * @brief Class for defining 3D parallel maps.
 */
class Mapping3c {

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
   Mapping3c();
 
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
   Mapping3c(const int, const int, const int, const int, const int, const int);
 
   /* copy constructor */
   // Mapping3(const Mapping3&);
 
   /* destructor */
   /**
    * @brief Destructor for Mapping3.
    */
   ~Mapping3c();
 
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
   int cijktoindex(const int i, const int j, const int k) 
   {
      if (maptype == 1)
         return (i + j * nx + qmap[0][k]*nx*ny);
      else
         return (k + qmap[2][i + j*nx] * nz);
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
   int cijktoindex1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx) + qmap[0][j] * (nx) * nz);
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
   int cijktoindex2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx) + qmap[0][j] * (nx) * ny);
     else
       return (i + qmap[0][j + k * ny] * (nx));
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
   int cijktoindex2t(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx) + qmap[0][j] * (nx) * ny);
     else
       return (i + qmap[0][j + k * ny] * (nx));
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
   int cijktop(const int i, const int j, const int k) 
   {
      if (maptype == 1)
         return (pmap[0][k]);
      else
         return (pmap[2][i + j * (nx)]);
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
   int cijktop1(const int i, const int j, const int k) {
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
   int cijktop2(const int i, const int j, const int k) {
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
   int cijktoq(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][k]);
     else
       return (qmap[2][i + j * (nx)]);
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
   int cijktoq1(const int i, const int j, const int k) {
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
   int cijktoq2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][j]);
     else
       return (qmap[0][j + k * ny]);
   }
 
   // real mapping
   int cijkrtoindex(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + j * (nx) + qmap[0][k] * (nx) * ny);
     else
       return (k + qmap[5][i + j * (nx)] * nz);
   }
   int cijkrtoindex1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx) + qmap[0][j] * (nx) * nz);
     else
       return (j + qmap[4][k + i * nz] * ny);
   }
   int cijkrtoindex2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx) + qmap[0][j] * (nx) * ny);
     else
       return (i + qmap[3][j + k * ny] * (nx));
   }
   int cijkrtoindex2t(const int i, const int j, const int k) {
     if (maptype == 1)
       return (i + k * (nx) + qmap[0][j] * (nx) * ny);
     else
       return (i + qmap[3][j + k * ny] * (nx));
   }
 
   int cijkrtop(const int i, const int j, const int k) {
     if (maptype == 1)
       return (pmap[0][k]);
     else
       return (pmap[5][i + j * (nx)]);
   }
   int cijkrtop1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (pmap[0][j]);
     else
       return (pmap[4][k + i * nz]);
   }
   int cijkrtop2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (pmap[0][j]);
     else
       return (pmap[3][j + k * ny]);
   }
 
   int cijkrtoq(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][k]);
     else
       return (qmap[5][i + j * (nx)]);
   }
   int cijkrtoq1(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][j]);
     else
       return (qmap[4][k + i * nz]);
   }
   int cijkrtoq2(const int i, const int j, const int k) {
     if (maptype == 1)
       return (qmap[0][j]);
     else
       return (qmap[3][j + k * ny]);
   }
};

} // namespace pwdft

#endif
