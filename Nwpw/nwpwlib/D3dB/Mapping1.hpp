#ifndef _MAPPING1d_HPP_
#define _MAPPING1d_HPP_

#pragma once

/* Mapping1.hpp
   Author - Eric Bylaska

        this class is used defining 1d parallel maps
*/

/**
 * @file Mapping1.hpp
 * @brief Defines the Mapping1 class for 1D parallel maps.
 * @author Eric Bylaska
 */


namespace pwdft {

/**
 * @class Mapping1
 * @brief Class for defining 1D parallel maps.
 */

class Mapping1 {

   int np, taskid;

   int *qmap[2], *pmap[2], *kmap[2], *nqarray[2];

public:
   int ispin, ne[2];
   int maptype1;
   int neq[2]; // Number of wave functions
 
   /* Constructors */
   /**
   * @brief Default constructor for Mapping1.
   */
   Mapping1();
 
   /**
    * @brief Constructor for Mapping1.
    *
    * @param ispin     Number of spins.
    * @param np        Total number of processors.
    * @param taskid    Processor ID.
    * @param inmaptype Mapping type.
    * @param inne      Array containing additional parameters.
    */
   Mapping1(const int, const int, const int, const int, const int *);
 
   /* destructor */
   /**
    * @brief Destructor for Mapping1.
    */
   ~Mapping1() {
     for (int ms = 0; ms < ispin; ++ms) {
       delete[] qmap[ms];
       delete[] pmap[ms];
       delete[] kmap[ms];
       delete[] nqarray[ms];
     }
   }
 
   /**
    * @brief Convert spin and index to an array index.
    *
    * @param ms Spin index.
    * @param n  Array index.
    * @return Array index corresponding to the given spin and index.
    */
   int msntoindex(const int ms, const int n) {
     return (qmap[ms][n] + ms * neq[0]);
   }
 
   /**
    * @brief Convert spin and index to a processor ID.
    *
    * @param ms Spin index.
    * @param n  Array index.
    * @return Processor ID corresponding to the given spin and index.
    */
   int msntop(const int ms, const int n) { return pmap[ms][n]; }
};

} // namespace pwdft

#endif
