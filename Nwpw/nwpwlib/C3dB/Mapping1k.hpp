#ifndef _MAPPING1k_HPP_
#define _MAPPING1k_HPP_

#pragma once

/* Mapping1k.hpp
   Author - Eric Bylaska

        this class is used defining 1d parallel maps
*/

/**
 * @file Mapping1k.hpp
 * @brief Defines the Mapping1 class for 1D parallel maps.
 * @author Eric Bylaska
 */


namespace pwdft {

/**
 * @class Mapping1
 * @brief Class for defining 1D parallel maps.
 */

class Mapping1k {

   int np, taskid;

   int *qmap, *pmap, *kmap, *nqarray;

public:
   int nbrillouin; // number of k-points
   int nbrillq;    // number of k-points
   int maptype1;
 
   /* Constructors */
   /**
   * @brief Default constructor for Mapping1k.
   */
   Mapping1k();
 
   /**
    * @brief Constructor for Mapping1,.
    *
    * @param nbrill    Number of k-points.
    * @param np        Total number of processors.
    * @param taskid    Processor ID.
    * @param inmaptype Mapping type.
    */
   Mapping1k(const int, const int, const int, const int);
 
   /* destructor */
   /**
    * @brief Destructor for Mapping1.
    */
   ~Mapping1k() {
       delete[] qmap;
       delete[] pmap;
       delete[] kmap;
       delete[] nqarray;
   }
 
   /**
    * @brief Convert k index to an array index.
    *
    * @param k vector index.
    * @return Array index corresponding to the given k vector.
    */
   int ktoindex(const int k) {
     return (qmap[k]);
   }
 
   /**
    * @brief Convert k index to a processor ID.
    *
    * @param k vector index.
    * @return Processor ID corresponding to the given k vector.
    */
   int ktop(const int k) { return pmap[k]; }
};

} // namespace pwdft

#endif
