/* k1db.cpp
   Author - Eric Bylaska

        this class is used for defining k point parallel maps
*/

/**
 * @class k1db
 * @brief Container for operations on distributed k point blocks.
 *
 * The `k1db` class is designed to handle operations related to distributed 3D blocks.
 * It provides functionality for performing matrix multiplications and other mathematical
 * operations on distributed data structures. This class is often used in conjunction
 * with the `Parallel` and `gdevice2` classes to manage parallelism and device operations.
 */

/*
#include        <cmath>
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>

*/

#include <iostream>
#include "Mapping1k.hpp"
#include "Parallel.hpp"
#include "blas.h"
#include "gdevice2.hpp"
#include "k1db.hpp"

#include "iofmt.hpp"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

namespace pwdft {


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/**
 * @brief Constructor for the d1db class.
 *
 * This constructor initializes an instance of the d1db class, which is used for
 * distributed 3D block operations. It sets up the necessary data structures and
 * communication channels for distributed matrix operations.
 *
 * @param inparall A pointer to the Parallel object managing parallelism.
 * @param inmaptype The mapping type for distributed operations.
 * @param inispin The initial spin value.
 * @param inne An array specifying the dimensions of the distributed block.
 *
 * @note This constructor also initiates asynchronous communication channels if
 *       the number of processes in the j-dimension (np_j) is greater than 1.
 *       These channels are used in distributed matrix operations.
 */
k1db::k1db(Parallel *inparall, const int inmaptype, const int innbrill)
     : Mapping1k(inmaptype, inparall->np_k(), inparall->taskid_k(), innbrill) 
{
  parall = inparall;

  int np_k = parall->np_k();
  if (np_k>1)
  {
     //request1_indx = parall->max_reqstat-2;
     //request2_indx = parall->max_reqstat-1;
     //parall->astart(request1_indx,2*np_k+2);
     //parall->astart(request2_indx,2*np_k+2);
  }
}



} // namespace pwdft

