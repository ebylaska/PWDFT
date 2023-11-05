/* Mapping1k.cpp
   Author - Eric Bylaska

        this class is used for defining 1d parallel maps
*/

/**
 * @class Mapping1k
 * @brief Class for defining 1D parallel maps.
 */

/*
#include        <cmath>
#include        <cstdlib>
#include        <iostream>

*/

#include "Mapping1k.hpp"
#include <iostream>

namespace pwdft {

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

Mapping1k::Mapping1k() {
  maptype1 = 0;
  nbrillouin = 0;
  nbrillq = 0;

  qmap = (int *)0;
  pmap = (int *)0;
  kmap = (int *)0;
  nqarray = (int *)0;
}

Mapping1k::Mapping1k(const int mapin, const int npin, const int taskidin, const int nbrillin) 
{
   // int k,p,q;

   np = npin;
   taskid = taskidin;
   maptype1 = mapin;

   nbrillouin = 0;
   nbrillq    = 0;

   nbrillouine = nbrillin;
   qmap = new (std::nothrow) int[nbrillouin];
   pmap = new (std::nothrow) int[nbrillouin];
   kmap = new (std::nothrow) int[nbrillouin];
   nqarray = new (std::nothrow) int[np];

   /* cyclic mapping */
   if (maptype1 == 0) 
   {
      /* cyclic mapping */
      int p = 0;
      int q = 0;
      for (auto k = 0; k < nbrillouin; ++k) 
      {
         qmap[k] = q;
         pmap[k] = p;
        
         if (p == taskid)
            nbrillq = q;
         p++;
         if (p >= np) 
         {
            p = 0;
            ++q;
         }
      }
   }

   /* block mapping */
   else 
   {
      int k, q;
      int p;
      for (p = 0; p < np; ++p)
         nqarray[p] = 0;
      p = 0;
      for (k = 0; k < nbrillouin; ++k) 
      {
         nqarray[p] += 1;
         p = (p + 1) % np;
      }
     
      k = 0;
      for (p = 0; p < np; ++p)
         for (q = 0; q < nqarray[p]; ++q) 
         {
            qmap[k] = q;
            pmap[k] = p;
            ++k;
         }
      nbrillq = nqarray[taskid];
   }

   for (auto k = 0; k < nbrillouin; ++k)
      if (pmap[k] == taskid)
         kmap[qmap[k]] = k;
  
}

} // namespace pwdft
