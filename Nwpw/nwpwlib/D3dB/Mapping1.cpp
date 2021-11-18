/* Mapping1.C
   Author - Eric Bylaska

	this class is used for defining 1d parallel maps
*/
/*
#include        <iostream>
#include        <cmath>
#include        <cstdlib>
using namespace std;
*/

#include        <iostream>
#include	"Mapping1.hpp"

namespace pwdft {
using namespace pwdft;


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

Mapping1::Mapping1()
{
   maptype1=0; 
   ne[0]=0; neq[0]=0;
   ne[1]=0; neq[1]=0;
   qmap[0] = (int *) 0; qmap[1] = (int *) 0;
   pmap[0] = (int *) 0; pmap[1] = (int *) 0;
   kmap[0] = (int *) 0; kmap[1] = (int *) 0;
   nqarray[0] = (int *) 0; nqarray[1] = (int *) 0;
}

Mapping1::Mapping1(const int mapin, const int npin, const int taskidin, 
                   const int ispinin, const int *nein)
{
   //int k,p,q;


   np      = npin;
   taskid = taskidin;
   maptype1 = mapin;
   ispin = ispinin;
   
   ne[0]=0; neq[0]=0;
   ne[1]=0; neq[1]=0;

   for (auto ms=0; ms<ispin; ++ms)
   {
      ne[ms]      = nein[ms];
      qmap[ms]    = new (std::nothrow) int[ne[ms]];
      pmap[ms]    = new (std::nothrow) int[ne[ms]];
      kmap[ms]    = new (std::nothrow) int[ne[ms]];
      nqarray[ms] = new (std::nothrow) int[np];

      /* cyclic mapping */
      if (maptype1==0)
      {

         /* cyclic mapping */
         int p = 0; 
         int q = 0;
         for (auto k=0; k<ne[ms]; ++k)
         {
             qmap[ms][k] = q;
             pmap[ms][k] = p;

             if (p==taskid) neq[ms] = q;
             p++; 
             if (p>=np) { p = 0; ++q; }
         }

      }

      /* block mapping */
      else
      {
         int k,q;
         int p;
         for (p=0; p<np; ++p) nqarray[ms][p] = 0;
         p = 0;
         for (k=0; k<ne[ms]; ++k)
         {
            nqarray[ms][p] += 1;
            p = (p+1)%np;
         }

         k = 0;
         for (p=0; p<np; ++p)
            for (q=0; q<nqarray[ms][p]; ++q)
            {
               qmap[ms][k] = q;
               pmap[ms][k] = p;
               ++k;
            }
         neq[ms] = nqarray[ms][taskid];
      }

      for (auto k=0; k<ne[ms]; ++k)
        if (pmap[ms][k]==taskid)
           kmap[ms][qmap[ms][k]] = k;

   }

}

}
