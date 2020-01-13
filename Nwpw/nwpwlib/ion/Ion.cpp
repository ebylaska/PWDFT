/* Ions.C - 
   Author - Eric Bylaska
*/

using namespace std;

#include	<string.h>


#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cmath>
#include        <cstdlib>

#include	"Ion.hpp"


/*******************************
 *                             *
 *        ion_find_nkatm       *
 *                             *
 *******************************/
static void ion_tmpsymbols(RTDB& myrtdb, const int nion, char *symb, char *symtmp)
{
   int j,ii,i,matype,nelem;
   char date[50];

   if (!myrtdb.get_info("geometry:geometry:tags",&matype,&nelem,date))
   { printf("rtdb.get_info error: geometry:geometry:tags not found\n"); exit(99);}

   else
   {
      if (!myrtdb.get("geometry:geometry:tags",matype,nelem,symtmp))
      { printf("rtdb_get error: geometry:geometry:tags not found\n"); exit(99);}

      j = 0;
      for (ii=0; ii<nion; ++ii)
      {
        i = 3*ii;
        while ((symtmp[j]!= ((char) 0)) && (symtmp[j]!= ((char) 10)))
           symb[i++] = symtmp[j++];
        symb[i] = '\0';
        if (symtmp[j]==((char) 10)) ++j;
      }
   }
}

static int ion_find_nkatm(RTDB& myrtdb, const int nion)
{
   int nkatm,i,ii,j,found;
   char *symtmp,*symb;
   symtmp = new char[3*nion];
   symb   = new char[3*nion];

   ion_tmpsymbols(myrtdb,nion,symb,symtmp);

   /* determine nkatm */
   nkatm = 0;
   for (ii=0; ii<nion; ++ii)
   {
      found = 0;
      for (i=0; i<nkatm; ++i)
         if (strcmp(&symb[3*ii],&symtmp[3*i])==0) found=1;
      if (!found)
      {
         strcpy(&symtmp[3*nkatm],&symb[3*ii]);
         ++nkatm;
      }
   }

   delete [] symb;
   delete [] symtmp;

   return nkatm;
}

/*******************************
 *                             *
 *    ion_find_atomkatmnatm    *
 *                             *
 *******************************/
static void ion_find_atomkatmnatm(RTDB& myrtdb, const int nion, const int nkatm,
                               char *atom, int *katm, int *natm)
{
   int i,ii,j,jj,found;
   char *symtmp,*symb;
   symtmp = new char[3*nion];
   symb   = new char[3*nion];
   ion_tmpsymbols(myrtdb,nion,symb,symtmp);

   /* determine atom */
   jj = 0;
   for (ii=0; ii<nion; ++ii)
   {
      found = 0;
      for (i=0; i<jj; ++i)
         if (strcmp(&symb[3*ii],&atom[3*i])==0) found=1;
      if (!found)
      {
         strcpy(&atom[3*jj],&symb[3*ii]);
         ++jj;
      }
   }

   /* determine katm */
   for (ii=0; ii<nion; ++ii)
   {
      j = 0;
      for (i=0; i<nkatm; ++i)
         if (strcmp(&symb[3*ii],&atom[3*i])==0) j=i;
      katm[ii] = j;
   }

   /* determine natm */
   for (i=0; i<nkatm; ++i)   natm[i] = 0;
   for (ii=0; ii<nion; ++ii) natm[katm[ii]] += 1;
}



/* Constructors */

/*********************************
 *                               *
 *          Ion::Ion             *
 *                               *
 *********************************/
Ion::Ion(RTDB& myrtdb) 
{

   int matype,nelem,ii,i,j,found;
   char *symtmp,*symb;
   char date[50];

   /* get parallel mappings */
   if (!myrtdb.get("geometry:geometry:ncenter",rtdb_int,1,&nion)) nion = 1;

   nkatm = ion_find_nkatm(myrtdb,nion);

   charge    = new double[nion];
   mass      = new double[nion];
   rion1     = new double[3*nion];
   rion2     = new double[3*nion];
   katm      = new int[nion];
   natm      = new int[nkatm];
   atomarray = new char[3*nkatm];

   ion_find_atomkatmnatm(myrtdb,nion,nkatm,atomarray,katm,natm);

   /*  read in ion positions */
   if (!myrtdb.get("geometry:geometry:coords",rtdb_double,3*nion,rion1)) 
   { rion1[0] = 0.0; rion1[1] = 0.0; rion1[2] = 0.0; }
   if (!myrtdb.get("geometry:geometry:coords",rtdb_double,3*nion,rion2)) 
   { rion2[0] = 0.0; rion2[1] = 0.0; rion2[2] = 0.0; }

   /*  read in masses and charges */
   if (!myrtdb.get("geometry:geometry:masses",rtdb_double,nion,mass)) mass[0] = 1.0;
   if (!myrtdb.get("geometry:geometry:charges",rtdb_double,nion,charge)) charge[0] = 1.0;

}

