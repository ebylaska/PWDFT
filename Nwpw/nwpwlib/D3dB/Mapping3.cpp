/* Mapping3.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/

/*
#include        <iostream>
#include        <cmath>
#include        <cstdlib>
using namespace std;
*/

#include	"hilbert.hpp"
#include	"Mapping3.hpp"

/***********************************
 *                                 *
 *        generate_map_indexes     *
 *                                 *
 ***********************************/
int  generate_map_indexes(int taskid,int np,
                          int ny,    int nz,
                          int *pmap, int *qmap)
{
   int nq,q,p;
   int *indx_proc;
   int *indx_q;
   int nq1        = (ny*nz)/np;
   int rmdr1      = (ny*nz)%np;
   int nq2        = nq1;


   indx_proc = new int[ny*nz];
   indx_q    = new int[ny*nz];

   if (rmdr1>0) ++nq2;
   nq=0; p=0; q=0;
   for (int i=0; i<(ny*nz); ++i)
   {
      indx_proc[i] = p;
      indx_q[i]    = q;
      if (taskid==p) ++nq;
      ++q;
      if (q>=nq2)
      {
          q=0;
          ++p;
          p = p%np;
          if (p>=rmdr1) nq2=nq1;
      }
   }

   for (int k=0; k<nz; ++k)
   for (int j=0; j<ny; ++j)
   {
      p = indx_proc[pmap[j+k*ny]];
      q = indx_q[pmap[j+k*ny]];
      pmap[j+k*ny] = p;
      qmap[j+k*ny] = q;
   }
   delete [] indx_proc;
   delete [] indx_q;

   return nq;
}


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

Mapping3::Mapping3()
{
   maptype=0; n2ft3d=0; nfft3d=0;
   nx=0; ny=0; nz=0;
   qmap[0] = 0; qmap[1] = 0; qmap[2] = 0;
   pmap[0] = 0; pmap[1] = 0; pmap[2] = 0;
   kmap = 0;
}

Mapping3::Mapping3(const int mapin, const int npin, const int taskidin, 
                   const int nxin,  const int nyin, const int nzin)
{
   int k,p,q;

   qmap[0] = 0; qmap[1] = 0; qmap[2] = 0;
   pmap[0] = 0; pmap[1] = 0; pmap[2] = 0;
   kmap = 0;

   np      = npin;
   taskid  = taskidin;
   nx      = nxin;
   ny      = nyin;
   nz      = nzin;
   maptype = mapin;

   /* slab mapping */
   if (maptype==1)
   {
      qmap[0] =  new int[nz];
      pmap[0] =  new int[nz];
      kmap    =  new int[nz];

       /* cyclic mapping */
       p = 0; q = 0;
       for (k=0; k<nz; ++k)
       {
          qmap[0][k] = q;
          pmap[0][k] = p;
          if (p==taskid) nq=q+1;
          ++p;
          if (p>=np) {p=0; ++q;}
       }
       nfft3d=(nx/2+1)*ny*nq;
       n2ft3d=2*nfft3d;

   }

   /* hilbert or hcurve  mapping */
   else
   {
      qmap[0] =  new int[ny*nz];       pmap[0] =  new int[ny*nz];
      qmap[1] =  new int[nz*(nx/2+1)]; pmap[1] =  new int[nz*(nx/2+1)];
      qmap[2] =  new int[(nx/2+1)*ny]; pmap[2] =  new int[(nx/2+1)*ny];

      if (maptype==2)
      {
         hilbert2d_map(ny,nz,pmap[0]);
         hilbert2d_map(nz,(nx/2+1),pmap[1]);
         hilbert2d_map((nx/2+1),ny,pmap[2]);
      }

      if (maptype==3)
      {
         hcurve2d_map(ny,nz,pmap[0]);
         hcurve2d_map(nz,(nx/2+1),pmap[1]);
         hcurve2d_map((nx/2+1),ny,pmap[2]);
      }

      nq1 = generate_map_indexes(taskid,np,ny,nz,      pmap[0],qmap[0]);
      nq2 = generate_map_indexes(taskid,np,nz,(nx/2+1),pmap[1],qmap[1]);
      nq3 = generate_map_indexes(taskid,np,(nx/2+1),ny,pmap[2],qmap[2]);
      nfft3d = (nx/2+1)*nq1;
      if ((ny*nq2) > nfft3d) nfft3d=ny*nq2;
      if ((nz*nq3) > nfft3d) nfft3d=nz*nq3;
      n2ft3d = 2*nfft3d;
   }
}


/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
Mapping3::~Mapping3()
{
   for (int i=0; i<3; ++i)
   {
      if (qmap[i]) delete [] pmap[i];
      if (pmap[i]) delete [] qmap[i];
   }
   if (kmap) delete [] kmap;

   if ((maptype==2) || (maptype==3))
   {
        int *h_iq_to_i1[6],*h_iq_to_i2[6];
        int *h_i1_start[6],*h_i2_start[6];
   }


}

