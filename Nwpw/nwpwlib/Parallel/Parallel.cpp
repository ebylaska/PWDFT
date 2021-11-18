/* Parallel.C
   Author - Eric Bylaska
	this class is used for defining 3d parallel maps
*/

#include        <cmath>
#include        <cstdlib>
using namespace std;

#include	"mpi.h"
#include	"Parallel.hpp"
#include	"Control2.hpp"

namespace pwdft {
using namespace pwdft;

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
//Parallel::Parallel(int argc, char *argv[]) 
Parallel::Parallel(MPI_Comm comm_world0)
{

   int indx,d,dd,ii,np0,done;

   //MPI::Init(argc,argv);
   //MPI_Init(&argc,&argv);

   dim = 1;

   //comm_i[0]   = MPI_COMM_WORLD;

   comm_world  = comm_world0;
   comm_i[0]   = comm_world;
   MPI_Comm_size(comm_i[0],&npi[0]);
   MPI_Comm_rank(comm_i[0],&taskidi[0]);

   comm_i[1]   = comm_i[0];
   comm_i[2]   = comm_i[0];
   npi[1]     = npi[0];
   npi[2]     = 1;
   taskidi[1] = taskidi[0];
   taskidi[2] = 0;

   procNd      = new int[npi[0]]();

   npi[1] = npi[0];
   npi[2] = 1;
   for (int i=0; i<npi[0]; ++i)
      procNd[i] = i;

}

void Parallel::init2d(const int ncolumns, const int pfft3_qsize)
{
   int ii;
   MPI_Group orig_group;

   MPI_Barrier(comm_world);
   if (ncolumns>1)
   {
      dim = 2;
      npi[1] = npi[0]/ncolumns;
      npi[2] = ncolumns;

      int icount = 0;
      for (int j=0; j<npi[2]; ++j)
      for (int i=0; i<npi[1]; ++i)
      {
         if (icount==taskidi[0])
         {
            taskidi[1] = i;
            taskidi[2] = j;
         }
         procNd[i+j*npi[1]] = icount;
         icount = (icount+1)%npi[0];
      }

      int *itmp = new int[npi[0]]();

      for (int i=0; i<npi[1]; ++i) itmp[i] = procNd[i+taskidi[2]*npi[1]];
      MPI_Comm_group(comm_world,&orig_group);
      MPI_Group_incl(orig_group,npi[1],itmp,&group_i[1]);
      MPI_Comm_create(comm_world,group_i[1],&comm_i[1]);


      for (int j=0; j<npi[2]; ++j) itmp[j] = procNd[taskidi[1]+j*npi[1]];
      MPI_Group_incl(orig_group,npi[2],itmp,&group_i[2]);
      MPI_Comm_create(comm_world,group_i[2],&comm_i[2]);

      delete [] itmp;
   }

   //ii = 3+control.pfft3_qsize();
   //request = new MPI::Request*[ii];
   ii = 3+pfft3_qsize;
   reqcnt   = new int[ii]();
   request  = new MPI_Request*[ii]();
   statuses = new MPI_Status*[ii]();

   MPI_Barrier(comm_world);

}

/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
Parallel::~Parallel()
{
   MPI_Barrier(comm_world);
   if (dim>1)
   {
      for (int d=1; d<=dim; ++d)
      {
         //group_i[d].Free();
        // comm_i[d].Free();
         MPI_Comm_free(&comm_i[d]);
         MPI_Group_free(&group_i[d]);
      }
   }

   delete [] procNd;

   delete [] reqcnt;
   delete [] request;
   delete [] statuses;


   MPI_Barrier(comm_world);
}
/********************************
 *                              *
 *          MaxAll              *
 *                              *
 ********************************/
double Parallel::MaxAll(const int d, const double sum)
{
   double sumout;
   if (npi[d]>1) 
      //comm_i[d].Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_MAX);
      MPI_Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_i[d]);
   else 
      sumout = sum;
   return sumout;
}

/********************************
 *                              *
 *          SumAll              *
 *                              *
 ********************************/
double Parallel::SumAll(const int d, const double sum)
{
   double sumout;

   if (npi[d]>1) 
      //comm_i[d].Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_SUM);
      MPI_Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm_i[d]);
   else
      sumout = sum;
   return sumout;
}

/********************************
 *                              *
 *         ISumAll              *
 *                              *
 ********************************/
int Parallel::ISumAll(const int d, const int sum)
{
   int sumout;

   if (npi[d]>1) 
      //comm_i[d].Allreduce(&sum,&sumout,1,MPI_INTEGER,MPI_SUM);
      MPI_Allreduce(&sum,&sumout,1,MPI_INTEGER,MPI_SUM,comm_i[d]);
   else
      sumout = sum;
   return sumout;
}


/********************************
 *                              *
 *       Vector_SumAll          *
 *                              *
 ********************************/
void Parallel::Vector_SumAll(const int d, const int n, double *sum)
{
   double *sumout;
   if (npi[d]>1)
   {
      sumout = new double [n];
      //comm_i[d].Allreduce(sum,sumout,n,MPI_DOUBLE_PRECISION,MPI_SUM);
      MPI_Allreduce(sum,sumout,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm_i[d]);
      for (int i=0; i<n; ++i) sum[i] = sumout[i];
      delete [] sumout;
   }
}

/********************************
 *                              *
 *       Vector_ISumAll         *
 *                              *
 ********************************/
void Parallel::Vector_ISumAll(const int d, const int n, int *sum)
{
   int *sumout;
   if (npi[d]>1)
   {
      sumout = new int [n];
      //comm_i[d].Allreduce(sum,sumout,n,MPI_INTEGER,MPI_SUM);
      MPI_Allreduce(sum,sumout,n,MPI_INTEGER,MPI_SUM,comm_i[d]);
      for (int i=0; i<n; ++i) sum[i] = sumout[i];
      delete [] sumout;
   }
}


/********************************
 *                              *
 *       Brdcst_Values          *
 *                              *
 ********************************/
void Parallel::Brdcst_Values(const int d, const int root, const int n, double *sum)
{
   //if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_DOUBLE_PRECISION,root);
   if (npi[d]>1) MPI_Bcast(sum,n,MPI_DOUBLE_PRECISION,root,comm_i[d]);
}

/********************************
 *                              *
 *       Brdcst_iValues         *
 *                              *
 ********************************/
void Parallel::Brdcst_iValues(const int d, const int root, const int n, int *sum)
{
   //if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_INTEGER,root);
   if (npi[d]>1) MPI_Bcast(sum,n,MPI_INTEGER,root,comm_i[d]);
}

/********************************
 *                              *
 *       Brdcst_iValue          *
 *                              *
 ********************************/
void Parallel::Brdcst_iValue(const int d, const int root, int *sum)
{
   //if (npi[d]>1) comm_i[d].Bcast(sum,1,MPI_INTEGER,root);
   if (npi[d]>1) MPI_Bcast(sum,1,MPI_INTEGER,root,comm_i[d]);
}

/********************************
 *                              *
 *       Brdcst_cValues         *
 *                              *
 ********************************/
void Parallel::Brdcst_cValues(const int d, const int root, const int n, void *sum)
{
   //if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_CHAR,root);
   if (npi[d]>1) MPI_Bcast(sum,n,MPI_CHAR,root,comm_i[d]);
}

/********************************
 *                              *
 *       Parallel::dsend        *
 *                              *
 ********************************/
void Parallel::dsend(const int d, const int tag, const int procto, const int n, double *sum)
{
   //if (npi[d]>1) comm_i[d].Send(sum,n,MPI_DOUBLE_PRECISION,procto,tag);
   if (npi[d]>1) MPI_Send(sum,n,MPI_DOUBLE_PRECISION,procto,tag,comm_i[d]);
}


/********************************
 *                              *
 *       Parallel::dreceive     *
 *                              *
 ********************************/
void Parallel::dreceive(const int d, const int tag, const int procfrom, const int n, double *sum)
{
   //MPI::Status status;
   //if (npi[d]>1) comm_i[d].Recv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag);
   MPI_Status status;
   if (npi[d]>1) MPI_Recv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag,comm_i[d],&status);
}



/********************************
 *                              *
 *       Parallel::isend        *
 *                              *
 ********************************/
void Parallel::isend(const int d, const int tag, const int procto, const int n, int *sum)
{
   //if (npi[d]>1) comm_i[d].Send(sum,n,MPI_INTEGER,procto,tag);
   if (npi[d]>1) MPI_Send(sum,n,MPI_INTEGER,procto,tag,comm_i[d]);
}  


/********************************
 *                              *
 *       Parallel::ireceive     *
 *                              *
 ********************************/
void Parallel::ireceive(const int d, const int tag, const int procfrom, const int n, int *sum)
{
   //MPI::Status status;
   //if (npi[d]>1) comm_i[d].Recv(sum,n,MPI_INTEGER,procfrom,tag);
   MPI_Status status;
   if (npi[d]>1) MPI_Recv(sum,n,MPI_INTEGER,procfrom,tag,comm_i[d],&status);
}  




/********************************
 *                              *
 *       Parallel::astart       *
 *                              *
 ********************************/
void Parallel::astart(const int d, const int sz)
{
   if (npi[d]>1)
   {
      reqcnt[d]  = 0;
      //request[d] = new MPI::Request[sz];
      request[d]  = new MPI_Request[sz];
      statuses[d] = new MPI_Status[sz];
   }
}
/********************************
 *                              *
 *       Parallel::aend         *
 *                              *
 ********************************/
void Parallel::aend(const int d)
{
   //MPI::Status status[reqcnt[d]];
   //request[d][0].Waitall(reqcnt[d],request[d],status);
   if (npi[d]>1)
   {
      //request[d][0].Waitall(reqcnt[d],request[d]);
      MPI_Waitall(reqcnt[d],request[d],statuses[d]);
      delete [] request[d];
      delete [] statuses[d];
   }
}

/********************************
 *                              *
 *       Parallel::adreceive     *
 *                              *
 ********************************/
void Parallel::adreceive(const int d, const int tag, const int procfrom, const int n, double *sum)
{
   //MPI::Status status;
   //if (npi[d]>1) request[d][reqcnt[d]++] = comm_i[d].Irecv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag);
   if (npi[d]>1) MPI_Irecv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag,comm_i[d],&request[d][reqcnt[d]++]);
}

/********************************
 *                              *
 *       Parallel::adsend        *
 *                              *
 ********************************/
void Parallel::adsend(const int d, const int tag, const int procto, const int n, double *sum)
{
   //if (npi[d]>1) request[d][reqcnt[d]++] = comm_i[d].Isend(sum,n,MPI_DOUBLE_PRECISION,procto,tag);
   if (npi[d]>1) MPI_Isend(sum,n,MPI_DOUBLE_PRECISION,procto,tag,comm_i[d],&request[d][reqcnt[d]++]);
}

}
