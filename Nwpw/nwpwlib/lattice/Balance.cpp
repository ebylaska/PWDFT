/* Balance.C
   Author - Eric Bylaska

*/


#include	"Balance.hpp"

namespace pwdft {


void nwave2_sort(const int n,  const int *f, int *indx)
{
   int i,j,idum;
   for (i=0; i<n; ++i) 
      indx[i] = i;

   for (i=0; i<(n-1); ++i)
      for (j=i+1; j<n; ++j)
         if (f[indx[j]] < f[indx[i]])
         {
            idum = indx[i];
            indx[i] = indx[j];
            indx[j] = idum;
         }
}

void balance_init_a(Parallel *parall, 
             const int nwavein, const int np, const int taskid, 
             int *nwave_out,
             int *npacket_out, int *receiver_out, int *sender_out,
             int *proc_to, int *proc_from,
             int *packet_size, int *indx_start)
{

   int i,j,nwave,ave,below,above,ishort,ilong,done,npacket,receiver,sender;
   int *nwave2 = new int[np];
   int *indx   = new int[np];

   /* define nwave2 */
   for (i=0; i<np; ++i)   
      nwave2[i] = 0;
   nwave2[taskid] = nwavein;
   parall->Vector_ISumAll(1,np,nwave2);


   /* get the sorting index */
   nwave2_sort(np,nwave2,indx);

   /* get the average */
   ave = 0;
   for (i=0; i<np; ++i)
      ave += nwave2[i];
   ave = ave/np;

   /* get below */
   below = -1;
   while (nwave2[indx[below+1]] < ave)
      ++below;

   /* get above */
   above = np;
   while (nwave2[indx[above-1]] > ave)
      --above;

   npacket  = 0;
   receiver = 0;
   sender   = 0;
   if (np > 1) 
   {
      i = 0;
      j = np-1;
      done = 0;
      if (i > below) done = 1;
      if (j < above) done = 1;
      while (!done)
      {
         ishort = ave - nwave2[indx[i]];
         ilong  = nwave2[indx[j]] - ave;

         if (taskid==indx[i])
         {
            ++npacket;
            proc_from[npacket-1] = indx[j];
            receiver = 1;
         }
         if (taskid==indx[j])
         {
            ++npacket;
            proc_to[npacket-1] = indx[i];
            sender = 1;
         }


         if (ishort==ilong)
         {
            if (taskid==indx[i])
            {
               packet_size[npacket-1] = ishort;
               indx_start[npacket-1]  = nwave2[indx[i]];
            }
            if (taskid==indx[j])
            {
               packet_size[npacket-1] = ilong;
               indx_start[npacket-1]  = nwave2[indx[j]] - ilong;
            }
            nwave2[indx[i]] += ishort;
            nwave2[indx[j]] -= ilong;
            ++i;
            --j;
         }
         else if (ishort<ilong)
         {
            if (taskid==indx[i])
            {
               packet_size[npacket-1] = ishort;
               indx_start[npacket-1]  = nwave2[indx[i]];
            }
            if (taskid==indx[j])
            {
               packet_size[npacket-1] = ishort;
               indx_start[npacket-1]  = nwave2[indx[j]] - ishort;
            }
            nwave2[indx[i]] += ishort;
            nwave2[indx[j]] -= ishort;
            ++i;
         }
         else if (ishort>ilong)
         {
            if (taskid==indx[i])
            {
               packet_size[npacket-1] = ilong;
               indx_start[npacket-1]  = nwave2[indx[i]];
            }
            if (taskid==indx[j])
            {
               packet_size[npacket-1] = ilong;
               indx_start[npacket-1]  = nwave2[indx[j]] - ilong;
            }
            nwave2[indx[i]] += ilong;
            nwave2[indx[j]] -= ilong;
            --j;
         }
 
         if (i > below) done = 1;
         if (j < above) done = 1;
      }
   }

   *nwave_out    = nwave2[taskid];
   *npacket_out  = npacket;
   *receiver_out = receiver;
   *sender_out   = sender;

   delete [] indx;
   delete [] nwave2;
}

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/


Balance::Balance(Parallel *inparall, const int maxsize0, const int *nidb, int *nidb_out)
{
   parall = inparall;
   maxsize = maxsize0;

   int nb,nwave,nwave_out;
   int np     = parall->np_i();
   int taskid = parall->taskid_i();

   for (nb=0; nb<maxsize; ++nb)
   {
      proc_to_list[nb]     = new int[np];
      proc_from_list[nb]   = new int[np];
      packet_size_list[nb] = new int[np];
      indx_start_list[nb]  = new int[np];
   }

   for (nb=0; nb<maxsize; ++nb)
   {
       nwave = nidb[nb];
       balance_init_a(parall,nwave,np,taskid,
               &nwave_out,
               &npacket_list[nb],
               &receiver_list[nb],
               &sender_list[nb],
               proc_to_list[nb],
               proc_from_list[nb],
               packet_size_list[nb],
               indx_start_list[nb]);
       nidb_out[nb] = nidb[nb] + (nwave_out-nwave);
   }
}



/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
Balance::~Balance()
{
   for (int nb=0; nb<maxsize; ++nb)
   {
      delete [] proc_to_list[nb];
      delete [] proc_from_list[nb];
      delete [] packet_size_list[nb];
      delete [] indx_start_list[nb];
   }
}


/***********************************
 *                                 *
 *    Balance::c_unbalance_start   *
 *                                 *
 ***********************************/
void Balance::c_unbalance_start(const int nb, double *a, const int request_indx, const int msgtype)
{
   int j,pto,pfrom,msglen,indx;

   if (sender_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pfrom  = proc_to_list[nb][j];
         msglen = 2*packet_size_list[nb][j];
         indx   = 2*indx_start_list[nb][j];
         if (msglen>0)
            parall->adreceive(request_indx,msgtype,pfrom,msglen,a+indx);
      }

   if (receiver_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pto    = proc_from_list[nb][j];
         msglen = 2*packet_size_list[nb][j];
         indx   = 2*indx_start_list[nb][j];
         if (msglen>0)
            parall->adsend(request_indx,msgtype,pto,msglen,a+indx);
      }
}
/***********************************
 *                                 *
 *    Balance::c_unbalance_end     *
 *                                 *
 ***********************************/
void Balance::c_unbalance_end(const int nb, double *a, const int request_indx)
{
   parall->awaitall(request_indx);
}


/********************************
 *                              *
 *    Balance::c_unbalance      *
 *                              *
 ********************************/

void Balance::c_unbalance(const int nb, double *a)
{
   int j,pto,pfrom,msglen,indx;

   if (sender_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pfrom  = proc_to_list[nb][j];
         msglen = 2*packet_size_list[nb][j];
         indx   = 2*indx_start_list[nb][j];
         if (msglen>0) 
            parall->dreceive(1,9,pfrom,msglen,&a[indx]);
      }
       
   if (receiver_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pto    = proc_from_list[nb][j];
         msglen = 2*packet_size_list[nb][j];
         indx   = 2*indx_start_list[nb][j];
         if (msglen>0) 
            parall->dsend(1,9,pto,msglen,&a[indx]);
      }
}


/********************************
 *                              *
 *    Balance::c_balance        *
 *                              *
 ********************************/

void Balance::c_balance(const int nb, double *a)
{
   int j,pto,pfrom,msglen,indx;

   if (sender_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pto    = proc_to_list[nb][j];
         msglen = 2*packet_size_list[nb][j];
         indx   = 2*indx_start_list[nb][j];
         if (msglen>0) 
            parall->dsend(1,9,pto,msglen,&a[indx]);
      }
       
   if (receiver_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pfrom  = proc_from_list[nb][j];
         msglen = 2*packet_size_list[nb][j];
         indx   = 2*indx_start_list[nb][j];
         if (msglen>0) 
            parall->dreceive(1,9,pfrom,msglen,&a[indx]);
      }
}


/********************************
 *                              *
 *    Balance::t_unbalance      *
 *                              *
 ********************************/

void Balance::t_unbalance(const int nb, double *a)
{
   int j,pto,pfrom,msglen,indx;

   if (sender_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pfrom  = proc_to_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx   = indx_start_list[nb][j];
         if (msglen>0)
            parall->dreceive(1,9,pfrom,msglen,&a[indx]);
      }

   if (receiver_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pto    = proc_from_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx   = indx_start_list[nb][j];
         if (msglen>0)
            parall->dsend(1,9,pto,msglen,&a[indx]);
      }
}


/********************************
 *                              *
 *    Balance::t_balance        *
 *                              *
 ********************************/

void Balance::t_balance(const int nb, double *a)
{
   int j,pto,pfrom,msglen,indx;

   if (sender_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pto    = proc_to_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx   = indx_start_list[nb][j];
         if (msglen>0)
            parall->dsend(1,9,pto,msglen,&a[indx]);
      }

   if (receiver_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pfrom  = proc_from_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx   = indx_start_list[nb][j];
         if (msglen>0)
            parall->dreceive(1,9,pfrom,msglen,&a[indx]);
      }
}



/********************************
 *                              *
 *    Balance::i_balance        *
 *                              *
 ********************************/

void Balance::i_balance(const int nb, int *a)
{
   int j,pto,pfrom,msglen,indx;

   if (sender_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pto    = proc_to_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx   = indx_start_list[nb][j];
         if (msglen>0)
            parall->isend(1,9,pto,msglen,&a[indx]);
      }

   if (receiver_list[nb])
      for (j=0; j<npacket_list[nb]; ++j)
      {
         pfrom  = proc_from_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx   = indx_start_list[nb][j];
         if (msglen>0)
            parall->ireceive(1,9,pfrom,msglen,&a[indx]);
      }
}

}
