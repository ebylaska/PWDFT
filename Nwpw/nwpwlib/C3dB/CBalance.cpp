/* CBalance.cpp
   Author - Eric Bylaska

*/

/**
 * @class CBalance
 * @brief Class for load balancing in parallel computing.
 *
 * The `Balance` class provides load balancing capabilities for parallel computations within
 * the parallel computing framework. It includes methods for load balancing, unbalancing, and
 * other related operations.
 *
 * @author Eric Bylaska
 * @date  2023-09-20
 *
 * @note Detailed function descriptions and implementation can be found in corresponding source files.
 *
 * @see Parallel.hpp
 */


#include "CBalance.hpp"
#include	<iostream>

namespace pwdft {

/********************************
 *                              *
 *         nwave2_sort          *
 *                              *
 ********************************/
/**
 * @brief Sorts an array of integers and returns the corresponding indices.
 *
 * This function sorts an array of integers in ascending order and returns the
 * indices that correspond to the sorted order.
 *
 * @param n The number of elements in the array.
 * @param f Pointer to the array of integers to be sorted.
 * @param indx Pointer to an integer array to store the sorted indices.
 */
static void nwave2_sort(const int n, const int *f, int *indx) {
  int i, j, idum;
  for (i = 0; i < n; ++i)
    indx[i] = i;

  for (i = 0; i < (n - 1); ++i)
    for (j = i + 1; j < n; ++j)
      if (f[indx[j]] < f[indx[i]]) {
        idum = indx[i];
        indx[i] = indx[j];
        indx[j] = idum;
      }
}


/********************************
 *                              *
 *       balance_init_a         *
 *                              *
 ********************************/
/**
 * @brief Initializes load balancing for a parallel application.
 *
 * This function initializes load balancing for a parallel application by
 * redistributing a certain workload among different processes.
 *
 * @param parall Pointer to the Parallel object for communication.
 * @param nwavein The workload to be redistributed.
 * @param np The number of processes.
 * @param taskid The ID of the current process.
 * @param nwave_out Output variable to store the redistributed workload.
 * @param npacket_out Output variable to store the number of packets.
 * @param receiver_out Output variable to indicate if the process is a receiver.
 * @param sender_out Output variable to indicate if the process is a sender.
 * @param proc_to Array to store process IDs to send packets to.
 * @param proc_from Array to store process IDs to receive packets from.
 * @param packet_size Array to store the size of each packet.
 * @param indx_start Array to store the starting index of packets.
 */
static void balance_init_a(Parallel *parall, const int nwavein, const int np,
                    const int taskid, int *nwave_out, int *npacket_out,
                    int *receiver_out, int *sender_out, int *proc_to,
                    int *proc_from, int *packet_size, int *indx_start) 
{

   int i, j, nwave, ave, below, above, ishort, ilong, done, npacket, receiver, sender;
   int *nwave2 = new int[np];
   int *indx = new int[np];

   /* define nwave2 */
   for (i=0; i<np; ++i)
      nwave2[i] = 0;
   nwave2[taskid] = nwavein;
   parall->Vector_ISumAll(1, np, nwave2);

   /* get the sorting index */
   nwave2_sort(np,nwave2,indx);

   /* get the average */
   ave = 0;
   for (i=0; i<np; ++i)
      ave += nwave2[i];
   ave = ave / np;

   /* get below */
   below = -1;
   while (nwave2[indx[below + 1]] < ave)
      ++below;

   /* get above */
   above = np;
   while (nwave2[indx[above - 1]] > ave)
      --above;

   npacket = 0;
   receiver = 0;
   sender = 0;
   if (np > 1) 
   {
      i = 0;
      j = np - 1;
      done = 0;
      if (i > below)
         done = 1;
      if (j < above)
         done = 1;
      while (!done) 
      {
         ishort = ave - nwave2[indx[i]];
         ilong = nwave2[indx[j]] - ave;
        
         if (taskid == indx[i]) 
         {
            ++npacket;
            proc_from[npacket - 1] = indx[j];
            receiver = 1;
         }
         if (taskid == indx[j]) 
         {
            ++npacket;
            proc_to[npacket - 1] = indx[i];
            sender = 1;
         }
        
         if (ishort == ilong) 
         {
            if (taskid == indx[i]) 
            {
               packet_size[npacket - 1] = ishort;
               indx_start[npacket - 1] = nwave2[indx[i]];
            }
            if (taskid == indx[j]) 
            {
               packet_size[npacket - 1] = ilong;
               indx_start[npacket - 1] = nwave2[indx[j]] - ilong;
            }
            nwave2[indx[i]] += ishort;
            nwave2[indx[j]] -= ilong;
            ++i;
            --j;
         } 
         else if (ishort < ilong) 
         {
            if (taskid == indx[i]) 
            {
               packet_size[npacket - 1] = ishort;
               indx_start[npacket - 1] = nwave2[indx[i]];
            }
            if (taskid == indx[j]) 
            {
               packet_size[npacket - 1] = ishort;
               indx_start[npacket - 1] = nwave2[indx[j]] - ishort;
            }
            nwave2[indx[i]] += ishort;
            nwave2[indx[j]] -= ishort;
            ++i;
         } 
         else if (ishort > ilong) 
         {
            if (taskid == indx[i]) 
            {
               packet_size[npacket - 1] = ilong;
               indx_start[npacket - 1] = nwave2[indx[i]];
            }
            if (taskid == indx[j]) 
            {
               packet_size[npacket - 1] = ilong;
               indx_start[npacket - 1] = nwave2[indx[j]] - ilong;
            }
            nwave2[indx[i]] += ilong;
            nwave2[indx[j]] -= ilong;
            --j;
         }
        
         if (i > below) done = 1;
         if (j < above) done = 1;
      }
   }

   *nwave_out = nwave2[taskid];
   *npacket_out = npacket;
   *receiver_out = receiver;
   *sender_out = sender;

   delete[] indx;
   delete[] nwave2;
}

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/**
 * @brief Constructor for the CBalance class.
 *
 * This constructor initializes the `CBalance` object by allocating memory for
 * various data structures used in balancing operations.
 *
 * @param inparall A pointer to the `Parallel` object.
 * @param maxsize0 The maximum size for balancing operations.
 * @param nidb An array containing wavefunction information.
 * @param nidb_out An array to store modified wavefunction information.
 */
CBalance::CBalance(Parallel *inparall, const int maxsize0, const int *nidb, int *nidb_out) 
{
   parall = inparall;
   maxsize = maxsize0;

   proc_to_list     = new (std::nothrow) int*[maxsize]();
   proc_from_list   = new (std::nothrow) int*[maxsize]();
   packet_size_list = new (std::nothrow) int*[maxsize]();
   indx_start_list  = new (std::nothrow) int*[maxsize]();

   npacket_list     = new (std::nothrow) int[maxsize]();
   receiver_list    = new (std::nothrow) int[maxsize]();
   sender_list      = new (std::nothrow) int[maxsize]();

 
   int nb, nwave, nwave_out;
   int np = parall->np_i();
   int taskid = parall->taskid_i();
 
   for (nb=0; nb<maxsize; ++nb) 
   {
      proc_to_list[nb] = new int[np];
      proc_from_list[nb] = new int[np];
      packet_size_list[nb] = new int[np];
      indx_start_list[nb] = new int[np];
   }
 
   for (nb=0; nb<maxsize; ++nb) 
   {
      nwave = nidb[nb];
      balance_init_a(parall, nwave, np, taskid, &nwave_out, &npacket_list[nb],
                     &receiver_list[nb], &sender_list[nb], proc_to_list[nb],
                     proc_from_list[nb], packet_size_list[nb],
                     indx_start_list[nb]);
      nidb_out[nb] = nidb[nb] + (nwave_out - nwave);
   }
}

/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
/**
 * @brief Destructor for the CBalance class.
 *
 * This destructor releases the memory allocated for various data structures
 * used for balancing operations.
 */
CBalance::~CBalance() 
{
   for (int nb=0; nb<maxsize; ++nb) 
   {
      delete[] proc_to_list[nb];
      delete[] proc_from_list[nb];
      delete[] packet_size_list[nb];
      delete[] indx_start_list[nb];
   }
   delete[] proc_to_list;
   delete[] proc_from_list;
   delete[] packet_size_list;
   delete[] indx_start_list;

   delete[] npacket_list;
   delete[] receiver_list;
   delete[] sender_list;
}


/***********************************
 *                                 *
 *   CBalance::c_unbalance_start   *
 *                                 *
 ***********************************/
/**
 * @brief Start continuous data unbalancing process.
 *
 * This function initiates the continuous data unbalancing process by sending
 * and receiving data based on the provided parameters.
 *
 * @param nb The number of packets for which data unbalancing is performed.
 * @param a An array of double precision data to unbalance.
 * @param request_indx The index of the communication request to track.
 * @param msgtype The message type for communication.
 */
void CBalance::c_unbalance_start(const int nffts, const int nb, double *a, const int request_indx, const int msgtype) 
{
   //int j, pto, pfrom, msglen, indx;
   if (sender_list[nb])
      for (auto j=0; j<npacket_list[nb]; ++j) 
      {
         int pfrom  = proc_to_list[nb][j];
         int msglen = 2*nffts*packet_size_list[nb][j];
         int indx   = 2*nffts*indx_start_list[nb][j];
         if (msglen > 0)
            parall->adreceive(request_indx, msgtype, pfrom, msglen, a + indx);
      }
 
   if (receiver_list[nb])
      for (auto j=0; j<npacket_list[nb]; ++j) 
      {
         int pto    = proc_from_list[nb][j];
         int msglen = 2*nffts*packet_size_list[nb][j];
         int indx   = 2*nffts*indx_start_list[nb][j];
         if (msglen > 0)
            parall->adsend(request_indx, msgtype, pto, msglen, a + indx);
      }
}


/***********************************
 *                                 *
 *    CBalance::c_unbalance_end    *
 *                                 *
 ***********************************/
/**
 * @brief Finish continuous data unbalancing and wait for completion.
 *
 * This function completes the continuous data unbalancing process and waits for
 * all outstanding communication requests to finish.
 *
 * @param nb The number of packets for which data unbalancing was performed.
 * @param a An array of double precision data used for unbalancing.
 * @param request_indx The index of the communication request to wait for.
 */
void CBalance::c_unbalance_end(const int nffts, const int nb, double *a, const int request_indx) {
  parall->awaitall(request_indx);
}


/********************************
 *                              *
 *    CBalance::c_unbalance     *
 *                              *
 ********************************/
/**
 * @brief Unbalance continuous data among parallel processes.
 *
 * This function unbalances continuous data among parallel processes by exchanging
 * data segments to achieve load unbalancing.
 *
 * @param nb The number of packets for which data unbalancing is performed.
 * @param a An array of double precision data to be unbalanced.
 */
void CBalance::c_unbalance(const int nb, double *a) {
  int j, pto, pfrom, msglen, indx;

  if (sender_list[nb])
    for (j = 0; j < npacket_list[nb]; ++j) {
      pfrom = proc_to_list[nb][j];
      msglen = 2 * packet_size_list[nb][j];
      indx = 2 * indx_start_list[nb][j];
      if (msglen > 0)
        parall->dreceive(1, 9, pfrom, msglen, &a[indx]);
    }

  if (receiver_list[nb])
    for (j = 0; j < npacket_list[nb]; ++j) {
      pto = proc_from_list[nb][j];
      msglen = 2 * packet_size_list[nb][j];
      indx = 2 * indx_start_list[nb][j];
      if (msglen > 0)
        parall->dsend(1, 9, pto, msglen, &a[indx]);
    }
}

/********************************
 *                              *
 *     CBalance::c_balance      *
 *                              *
 ********************************/
/**
 * @brief CBalance continuous data among parallel processes.
 *
 * This function balances continuous data among parallel processes by exchanging
 * data segments to achieve load balancing.
 *
 * @param nb The number of packets for which data balancing is performed.
 * @param a An array of double precision data to be balanced.
 */
void CBalance::c_balance(const int nb, double *a) 
{
   int j, pto, pfrom, msglen, indx;
 
   if (sender_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pto = proc_to_list[nb][j];
         msglen = 2 * packet_size_list[nb][j];
         indx = 2 * indx_start_list[nb][j];
         if (msglen > 0)
            parall->dsend(1, 9, pto, msglen, a+indx);
      }
 
   if (receiver_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pfrom = proc_from_list[nb][j];
         msglen = 2 * packet_size_list[nb][j];
         indx = 2 * indx_start_list[nb][j];
         if (msglen > 0)
            parall->dreceive(1, 9, pfrom, msglen, a+indx);
      }
}

/***********************************
 *                                 *
 *   CBalance::c_balance_start     *
 *                                 *
 ***********************************/
/**
 * @brief Initiate communication for balancing continuous data.
 *
 * This function initiates communication operations for balancing continuous data.
 * It is used to send and receive data segments to achieve load balancing among
 * parallel processes.
 *
 * @param nb The number of packets for which communication is initiated.
 * @param a An array of double precision data to be balanced.
 * @param request_indx The index of the communication request to track the progress.
 * @param msgtype The message type for communication.
 */
void CBalance::c_balance_start(const int nffts, const int nb, double *a, const int request_indx, const int msgtype) 
{
   int j, pto, pfrom, msglen, indx;
 
   if (sender_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pto = proc_to_list[nb][j];
         msglen = 2 * nffts * packet_size_list[nb][j];
         indx = 2 * nffts * indx_start_list[nb][j];
         if (msglen > 0)
            parall->adsend(request_indx, msgtype, pto, msglen, a + indx);
      }
 
   if (receiver_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pfrom = proc_from_list[nb][j];
         msglen = 2 * nffts *  packet_size_list[nb][j];
         indx = 2 *  nffts * indx_start_list[nb][j];
         if (msglen > 0)
            parall->adreceive(request_indx, msgtype, pfrom, msglen, a + indx);
      }
}

/***********************************
 *                                 *
 *    CBalance::c_balance_end      *
 *                                 *
 ***********************************/
/**
 * @brief Wait for the completion of communication for balancing continuous data.
 *
 * This function waits for the completion of communication operations initiated
 * by `c_balance_start` for balancing continuous data. It is used to ensure that
 * the data distribution is properly balanced among parallel processes before
 * continuing with further computations.
 *
 * @param nb The number of packets for which communication was initiated.
 * @param a An array of double precision data to be balanced.
 * @param request_indx The index of the communication request to await.
 */
void CBalance::c_balance_end(const int nffts, const int nb, double *a, const int request_indx) 
{
  parall->awaitall(request_indx);
}

/********************************
 *                              *
 *    CBalance::r_unbalance     *
 *                              *
 ********************************/
/**
 * @brief Unbalance the distribution of double precision data among parallel processes.
 *
 * This function unbalances the distribution of double precision data among parallel processes.
 * It sends and receives data between processes to achieve an unbalanced distribution.
 *
 * @param nb The number of packets to unbalance.
 * @param a An array of double precision data to unbalance.
 *
 * @note The sender and receiver lists for each packet are determined by the
 *       previous configuration.
 */
void CBalance::r_unbalance(const int nb, double *a) 
{
   int j, pto, pfrom, msglen, indx;
 
   if (sender_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pfrom = proc_to_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->dreceive(1, 9, pfrom, msglen, &a[indx]);
      }
 
   if (receiver_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pto = proc_from_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->dsend(1, 9, pto, msglen, &a[indx]);
      }
}


/********************************
 *                              *
 *    CBalance::t_unbalance     *
 *                              *
 ********************************/
/**
 * @brief Unbalance the distribution of double precision data among parallel processes.
 *
 * This function unbalances the distribution of double precision data among parallel processes.
 * It sends and receives data between processes to achieve an unbalanced distribution.
 *
 * @param nb The number of packets to unbalance.
 * @param a An array of double precision data to unbalance.
 *
 * @note The sender and receiver lists for each packet are determined by the
 *       previous configuration.
 */
void CBalance::t_unbalance(const int nb, double *a)
{
   int j, pto, pfrom, msglen, indx;

   if (sender_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j)
      {
         pfrom = proc_to_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->dreceive(1, 9, pfrom, msglen, &a[indx]);
      }

   if (receiver_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j)
      {
         pto = proc_from_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->dsend(1, 9, pto, msglen, &a[indx]);
      }
}



/********************************
 *                              *
 *    CBalance::r_balance       *
 *                              *
 ********************************/
/**
 * @brief CBalance the distribution of double precision data among parallel processes.
 *
 * This function balances the distribution of double precision data among parallel processes.
 * It sends and receives data between processes to achieve a balanced distribution.
 *
 * @param nb The number of packets to balance.
 * @param a An array of double precision data to balance.
 *
 * @note The sender and receiver lists for each packet are determined by the
 *       previous configuration.
 */
void CBalance::r_balance(const int nb, double *a) 
{
   int j, pto, pfrom, msglen, indx;
 
   if (sender_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pto = proc_to_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->dsend(1, 9, pto, msglen, &a[indx]);
      }
 
   if (receiver_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pfrom = proc_from_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->dreceive(1, 9, pfrom, msglen, &a[indx]);
      }
}



/********************************
 *                              *
 *    CBalance::t_balance       *
 *                              *
 ********************************/
/**
 * @brief CBalance the distribution of double precision data among parallel processes.
 *
 * This function balances the distribution of double precision data among parallel processes.
 * It sends and receives data between processes to achieve a balanced distribution.
 *
 * @param nb The number of packets to balance.
 * @param a An array of double precision data to balance.
 *
 * @note The sender and receiver lists for each packet are determined by the
 *       previous configuration.
 */
void CBalance::t_balance(const int nb, double *a)
{
   int j, pto, pfrom, msglen, indx;

   if (sender_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j)
      {
         pto = proc_to_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->dsend(1, 9, pto, msglen, &a[indx]);
      }

   if (receiver_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j)
      {
         pfrom = proc_from_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->dreceive(1, 9, pfrom, msglen, &a[indx]);
      }
}


/********************************
 *                              *
 *    CBalance::i_balance        *
 *                              *
 ********************************/
/**
 * @brief CBalance the distribution of integers among parallel processes.
 *
 * This function balances the distribution of integers among parallel processes.
 * It sends and receives data between processes to achieve a balanced distribution.
 *
 * @param nb The number of packets to balance.
 * @param a An array of integers to balance.
 *
 * @note The sender and receiver lists for each packet are determined by the
 *       previous configuration.
 */
void CBalance::i_balance(const int nb, int *a) 
{
   int j, pto, pfrom, msglen, indx;
 
   if (sender_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pto = proc_to_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->isend(1, 9, pto, msglen, &a[indx]);
      }
 
   if (receiver_list[nb])
      for (j = 0; j < npacket_list[nb]; ++j) 
      {
         pfrom = proc_from_list[nb][j];
         msglen = packet_size_list[nb][j];
         indx = indx_start_list[nb][j];
         if (msglen > 0)
            parall->ireceive(1, 9, pfrom, msglen, &a[indx]);
      }
}

} // namespace pwdft
