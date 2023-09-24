/* Parallel.cpp
   Author - Eric Bylaska
        this class is used for defining 3d parallel maps
*/


/**
 * @class Parallel
 * @brief Provides MPI-based parallelism and communication management.
 *
 * The `Parallel` class encapsulates functionality for managing MPI-based parallel
 * computations. It allows users to distribute tasks across multiple MPI processes,
 * perform collective operations, and handle data communication between processes.
 *
 * This class is designed to work in a multi-dimensional parallel computing environment,
 * supporting both 1D and 2D process grids.
 */


#include <iostream> // DEBUG

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "Control2.hpp"
#include "Parallel.hpp"
#include "mpi.h"

#define MASTER 0

namespace pwdft {

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
// Parallel::Parallel(int argc, char *argv[])
/**
 * @brief Construct a Parallel object using a specified MPI communicator.
 *
 * This constructor initializes a Parallel object with the provided MPI communicator `comm_world0`.
 * The communicator is used for parallel communication among processes.
 *
 * @param[in] comm_world0 The MPI communicator to use for parallel communication.
 *
 * @note
 * - This constructor initializes the Parallel object and sets up the MPI environment.
 * - The provided MPI communicator is used to determine the number of processes and their ranks.
 * - The constructor sets the dimension `dim` to 1 and initializes other communication-related variables.
 * - The variable `base_stdio_print` is set to `is_master` based on the rank of the process.
 * - Use this constructor to create a Parallel object for parallelized computations.
 */
Parallel::Parallel(MPI_Comm comm_world0) {

  int indx, d, dd, ii, np0, done;

  // MPI::Init(argc,argv);
  // MPI_Init(&argc,&argv);

  dim = 1;

  // comm_i[0]   = MPI_COMM_WORLD;

  comm_world = comm_world0;
  comm_i[0] = comm_world;
  MPI_Comm_size(comm_i[0], &npi[0]);
  MPI_Comm_rank(comm_i[0], &taskidi[0]);

  comm_i[1] = comm_i[0];
  comm_i[2] = comm_i[0];
  npi[1] = npi[0];
  npi[2] = 1;
  taskidi[1] = taskidi[0];
  taskidi[2] = 0;

  procNd = new int[npi[0]]();

  npi[1] = npi[0];
  npi[2] = 1;
  for (int i = 0; i < npi[0]; ++i)
    procNd[i] = i;

  /* set initial base_stdio_print to is_master */
  base_stdio_print = (taskidi[0] == MASTER);
}

/********************************
 *                              *
 *       Parallel::init2d       *
 *                              *
 ********************************/
/**
 * @brief Initialize a 2D parallel configuration and communication setup.
 *
 * This function initializes a 2D parallel configuration for parallel computations.
 * It sets up the MPI communicator, determines the number of processes, and assigns
 * ranks based on a 2D grid. Additionally, it allocates memory for MPI requests and
 * statuses.
 *
 * @param[in] ncolumns   The number of columns in the 2D grid configuration.
 * @param[in] pfft3_qsize The size parameter used in memory allocation.
 *
 * @note
 * - This function is typically used to set up a 2D parallel configuration when needed.
 * - It assigns ranks to processes based on a 2D grid with `ncolumns` columns.
 * - Memory for MPI requests and statuses is allocated for communication.
 */
void Parallel::init2d(const int ncolumns, const int pfft3_qsize) 
{
   int ii;
   MPI_Group orig_group;
 
   MPI_Barrier(comm_world);
   if (ncolumns > 1) {
     dim = 2;
     npi[1] = npi[0] / ncolumns;
     npi[2] = ncolumns;
 
     int icount = 0;
     for (int j = 0; j < npi[2]; ++j)
       for (int i = 0; i < npi[1]; ++i) {
         if (icount == taskidi[0]) {
           taskidi[1] = i;
           taskidi[2] = j;
         }
         procNd[i + j * npi[1]] = icount;
         icount = (icount + 1) % npi[0];
       }
 
     int *itmp = new int[npi[0]]();
 
     for (int i = 0; i < npi[1]; ++i)
       itmp[i] = procNd[i + taskidi[2] * npi[1]];
     MPI_Comm_group(comm_world, &orig_group);
     MPI_Group_incl(orig_group, npi[1], itmp, &group_i[1]);
     MPI_Comm_create(comm_world, group_i[1], &comm_i[1]);
 
     for (int j = 0; j < npi[2]; ++j)
       itmp[j] = procNd[taskidi[1] + j * npi[1]];
     MPI_Group_incl(orig_group, npi[2], itmp, &group_i[2]);
     MPI_Comm_create(comm_world, group_i[2], &comm_i[2]);
 
     delete[] itmp;
   }
 
   // ii = 3+control.pfft3_qsize();
   // request = new MPI::Request*[ii];
   //          ii = 3 + (pfft3_qsize+6);
   max_reqstat = 4 + (pfft3_qsize+6);
   reqcnt = new int[max_reqstat]();
   request = new MPI_Request *[max_reqstat]();
   statuses = new MPI_Status *[max_reqstat]();
 
   MPI_Barrier(comm_world);
}


/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
/**
 * @brief Destructor for the Parallel class.
 *
 * This destructor cleans up resources associated with the Parallel class,
 * including MPI communicators and groups, memory allocated for MPI requests
 * and statuses, and the `procNd` array.
 *
 * @note
 * - This destructor should be called to release resources when the Parallel
 *   object is no longer needed.
 */
Parallel::~Parallel() 
{
   MPI_Barrier(comm_world);
   if (dim > 1) {
     for (int d = 1; d <= dim; ++d) {
       // group_i[d].Free();
       // comm_i[d].Free();
       MPI_Comm_free(&comm_i[d]);
       MPI_Group_free(&group_i[d]);
     }
   }
 
   delete[] procNd;
 
   delete[] reqcnt;
   delete[] request;
   delete[] statuses;
 
   MPI_Barrier(comm_world);
}


/********************************
 *                              *
 *          MaxAll              *
 *                              *
 ********************************/
/**
 * @brief Compute the maximum value across all processes in the specified communicator.
 *
 * This function computes the maximum value of a double across all processes in
 * the specified communicator `d` and returns the result.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param sum The input double value for which the maximum is to be computed locally.
 * @return The maximum value of `sum` across all processes in communicator `d`.
 *
 * @note
 * - This function performs an MPI_Allreduce operation to obtain the maximum value
 *   across all processes in the communicator.
 */
double Parallel::MaxAll(const int d, const double sum) {
  double sumout;
  if (npi[d] > 1)
    // comm_i[d].Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_MAX);
    MPI_Allreduce(&sum, &sumout, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_i[d]);
  else
    sumout = sum;
  return sumout;
}


/********************************
 *                              *
 *          SumAll              *
 *                              *
 ********************************/
/**
 * @brief Compute the sum of values across all processes in the specified communicator.
 *
 * This function computes the sum of a double value across all processes in the specified
 * communicator `d` and returns the result.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param sum The input double value for which the sum is to be computed locally.
 * @return The sum of `sum` across all processes in communicator `d`.
 *
 * @note
 * - This function performs an MPI_Allreduce operation to obtain the sum of values
 *   across all processes in the communicator.
 */
double Parallel::SumAll(const int d, const double sum) {
  double sumout;

  if (npi[d] > 1)
    // comm_i[d].Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_SUM);
    MPI_Allreduce(&sum, &sumout, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_i[d]);
  else
    sumout = sum;
  return sumout;
}


/********************************
 *                              *
 *         ISumAll              *
 *                              *
 ********************************/
/**
 * @brief Compute the integer sum of values across all processes in the specified communicator.
 *
 * This function computes the sum of an integer value across all processes in the specified
 * communicator `d` and returns the result.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param sum The input integer value for which the sum is to be computed locally.
 * @return The sum of `sum` across all processes in communicator `d`.
 *
 * @note
 * - This function performs an MPI_Allreduce operation to obtain the sum of integer values
 *   across all processes in the communicator.
 */
int Parallel::ISumAll(const int d, const int sum) {
  int sumout;

  if (npi[d] > 1)
    // comm_i[d].Allreduce(&sum,&sumout,1,MPI_INTEGER,MPI_SUM);
    MPI_Allreduce(&sum, &sumout, 1, MPI_INTEGER, MPI_SUM, comm_i[d]);
  else
    sumout = sum;
  return sumout;
}


/********************************
 *                              *
 *       Vector_SumAll          *
 *                              *
 ********************************/
/**
 * @brief Compute the sum of a vector of double values across all processes in the specified communicator.
 *
 * This function computes the element-wise sum of a vector of double values across all processes
 * in the specified communicator `d` and stores the result back in the `sum` array.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param n The number of elements in the `sum` vector.
 * @param sum An array of double values containing the input values to be summed. On output, it will contain the sum.
 *
 * @note
 * - This function performs an MPI_Allreduce operation to obtain the sum of double values
 *   across all processes in the communicator.
 * - The `sum` array will be modified to contain the element-wise sum after the function call.
 */
void Parallel::Vector_SumAll(const int d, const int n, double *sum) {
  double *sumout;
  if (npi[d] > 1)
  {
     sumout = new double[n];
     // comm_i[d].Allreduce(sum,sumout,n,MPI_DOUBLE_PRECISION,MPI_SUM);
     MPI_Allreduce(sum, sumout, n, MPI_DOUBLE_PRECISION, MPI_SUM, comm_i[d]);
     for (int i = 0; i < n; ++i)
        sum[i] = sumout[i];
     delete[] sumout;
  }
}


/********************************
 *                              *
 *       Vector_ISumAll         *
 *                              *
 ********************************/
/**
 * @brief Compute the element-wise sum of an integer vector across all processes in the specified communicator.
 *
 * This function computes the element-wise sum of an integer vector across all processes
 * in the specified communicator `d` and stores the result back in the `sum` array.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param n The number of elements in the `sum` vector.
 * @param sum An array of integer values containing the input values to be summed. On output, it will contain the sum.
 *
 * @note
 * - This function performs an MPI_Allreduce operation to obtain the sum of integer values
 *   across all processes in the communicator.
 * - The `sum` array will be modified to contain the element-wise sum after the function call.
 */
void Parallel::Vector_ISumAll(const int d, const int n, int *sum) {
  int *sumout;
  if (npi[d] > 1) {
    sumout = new int[n];
    // comm_i[d].Allreduce(sum,sumout,n,MPI_INTEGER,MPI_SUM);
    MPI_Allreduce(sum, sumout, n, MPI_INTEGER, MPI_SUM, comm_i[d]);
    for (int i = 0; i < n; ++i)
      sum[i] = sumout[i];
    delete[] sumout;
  }
}


/********************************
 *                              *
 *       Brdcst_Values          *
 *                              *
 ********************************/
/**
 * @brief Broadcast values from the root process to all processes in the specified communicator.
 *
 * This function broadcasts a set of double-precision values from the root process to all
 * processes in the specified communicator `d`. The `root` process is responsible for providing
 * the values to be broadcast, and all other processes receive the broadcasted values.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param root The rank of the root process responsible for providing the values.
 * @param n The number of double-precision values to broadcast.
 * @param sum An array of double values containing the values to be broadcasted from the root process.
 *            On non-root processes, this array will receive the broadcasted values.
 *
 * @note
 * - This function performs an MPI_Bcast operation to broadcast the values from the root process
 *   to all other processes in the communicator.
 * - On the root process, this function sends the values. On non-root processes, it receives the values.
 * - After the function call, all processes in the communicator will have the same values in the `sum` array.
 */
void Parallel::Brdcst_Values(const int d, const int root, const int n,
                             double *sum) {
  // if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_DOUBLE_PRECISION,root);
  if (npi[d] > 1)
    MPI_Bcast(sum, n, MPI_DOUBLE_PRECISION, root, comm_i[d]);
}


/********************************
 *                              *
 *       Brdcst_iValues         *
 *                              *
 ********************************/
/**
 * @brief Broadcast integer values from the root process to all processes in the specified communicator.
 *
 * This function broadcasts a set of integer values from the root process to all processes in the
 * specified communicator `d`. The `root` process is responsible for providing the values to be
 * broadcast, and all other processes receive the broadcasted values.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param root The rank of the root process responsible for providing the integer values.
 * @param n The number of integer values to broadcast.
 * @param sum An array of integer values containing the values to be broadcasted from the root process.
 *            On non-root processes, this array will receive the broadcasted values.
 *
 * @note
 * - This function performs an MPI_Bcast operation to broadcast the integer values from the root process
 *   to all other processes in the communicator.
 * - On the root process, this function sends the values. On non-root processes, it receives the values.
 * - After the function call, all processes in the communicator will have the same values in the `sum` array.
 */
void Parallel::Brdcst_iValues(const int d, const int root, const int n,
                              int *sum) {
  // if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_INTEGER,root);
  if (npi[d] > 1)
    MPI_Bcast(sum, n, MPI_INTEGER, root, comm_i[d]);
}


/********************************
 *                              *
 *       Brdcst_iValue          *
 *                              *
 ********************************/
/**
 * @brief Broadcast an integer value from the root process to all processes in the specified communicator.
 *
 * This function broadcasts a single integer value from the root process to all processes in the specified
 * communicator `d`. The `root` process is responsible for providing the integer value, and all other
 * processes receive the broadcasted value.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param root The rank of the root process responsible for providing the integer value.
 * @param sum A pointer to the integer value to be broadcasted. On non-root processes, this variable will
 *            receive the broadcasted value.
 *
 * @note
 * - This function performs an MPI_Bcast operation to broadcast the integer value from the root process
 *   to all other processes in the communicator.
 * - On the root process, this function sends the value. On non-root processes, it receives the value.
 * - After the function call, all processes in the communicator will have the same integer value in the `sum` variable.
 */
void Parallel::Brdcst_iValue(const int d, const int root, int *sum) {
  // if (npi[d]>1) comm_i[d].Bcast(sum,1,MPI_INTEGER,root);
  if (npi[d] > 1)
    MPI_Bcast(sum, 1, MPI_INTEGER, root, comm_i[d]);
}


/********************************
 *                              *
 *       Brdcst_cValues         *
 *                              *
 ********************************/
/**
 * @brief Broadcast a block of characters from the root process to all processes in the specified communicator.
 *
 * This function broadcasts a block of characters (bytes) from the root process to all processes in the specified
 * communicator `d`. The `root` process is responsible for providing the character data, and all other
 * processes receive the broadcasted data.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param root The rank of the root process responsible for providing the character data.
 * @param n The number of characters (bytes) to be broadcasted.
 * @param sum A pointer to the character data to be broadcasted. On non-root processes, this variable will
 *            receive the broadcasted data.
 *
 * @note
 * - This function performs an MPI_Bcast operation to broadcast the character data from the root process
 *   to all other processes in the communicator.
 * - On the root process, this function sends the data. On non-root processes, it receives the data.
 * - After the function call, all processes in the communicator will have the same character data in the `sum` variable.
 */
void Parallel::Brdcst_cValues(const int d, const int root, const int n,
                              void *sum) {
  // if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_CHAR,root);
  if (npi[d] > 1)
    MPI_Bcast(sum, n, MPI_CHAR, root, comm_i[d]);
}


/********************************
 *                              *
 *       Reduce_Values          *
 *                              *
 ********************************/
/**
 * @brief Perform a reduction operation on an array of double values across all processes in the specified communicator.
 *
 * This function performs a reduction operation on an array of double values across all processes in the specified
 * communicator `d`. The result of the reduction is stored in the `sumout` array on the root process.
 *
 * @param d The communicator index (0 for the primary communicator, 1 or 2 for derived communicators).
 * @param root The rank of the root process where the reduction result will be stored.
 * @param n The number of double values in the arrays `sumin` and `sumout`.
 * @param sumin A pointer to the input array of double values to be reduced.
 * @param sumout A pointer to the output array of double values that will store the reduction result. This array
 *               should be allocated and initialized on the root process and can be used to retrieve the result
 *               on all processes.
 *
 * @note
 * - This function performs an MPI_Reduce operation to combine and reduce the values in the `sumin` arrays across
 *   all processes into the `sumout` array on the root process.
 * - On the root process, the result of the reduction operation will be stored in the `sumout` array. On non-root
 *   processes, the `sumout` array remains unchanged.
 * - After the function call, the `sumout` array on the root process contains the reduction result, which can be
 *   accessed by all processes in the communicator.
 */
void Parallel::Reduce_Values(const int d, const int root, const int n,
                             double *sumin, double *sumout) {
   // if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_DOUBLE_PRECISION,root);
   if (npi[d] > 1)
      MPI_Reduce(sumin,sumout,n,MPI_DOUBLE_PRECISION,MPI_SUM,root,comm_i[d]);
   else
      std::memcpy(sumout,sumin,n*sizeof(double));
}


/********************************
 *                              *
 *       Parallel::dsend        *
 *                              *
 ********************************/
/**
 * @brief Perform a double-precision send operation in a specified dimension.
 *
 * This function initiates a double-precision send operation in a specified dimension `d`.
 * It allows the sending of double-precision floating-point data to a specified process.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] tag An integer tag for identifying the communication.
 * @param[in] procto The rank of the destination process receiving the data.
 * @param[in] n The number of double-precision values to send.
 * @param[in] sum A pointer to the memory location containing the double-precision values to be sent.
 *
 * @note
 * - Use this function to initiate double-precision sends in a specific communication dimension.
 * - The dimension `d` specifies the type of communication and is typically used for conditional checks.
 * - The `tag` parameter is used to identify the communication.
 * - `procto` specifies the rank of the destination process receiving the data.
 * - `n` indicates the number of double-precision values to send.
 * - The double-precision values to be sent are provided in the memory location pointed to by `sum`.
 * - If `npi[d] > 1`, the function initiates the send operation.
 * - The function does not block, and further completion checks may be required using related functions.
 */
void Parallel::dsend(const int d, const int tag, const int procto, const int n, double *sum) 
{
  // if (npi[d]>1) comm_i[d].Send(sum,n,MPI_DOUBLE_PRECISION,procto,tag);
   if (npi[d] > 1)
      MPI_Send(sum, n, MPI_DOUBLE_PRECISION, procto, tag, comm_i[d]);
}


/********************************
 *                              *
 *       Parallel::dreceive     *
 *                              *
 ********************************/
/**
 * @brief Perform a double-precision receive operation in a specified dimension.
 *
 * This function initiates a double-precision receive operation in a specified dimension `d`.
 * It allows the receiving of double-precision floating-point data from a specified process.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] tag An integer tag for identifying the communication.
 * @param[in] procfrom The rank of the source process sending the data.
 * @param[in] n The number of double-precision values to receive.
 * @param[out] sum A pointer to the memory location where the received double-precision values will be stored.
 *
 * @note
 * - Use this function to perform double-precision receives in a specific communication dimension.
 * - The dimension `d` specifies the type of communication and is typically used for conditional checks.
 * - The `tag` parameter is used to identify the communication.
 * - `procfrom` specifies the rank of the source process sending the data.
 * - `n` indicates the number of double-precision values to receive.
 * - The received double-precision values will be stored in the memory location pointed to by `sum`.
 * - If `npi[d] > 1`, the function initiates the receive operation.
 * - The function does not block, and further completion checks may be required using related functions.
 */
void Parallel::dreceive(const int d, const int tag, const int procfrom,
                        const int n, double *sum) {
  // MPI::Status status;
  // if (npi[d]>1) comm_i[d].Recv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag);
  MPI_Status status;
  if (npi[d] > 1)
    MPI_Recv(sum, n, MPI_DOUBLE_PRECISION, procfrom, tag, comm_i[d], &status);
}


/********************************
 *                              *
 *       Parallel::isend        *
 *                              *
 ********************************/
/**
 * @brief Perform an integer send operation in a specified dimension.
 *
 * This function initiates an integer send operation in a specified dimension `d`.
 * It allows the sending of integer data to a specified process.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] tag An integer tag for identifying the communication.
 * @param[in] procto The rank of the destination process receiving the data.
 * @param[in] n The number of integers to send.
 * @param[in] sum A pointer to the memory location containing the integers to send.
 *
 * @note
 * - Use this function to perform integer sends in a specific communication dimension.
 * - The dimension `d` specifies the type of communication and is typically used for conditional checks.
 * - The `tag` parameter is used to identify the communication.
 * - `procto` specifies the rank of the destination process receiving the data.
 * - `n` indicates the number of integers to send.
 * - The integer data to be sent is provided from the memory location pointed to by `sum`.
 * - If `npi[d] > 1`, the function initiates the send operation.
 * - The function does not block, and further completion checks may be required using related functions.
 */
void Parallel::isend(const int d, const int tag, const int procto, const int n,
                     int *sum) {
  // if (npi[d]>1) comm_i[d].Send(sum,n,MPI_INTEGER,procto,tag);
  if (npi[d] > 1)
    MPI_Send(sum, n, MPI_INTEGER, procto, tag, comm_i[d]);
}


/********************************
 *                              *
 *       Parallel::ireceive     *
 *                              *
 ********************************/
/**
 * @brief Perform a non-blocking integer receive operation in a specified dimension.
 *
 * This function initiates a non-blocking integer receive operation in a specified dimension `d`.
 * It allows the receiving of integer data from a specified process.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] tag An integer tag for identifying the communication.
 * @param[in] procfrom The rank of the source process sending the data.
 * @param[in] n The number of integers to receive.
 * @param[out] sum A pointer to the memory location where the received integers will be stored.
 *
 * @note
 * - Use this function to perform non-blocking integer receives in a specific communication dimension.
 * - The dimension `d` specifies the type of communication and is typically used for conditional checks.
 * - The `tag` parameter is used to identify the communication.
 * - `procfrom` specifies the rank of the source process sending the data.
 * - `n` indicates the number of integers to receive.
 * - The received integer data is stored in the memory location pointed to by `sum`.
 * - If `npi[d] > 1`, the function initiates the non-blocking receive operation.
 * - After initiating the receive, use a corresponding `Parallel::wait` or `Parallel::test` function
 *   to check for completion and retrieve the received data.
 */
void Parallel::ireceive(const int d, const int tag, const int procfrom,
                        const int n, int *sum) {
  // MPI::Status status;
  // if (npi[d]>1) comm_i[d].Recv(sum,n,MPI_INTEGER,procfrom,tag);
  MPI_Status status;
  if (npi[d] > 1)
    MPI_Recv(sum, n, MPI_INTEGER, procfrom, tag, comm_i[d], &status);
}


/********************************
 *                              *
 *       Parallel::astart       *
 *                              *
 ********************************/
/**
 * @brief Initialize MPI communication requests and statuses for asynchronous communication.
 *
 * This function initializes MPI communication requests and statuses for asynchronous communication.
 * It is typically used to prepare for non-blocking MPI communication operations.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] sz The size of arrays to allocate for communication requests and statuses.
 *
 * @note
 * - Use this function to prepare for asynchronous MPI communication operations.
 * - It initializes arrays for MPI requests and statuses based on the provided size.
 * - The dimension `d` specifies the type of communication, and it is typically used for conditional checks.
 * - Communication requests and statuses are allocated dynamically.
 * - After the communication is completed, deallocate these arrays to avoid memory leaks.
 */
void Parallel::astart(const int d, const int sz) 
{
   //this will need to be changed to  d>3 when k-points added
   if ((d > 2) ? true : (npi[d] > 1)) 
   {
      reqcnt[d]   = 0;
      request[d]  = new MPI_Request[sz];
      statuses[d] = new MPI_Status[sz];
   }
}


/********************************
 *                              *
 *       Parallel::awaitall     *
 *                              *
 ********************************/
/**
 * @brief Wait for the completion of MPI communication requests.
 *
 * This function waits for the completion of MPI communication requests, typically
 * used for synchronization between MPI processes.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 *
 * @note
 * - Use this function to wait for the completion of pending MPI communication requests
 *   in multi-process MPI applications.
 * - It checks for pending requests and waits for their completion using MPI_Waitall.
 * - If an error occurs during communication, an error message with the error code and task ID is
 *   printed to the standard output.
 * - After waiting for completion, the request count is reset to zero.
 */
void Parallel::awaitall(const int d) 
{
   //this will need to be changed to  d>3 when k-points added
   if ((d > 2) ? true : (npi[d] > 1)) 
   {
      if (reqcnt[d] > 0)
      {
         int ierr = MPI_Waitall(reqcnt[d], request[d], statuses[d]);
         if (ierr!=0) 
            std::cout << "Parallel::awaitall ierr=" << ierr << "taskid=" << taskidi[0] << std::endl;
      }
      reqcnt[d] = 0;
   }
}


/********************************
 *                              *
 *       Parallel::aend         *
 *                              *
 ********************************/
/**
 * @brief Wait for and complete asynchronous communications in a 2D communication space.
 *
 * This function waits for and completes asynchronous communications in a 2D communication space.
 * It ensures that all previously initiated non-blocking communication requests are finished.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 *
 * @note
 * - Use this function to wait for the completion of asynchronous communication requests.
 * - The dimension `d` specifies the type of communication and is typically used for conditional checks.
 * - This function should be called to ensure that asynchronous communications are properly completed.
 * - It waits for all previously initiated non-blocking communication requests to finish.
 * - After calling this function, it is safe to access and manipulate the communicated data.
 * - This function cleans up any allocated resources related to the communication.
 */
void Parallel::aend(const int d) 
{
   // MPI::Status status[reqcnt[d]];
   // request[d][0].Waitall(reqcnt[d],request[d],status);
   // if (npi[d]>1)

   //this will need to be changed to  d>3 when k-points added
   if ((d > 2) ? true : (npi[d] > 1)) 
   {
      // request[d][0].Waitall(reqcnt[d],request[d]);
      if (reqcnt[d] > 0)
         MPI_Waitall(reqcnt[d], request[d], statuses[d]);
      delete[] request[d];
      delete[] statuses[d];
      reqcnt[d] = 0;
   }
}


/********************************
 *                              *
 *       Parallel::adreceive    *
 *                              *
 ********************************/
/**
 * @brief Asynchronously receive a double-precision array from a specific processor in a 2D communication space.
 *
 * This function asynchronously receives a double-precision array from a specified processor
 * in a 2D communication space. It uses non-blocking MPI communication.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] tag An integer tag to identify the communication.
 * @param[in] procfrom The processor rank from which the data should be received.
 * @param[in] n The number of elements in the array to receive.
 * @param[out] sum Pointer to the array of double-precision values where the received data should be stored.
 *
 * @note
 * - Use this function to asynchronously receive data from a specific processor.
 * - The dimension `d` specifies the type of communication and is typically used for conditional checks.
 * - If `d > 2`, this function uses communicator 1 by default for communication.
 * - It uses non-blocking MPI communication, allowing the program to continue execution while the communication is in progress.
 * - The `tag` parameter is used to identify the communication and should match with the sending end's tag.
 * - After calling this function, you may need to use `Parallel::awaitall` to wait for the communication to complete.
 */
void Parallel::adreceive(const int d, const int tag, const int procfrom, const int n, double *sum) 
{
   // MPI::Status status;
   // if (npi[d]>1) request[d][reqcnt[d]++] =
   // comm_i[d].Irecv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag);

   //this will need to be changed to  d>3 when k-points added
   if ((d > 2) ? true : (npi[d] > 1))
     MPI_Irecv(sum, n, MPI_DOUBLE_PRECISION, procfrom, tag,
               comm_i[((d > 2) ? 1 : d)], &request[d][reqcnt[d]++]);
}


/********************************
 *                              *
 *       Parallel::adsend       *
 *                              *
 ********************************/
/**
 * @brief Asynchronously send a double-precision array to a specific processor in a 2D communication space.
 *
 * This function asynchronously sends a double-precision array to a specified processor
 * in a 2D communication space. It uses non-blocking MPI communication.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] tag An integer tag to identify the communication.
 * @param[in] procto The processor rank to which the data should be sent.
 * @param[in] n The number of elements in the array to send.
 * @param[in] sum Pointer to the array of double-precision values to be sent.
 *
 * @note
 * - Use this function to asynchronously send data to a specific processor.
 * - The dimension `d` specifies the type of communication, and it is typically used for conditional checks.
 * - If `d > 2`, this function uses communicator 1 by default for communication.
 * - It uses non-blocking MPI communication, allowing the program to continue execution while the communication is in progress.
 * - The `tag` parameter is used to identify the communication and should match with the receiving end's tag.
 * - After calling this function, you may need to use `Parallel::awaitall` to wait for the communication to complete.
 */
void Parallel::adsend(const int d, const int tag, const int procto, const int n, double *sum) 
{
   // if (npi[d]>1) request[d][reqcnt[d]++] =
   // comm_i[d].Isend(sum,n,MPI_DOUBLE_PRECISION,procto,tag);

   //this will need to be changed to  d>3 when k-points added
   if ((d > 2) ? true : (npi[d] > 1))
     MPI_Isend(sum, n, MPI_DOUBLE_PRECISION, procto, tag,
               comm_i[((d > 2) ? 1 : d)], &request[d][reqcnt[d]++]);
}


/********************************
 *                              *
 *       Parallel::a2dreceive   *
 *                              *
 ********************************/
/**
 * @brief Asynchronously receive a double-precision array from a specific processor in a 2D communication space.
 *
 * This function asynchronously receives a double-precision array from a specified processor
 * in a 2D communication space. It uses non-blocking MPI communication.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] tag An integer tag to identify the communication.
 * @param[in] procfrom The processor rank from which the data should be received.
 * @param[in] n The number of elements in the array to receive.
 * @param[out] sum Pointer to the array of double-precision values where the received data will be stored.
 *
 * @note
 * - Use this function to asynchronously receive data from a specific processor.
 * - The dimension `d` specifies the type of communication, and it is typically used for conditional checks.
 * - If `d > 2`, this function uses communicator 2 by default for communication.
 * - It uses non-blocking MPI communication, allowing the program to continue execution while the communication is in progress.
 * - The `tag` parameter is used to identify the communication and should match with the sending end's tag.
 * - After calling this function, you may need to use `Parallel::awaitall` to wait for the communication to complete.
 */
void Parallel::a2dreceive(const int d, const int tag, const int procfrom, const int n, double *sum)
{
   // MPI::Status status;
   // if (npi[d]>1) request[d][reqcnt[d]++] =
   // comm_i[d].Irecv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag);

   //this will need to be changed to  d>3 when k-points added
   if ((d > 2) ? true : (npi[d] > 1))
     MPI_Irecv(sum, n, MPI_DOUBLE_PRECISION, procfrom, tag,
               comm_i[((d > 2) ? 2 : d)], &request[d][reqcnt[d]++]);
}


/********************************
 *                              *
 *       Parallel::a2dsend      *
 *                              *
 ********************************/
/**
 * @brief Asynchronously send a double-precision array to a specific processor in a 2D communication space.
 *
 * This function sends a double-precision array to a specified processor asynchronously
 * in a 2D communication space. It uses non-blocking MPI communication.
 *
 * @param[in] d An integer indicating the dimension of communication (e.g., dimension of the communication space).
 * @param[in] tag An integer tag to identify the communication.
 * @param[in] procto The processor rank to which the data should be sent.
 * @param[in] n The number of elements in the array to send.
 * @param[in] sum Pointer to the array of double-precision values to send.
 *
 * @note
 * - Use this function to asynchronously send data to a specific processor.
 * - The dimension `d` specifies the type of communication, and it is typically used for conditional checks.
 * - If `d > 2`, this function uses communicator 2 by default for communication.
 * - It uses non-blocking MPI communication, allowing the program to continue execution while the communication is in progress.
 * - The `tag` parameter is used to identify the communication and should match with the receiving end's tag.
 * - After calling this function, you may need to use `Parallel::awaitall` to wait for the communication to complete.
 */
void Parallel::a2dsend(const int d, const int tag, const int procto, const int n, double *sum)
{
   // if (npi[d]>1) request[d][reqcnt[d]++] =
   // comm_i[d].Isend(sum,n,MPI_DOUBLE_PRECISION,procto,tag);

   //this will need to be changed to  d>3 when k-points added
   if ((d > 2) ? true : (npi[d] > 1))
     MPI_Isend(sum, n, MPI_DOUBLE_PRECISION, procto, tag,
               comm_i[((d > 2) ? 2 : d)], &request[d][reqcnt[d]++]);
}


} // namespace pwdft
