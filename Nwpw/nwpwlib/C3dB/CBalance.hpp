#ifndef _CBALANCE_HPP_
#define _CBALANCE_HPP_

#pragma once

/* CBalance.hpp
   Author - Eric Bylaska

*/

#include "Parallel.hpp"
namespace pwdft {

#define MAXSIZE_MAX 2

/**
 * @class Balance
 * @brief Class for load balancing in parallel computing.
 *
 * The `Balance` class provides load balancing capabilities for parallel computations within
 * the parallel computing framework. It includes methods for load balancing, unbalancing, and
 * other related operations.
 *
 * @author Eric Bylaska
 * @date 2023-09-20
 *
 * @note Detailed function descriptions and implementation can be found in corresponding source files.
 *
 * @see Parallel.hpp
 */
class CBalance {

  Parallel *parall;

  int maxsize;
  int *proc_to_list[MAXSIZE_MAX];
  int *proc_from_list[MAXSIZE_MAX];
  int *packet_size_list[MAXSIZE_MAX];
  int *indx_start_list[MAXSIZE_MAX];

  int npacket_list[MAXSIZE_MAX];
  int receiver_list[MAXSIZE_MAX];
  int sender_list[MAXSIZE_MAX];

public:
  /* Constructors */
  CBalance(Parallel *, const int, const int *, int *);

  /* destructor */
  ~CBalance();

  void c_unbalance(const int, double *);
  void c_balance(const int, double *);
  void r_unbalance(const int, double *);
  void r_balance(const int, double *);
  void i_balance(const int, int *);

  void c_unbalance_start(const int, double *, const int, const int);
  void c_unbalance_end(const int, double *, const int);

  void c_balance_start(const int, double *, const int, const int);
  void c_balance_end(const int, double *, const int);
};

} // namespace pwdft

#endif
