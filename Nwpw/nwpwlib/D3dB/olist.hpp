#ifndef _OLIST_HPP_
#define _OLIST_HPP_

#pragma once

/* olist.hpp -
   Author - Eric Bylaska

*/

#include <cstdlib>
#include <iostream>

namespace pwdft {

class OList {

  int max_index;
  int *list;

public:
  /* Constructor */
  OList(const int size) {
    list = new int[size];
    max_index = 0;
    for (int i = 0; i < size; ++i)
      list[i] = 0;
  }

  /* Destructor */
  ~OList() { delete[] list; }

  /* insertion operator */
  void insert(const int item) {
    int ii, j;
    ++max_index;
    ii = 0;
    while ((list[ii] < item) && (ii < (max_index - 1)))
      ++ii;
    for (j = (max_index - 1); j > ii; --j)
      list[j] = list[j - 1];
    list[ii] = item;
  }

  int index(const int item) {
    int ii = 0;
    while (list[ii] < item)
      ++ii;
    return ii;
  }

  void print() {
    int i;
    std::cout << max_index << ": ";
    for (i = 0; i < (max_index); ++i)
      std::cout << list[i] << " ";
    std::cout << "\n";
  }
};

// extern void    create_olist(OList_Type*, const int);
// extern void    insert_olist(OList_Type*, const int);
// extern	int	index_olist();
// extern	void	destroy_olist();
// extern	void	print_olist();

} // namespace pwdft

#endif
