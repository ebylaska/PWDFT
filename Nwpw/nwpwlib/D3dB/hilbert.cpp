/* hilbert.cpp -
$Id: hilbert.c,v 1.6 2005-10-07 19:34:21 bylaska Exp $
   Author - Eric Bylaska
   This file contains 2d hilbert mapping routines
*/

#include "olist.hpp"
#include <cmath>

namespace pwdft {

#define bottom_left 0
#define bottom_right 1
#define top_left 2
#define top_right 3

#define right 0
#define left 1
#define up 2
#define down 3


/***********************************
 *                                 *
 *              parent             *
 *                                 *
 ***********************************/
/**
 * @brief Calculate the parent node of a binary tree given a node index.
 *
 * This function calculates and returns the index of the parent node of a node
 * in a binary tree. It divides the given node index by 2 to find the parent
 * node index.
 *
 * @param i  The index of the node for which the parent node is calculated.
 * @return   The index of the parent node.
 */
int parent(int i) { return (i / 2); }


/***********************************
 *                                 *
 *              corner             *
 *                                 *
 ***********************************/
/**
 * @brief Determine the corner type of a cell or grid element.
 *
 * This function calculates and returns an integer representing the corner type
 * of a cell or grid element based on its row and column indices.
 *
 * @param i  The row index of the cell or grid element.
 * @param j  The column index of the cell or grid element.
 * @return   An integer representing the corner type.
 */
int corner(int i, int j) { return (2 * (j % 2) + (i % 2)); }


/***********************************
 *                                 *
 *          hilbert_dir            *
 *                                 *
 ***********************************/
/**
 * @brief Calculate the direction of a point in a Hilbert curve.
 *
 * This function calculates the direction of a point in a Hilbert curve based on
 * its coordinates (i, j), the current recursion level (level), and the highest
 * level (high). It also updates the starting position (start) of the curve segment.
 *
 * @param i       The x-coordinate of the point.
 * @param j       The y-coordinate of the point.
 * @param level   The current recursion level.
 * @param high    The highest recursion level.
 * @param start   Pointer to the starting position, which is updated by the function.
 * @return        The direction of the point in the Hilbert curve (e.g., up, down, left, right).
 */
int hilbert_dir(int i, int j, int level, int high, int *start) {
  int direction, parent_direction, crnr, length, count;

  length = 1;
  for (count = 0; count < (high - level); ++count)
    length = length * 4;

  if (level == 0) {
    direction = right;
    *start = 0;
  } else {
    parent_direction =
        hilbert_dir(parent(i), parent(j), level - 1, high, start);
    crnr = corner(i, j);

    if (parent_direction == right) {
      if (crnr == bottom_left) {
        direction = up;
        *start = *start + 0 * length;
      }
      if (crnr == bottom_right) {
        direction = down;
        *start = *start + 3 * length;
      }
      if (crnr == top_left) {
        direction = right;
        *start = *start + 1 * length;
      }
      if (crnr == top_right) {
        direction = right;
        *start = *start + 2 * length;
      }
    }

    if (parent_direction == left) {
      if (crnr == bottom_left) {
        direction = left;
        *start = *start + 2 * length;
      }
      if (crnr == bottom_right) {
        direction = left;
        *start = *start + 1 * length;
      }
      if (crnr == top_left) {
        direction = up;
        *start = *start + 3 * length;
      }
      if (crnr == top_right) {
        direction = down;
        *start = *start + 0 * length;
      }
    }

    if (parent_direction == up) {
      if (crnr == bottom_left) {
        direction = right;
        *start = *start + 0 * length;
      }
      if (crnr == bottom_right) {
        direction = up;
        *start = *start + 1 * length;
      }
      if (crnr == top_left) {
        direction = left;
        *start = *start + 3 * length;
      }
      if (crnr == top_right) {
        direction = up;
        *start = *start + 2 * length;
      }
    }

    if (parent_direction == down) {
      if (crnr == bottom_left) {
        direction = down;
        *start = *start + 2 * length;
      }
      if (crnr == bottom_right) {
        direction = right;
        *start = *start + 3 * length;
      }
      if (crnr == top_left) {
        direction = down;
        *start = *start + 1 * length;
      }
      if (crnr == top_right) {
        direction = left;
        *start = *start + 0 * length;
      }
    }
  }
  return direction;
}


/***********************************
 *                                 *
 *            hilbert2d            *
 *                                 *
 ***********************************/
/**
 * @brief Compute the starting position of a Hilbert curve traversal in 2D.
 *
 * This function calculates the starting position of a Hilbert curve traversal
 * in a 2D grid based on the given coordinates (`i`, `j`) and the level of the
 * Hilbert curve.
 *
 * @param i      The x-coordinate in the grid.
 * @param j      The y-coordinate in the grid.
 * @param level  The level of the Hilbert curve.
 * @return       The starting position of the Hilbert curve traversal.
 */
int hilbert2d(int i, int j, int level) {
  int start, direction;
  direction = hilbert_dir(i, j, level, level, &start);
  return start;
}


/***********************************
 *                                 *
 *         hilbert2d_map           *
 *                                 *
 ***********************************/
/**
 * @brief Generate a Hilbert curve-based mapping for a 2D grid.
 *
 * This function generates a Hilbert curve-based mapping for a 2D grid of size
 * `sizex` by `sizey`. The resulting mapping is stored in the `map` array, where
 * each element corresponds to a cell in the original grid, and it represents
 * the order of the cell in the Hilbert curve traversal.
 *
 * @param sizex The size of the grid in the x-direction.
 * @param sizey The size of the grid in the y-direction.
 * @param map   An array to store the Hilbert curve-based mapping.
 */
void hilbert2d_map(const int sizex, const int sizey, int *map) {
  int i, j, size;
  int ii, jj, iii, jjj;
  int level, count;
  double dx, dy, x2, y2, dx2, dy2;
  int *map2, *mapii, *mapjj;

  size = sizex;
  if (sizey > size)
    size = sizey;

  /* get the level of map */
  count = 1;
  level = 0;
  while (count < size) {
    ++level;
    count = count * 2;
  }

  map2 = new int[count * count];
  mapii = new int[count * count];
  mapjj = new int[count * count];
  OList olist2(count * count);

  for (jj = 0; jj < count; ++jj)
    for (ii = 0; ii < count; ++ii) {

      map2[ii + jj * count] = hilbert2d(ii, jj, level);
      olist2.insert(map2[ii + jj * count]);
    }
  for (jj = 0; jj < count; ++jj)
    for (ii = 0; ii < count; ++ii)
      map2[ii + jj * count] = olist2.index(map2[ii + jj * count]);

  for (jj = 0; jj < count; ++jj)
    for (ii = 0; ii < count; ++ii) {
      mapii[map2[ii + jj * count]] = ii;
      mapjj[map2[ii + jj * count]] = jj;
    }

  dx2 = 1.0 / ((double)count);
  dy2 = 1.0 / ((double)count);
  dx = 1.0 / ((double)sizex);
  dy = 1.0 / ((double)sizey);

  for (j = 0; j < sizex * sizey; ++j)
    map[j] = -9;

  iii = 0;
  for (jjj = 0; jjj < count * count; ++jjj) {
    ii = mapii[jjj];
    jj = mapjj[jjj];

    x2 = dx2 * (ii + 0.5);
    y2 = dy2 * (jj + 0.5);
    i = ((int)rint((x2 / dx) - 0.5));
    j = ((int)rint((y2 / dy) - 0.5));

    if (map[i + j * sizex] == -9) {
      map[i + j * sizex] = iii;
      ++iii;
    }
  }
  delete[] mapii;
  delete[] mapjj;
  delete[] map2;
}
} // namespace pwdft
