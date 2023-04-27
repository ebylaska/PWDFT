/* hilbert3.C -
   Author - Eric Bylaska

   This file contains 3d hilbert mapping routines

*/

#include "olist.hpp"
#include <cmath>

static int length;

#define front_bottom_left 0
#define front_bottom_right 1
#define front_top_left 2
#define front_top_right 3
#define back_bottom_left 4
#define back_bottom_right 5
#define back_top_left 6
#define back_top_right 7

#define right 0
#define left 1
#define up 2
#define down 3

#define forward 4
#define back 5

int hilbert3d(int *p0, int *p1, int *p2, int *p3) {
  int d1[3], d2[3], d3[3];

  d1[0] = p1[0] - p0[0];
  d1[1] = p1[1] - p0[1];
  d1[2] = p1[2] - p0[2];
  d2[0] = p2[0] - p0[0];
  d2[1] = p2[1] - p0[1];
  d2[2] = p2[2] - p0[2];
  d3[0] = p3[0] - p0[0];
  d3[1] = p3[1] - p0[1];
  d3[2] = p3[2] - p0[2];

  if (d1[0] == 1) {
    int indx = p0[0] + p0[1] * count + p0[2] * count * count;
    map3[indx] = length;
    ++length;
  } else {
    q0 = p0;
    q1 = p0 + d2 / 2;
    q2 = p0 + d3 / 2;
    q3 = p0 + d1 / 2;
    hilbert3d(q0, q1, q2, q3);

    q0 = p0 + d2 / 2;
    q1 = p0 + d2 / 2 + d3 / 2;
    q2 = p0 + d2;
    q3 = p0 + d1 / 2 + d2 / 2;
    hilbert3d(q0, q1, q2, q3);

    q0 = p0 + d2 / 2 + d3 / 2;
    q1 = p0 + d2 / 2 + d3;
    q2 = p0 + d2 + d3 / 2;
    q3 = p0 + d1 / 2 + d2 / 2 + d3 / 2;
    hilbert3d(q0, q1, q2, q3);
  }

  length = 1;
  for (count = 0; count < (high - level); ++count)
    length = length * 8;

  if (level == 0) {
    direction = right;
    *start = 0;
  } else {
    parent_direction =
        hilbert3_dir(parent(i), parent(j), parent(k), level - 1, high, start);
    crnr = corner(i, j, k);

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

    if (parent_direction == forward) {
    }
    if (parent_direction == back) {
    }
  }
  return direction;
}

int hilbert3d(int i, int j, int k, int level) {
  int start, direction;
  direction = hilbert3_dir(i, j, k, level, level, &start);
  return start;
}

void hilbert3d_map(const int sizex, const int sizey, const int sizez,
                   int *map) {
  int i, j, size;
  int ii, jj, iii, jjj;
  int level, count;
  double dx, dy, x2, y2, dx2, dy2;
  int *map3, *mapii, *mapjj, *mapkk;

  size = sizex;
  if (sizey > size)
    size = sizey;
  if (sizez > size)
    size = sizez;

  /* get the level of map */
  count = 1;
  level = 0;
  while (count < size) {
    ++level;
    count = count * 2;
  }

  map3 = new int[count * count * count];
  mapii = new int[count * count * count];
  mapjj = new int[count * count * count];
  mapkk = new int[count * count * count];
  OList olist3(count * count * count);

  for (kk = 0; kk < count; ++kk)
    for (jj = 0; jj < count; ++jj)
      for (ii = 0; ii < count; ++ii) {

        map3[ii + jj * count + kk * count * count] =
            hilbert3d(ii, jj, kk, level);
        olist3.insert(map3[ii + jj * count + kk * count * count]);
      }
  for (kk = 0; kk < count; ++kk)
    for (jj = 0; jj < count; ++jj)
      for (ii = 0; ii < count; ++ii)
        map3[ii + jj * count + kk * count * count] =
            olist3.index(map3[ii + jj * count + kk * count * count]);

  for (kk = 0; kk < count; ++kk)
    for (jj = 0; jj < count; ++jj)
      for (ii = 0; ii < count; ++ii) {
        mapii[map3[ii + jj * count + kk * count * count]] = ii;
        mapjj[map3[ii + jj * count + kk * count * count]] = jj;
        mapkk[map3[ii + jj * count + kk * count * count]] = kk;
      }

  dx2 = 1.0 / ((double)count);
  dy2 = 1.0 / ((double)count);
  dz2 = 1.0 / ((double)count);
  dx = 1.0 / ((double)sizex);
  dy = 1.0 / ((double)sizey);
  dz = 1.0 / ((double)sizez);

  for (j = 0; j < sizex * sizey * sizez; ++j)
    map[j] = -9;

  iii = 0;
  for (jjj = 0; jjj < count * count * count; ++jjj) {
    ii = mapii[jjj];
    jj = mapjj[jjj];
    kk = mapkk[jjj];

    x2 = dx2 * (ii + 0.5);
    y2 = dy2 * (jj + 0.5);
    z2 = dz2 * (kk + 0.5);
    i = ((int)rint((x2 / dx) - 0.5));
    j = ((int)rint((y2 / dy) - 0.5));
    k = ((int)rint((z2 / dz) - 0.5));

    if (map[i + j * sizex + k * sizex * sizey] == -9) {
      map[i + j * sizex + k * sizex * sizey] = iii;
      ++iii;
    }
  }
  delete[] mapii;
  delete[] mapjj;
  delete[] mapkk;
  delete[] map3;
}
