#ifndef _MAPPING3_HPP_
#define _MAPPING3_HPP_

#pragma once

/* mapping.h
   Author - Eric Bylaska

        this class is used defining 3d parallel maps
*/

namespace pwdft {

class Mapping3 {

  // int *qmap[3],*pmap[3],*kmap;
  int *qmap[6], *pmap[6], *kmap;

public:
  int np, taskid;
  int maptype, n2ft3d, nfft3d, n2ft3d_map, nfft3d_map, nrft3d, nrft3d_map;
  int nx, ny, nz;
  int nq, nq1, nq2, nq3, nqr, nqr1, nqr2, nqr3;

  /* Constructors */
  Mapping3();
  Mapping3(const int, const int, const int, const int, const int, const int);

  /* copy constructor */
  // Mapping3(const Mapping3&);

  /* destructor */
  ~Mapping3();

  // complex mapping
  int ijktoindex(const int i, const int j, const int k) 
  {
     if (maptype == 1)
        return (i + j * (nx/2+1) + qmap[0][k]*(nx/2+1)*ny);
     else
        return (k + qmap[2][i + j * (nx / 2 + 1)] * nz);
  }

  int ijktoindex1(const int i, const int j, const int k) {
    if (maptype == 1)
      return (i + k * (nx / 2 + 1) + qmap[0][j] * (nx / 2 + 1) * nz);
    else
      return (j + qmap[1][k + i * nz] * ny);
  }
  int ijktoindex2(const int i, const int j, const int k) {
    if (maptype == 1)
      return (i + k * (nx + 2) + qmap[0][j] * (nx + 2) * ny);
    else
      return (i + qmap[0][j + k * ny] * (nx + 2));
  }
  int ijktoindex2t(const int i, const int j, const int k) {
    if (maptype == 1)
      return (i + k * (nx + 2) + qmap[0][j] * (nx / 2 + 1) * ny);
    else
      return (i + qmap[0][j + k * ny] * (nx / 2 + 1));
  }

  int ijktop(const int i, const int j, const int k) 
  {
     if (maptype == 1)
        return (pmap[0][k]);
     else
        return (pmap[2][i + j * (nx/2+1)]);
  }

  int ijktop1(const int i, const int j, const int k) {
    if (maptype == 1)
      return (pmap[0][j]);
    else
      return (pmap[1][k + i * nz]);
  }
  int ijktop2(const int i, const int j, const int k) {
    if (maptype == 1)
      return (pmap[0][j]);
    else
      return (pmap[0][j + k * ny]);
  }

  int ijktoq(const int i, const int j, const int k) {
    if (maptype == 1)
      return (qmap[0][k]);
    else
      return (qmap[2][i + j * (nx / 2 + 1)]);
  }
  int ijktoq1(const int i, const int j, const int k) {
    if (maptype == 1)
      return (qmap[0][j]);
    else
      return (qmap[1][k + i * nz]);
  }
  int ijktoq2(const int i, const int j, const int k) {
    if (maptype == 1)
      return (qmap[0][j]);
    else
      return (qmap[0][j + k * ny]);
  }

  // real mapping
  int ijkrtoindex(const int i, const int j, const int k) {
    if (maptype == 1)
      return (i + j * (nx + 2) + qmap[0][k] * (nx + 2) * ny);
    else
      return (k + qmap[5][i + j * (nx + 2)] * nz);
  }
  int ijkrtoindex1(const int i, const int j, const int k) {
    if (maptype == 1)
      return (i + k * (nx + 2) + qmap[0][j] * (nx + 2) * nz);
    else
      return (j + qmap[4][k + i * nz] * ny);
  }
  int ijkrtoindex2(const int i, const int j, const int k) {
    if (maptype == 1)
      return (i + k * (nx + 2) + qmap[0][j] * (nx + 2) * ny);
    else
      return (i + qmap[3][j + k * ny] * (nx + 2));
  }
  int ijkrtoindex2t(const int i, const int j, const int k) {
    if (maptype == 1)
      return (i + k * (nx + 2) + qmap[0][j] * (nx + 2) * ny);
    else
      return (i + qmap[3][j + k * ny] * (nx + 2));
  }

  int ijkrtop(const int i, const int j, const int k) {
    if (maptype == 1)
      return (pmap[0][k]);
    else
      return (pmap[5][i + j * (nx + 2)]);
  }
  int ijkrtop1(const int i, const int j, const int k) {
    if (maptype == 1)
      return (pmap[0][j]);
    else
      return (pmap[4][k + i * nz]);
  }
  int ijkrtop2(const int i, const int j, const int k) {
    if (maptype == 1)
      return (pmap[0][j]);
    else
      return (pmap[3][j + k * ny]);
  }

  int ijkrtoq(const int i, const int j, const int k) {
    if (maptype == 1)
      return (qmap[0][k]);
    else
      return (qmap[5][i + j * (nx + 2)]);
  }
  int ijkrtoq1(const int i, const int j, const int k) {
    if (maptype == 1)
      return (qmap[0][j]);
    else
      return (qmap[4][k + i * nz]);
  }
  int ijkrtoq2(const int i, const int j, const int k) {
    if (maptype == 1)
      return (qmap[0][j]);
    else
      return (qmap[3][j + k * ny]);
  }
};

} // namespace pwdft

#endif
