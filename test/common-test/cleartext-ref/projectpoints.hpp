#pragma once
#include <printutil.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include <cleartext-ref/matmult.hpp>
#include <cleartext-ref/rodrigues.hpp>

using namespace std;

// Note: All 2D matrix inputs passed as single dimension array
// P=4xnumPoints [x1, x2, ...;  y1, y2, ...;  z1, z2, ...; 1, 1, ... ]
// x=6x1 [ rotation angles; translation] (column)
// K=3x3 camera intrinsic
// result=3xnumPoints preallocated, only first two rows are x,y.
//      last row needed for intermediate calculation (homogeneous)
//      [ u1, u2, ...;  v1, v2, ...;  1, 1, ... ]
template <typename T>
void projectPoints(T* _P, T* _x, T* _K, T* _res, int numPoints) {
  // rodreigues operates on a 3x4 matrix
  T _R[12];  // 3x4
  rodrigues(_x, _R);

  // copy translation into last column of R
  _R[3] = _x[3];
  _R[7] = _x[4];
  _R[11] = _x[5];

  // printVector("R", _R, 12);
  // printVector("K", _K, 9);

  // res = K*R*P;
  T* scratch = new T[3 * 4];
  matmult(_K, 3, 3, _R, 3, 4, scratch);
  // printVector("scratch", scratch, 3*4);
  matmult(scratch, 3, 4, _P, 4, numPoints, _res);
  // printVector("_res", _res, 3*numPoints);
  delete[] scratch;

  // homogeneous coords to x,y
  for (int i = 0; i < numPoints; i++) {
    _res[i] = _res[i] / _res[2 * numPoints + i];
    _res[numPoints + i] = _res[numPoints + i] / _res[2 * numPoints + i];
    _res[2 * numPoints + i] = 1;
  }
  // printVector("_resHomog", _res, 3*numPoints);
}
