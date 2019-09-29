#pragma once
#include <algorithm>
#include <cassert>

// M mxn
// N mmxnn
// res mxnn
template <typename T>
void matmult(T* M, int m, int n, T* N, int mm, int nn, T* res) {
  (void)mm;
  assert(n == mm);
  for (int i = 0; i < m; i++) {     // row
    for (int j = 0; j < nn; j++) {  // col
      res[i * nn + j] = 0;
      for (int k = 0; k < n; k++) {
        res[i * nn + j] += M[i * n + k] * N[k * nn + j];
      }
    }
  }
}

// M mxn
// N mmxnn
// res mxnn
template <typename T>
void matmult2DwTranspose(T** M, int m, int n, const bool transposeM, T** N,
                         int mm, int nn, const bool transposeN, T** res) {
  if (transposeM)
    std::swap(m, n);
  if (transposeN)
    std::swap(mm, nn);
  assert(n == mm);
  for (int i = 0; i < m; i++) {     // row
    for (int j = 0; j < nn; j++) {  // col
      res[i][j] = 0;
      for (int k = 0; k < n; k++) {
        if (!transposeM && !transposeN)
          res[i][j] += M[i][k] * N[k][j];
        else if (!transposeM && transposeN)
          res[i][j] += M[i][k] * N[j][k];
        else if (transposeM && !transposeN)
          res[i][j] += M[k][i] * N[k][j];
        else if (transposeM && transposeN)
          res[i][j] += M[k][i] * N[j][k];
      }
    }
  }
}
