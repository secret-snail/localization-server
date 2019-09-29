#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

using namespace emp;
using namespace std;

void BuildMatmultCircuit(Float M[], int m, int n, Float N[], int mm, int nn,
                         Float res[]) {
  (void)mm;
  assert(n == mm);
  for (int i = 0; i < m; i++) {     // row
    for (int j = 0; j < nn; j++) {  // col
      // res[i*nn +j] = Float(0.0, ALICE);
      new (&res[i * nn + j]) Float(0.0, PUBLIC);
      for (int k = 0; k < n; k++) {
        Float temp = M[i * n + k] * N[k * nn + j];
        res[i * nn + j] = res[i * nn + j] + temp;
      }
    }
  }
}

void BuildMatmult2DwTransposeCircuit(Float** M, int m, int n,
                                     const bool transposeM, Float** N, int mm,
                                     int nn, const bool transposeN,
                                     Float** res) {
  if (transposeM)
    std::swap(m, n);
  if (transposeN)
    std::swap(mm, nn);
  assert(n == mm);

  for (int i = 0; i < m; i++) {     // row
    for (int j = 0; j < nn; j++) {  // col
      new (&res[i][j]) Float(0.0, PUBLIC);
      for (int k = 0; k < n; k++) {
        Float temp = Float(0.0, PUBLIC);
        if (!transposeM && !transposeN)
          temp = M[i][k] * N[k][j];
        else if (!transposeM && transposeN)
          temp = M[i][k] * N[j][k];
        else if (transposeM && !transposeN)
          temp = M[k][i] * N[k][j];
        else if (transposeM && transposeN)
          temp = M[k][i] * N[j][k];
        res[i][j] = res[i][j] + temp;
      }
    }
  }
}
