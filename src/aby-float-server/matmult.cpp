#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

using namespace std;

const uint32_t bitlen = 32;

// M mxn (M is secret, m n are public)
// N mmxnn
// res mxnn
void BuildMatmultCircuit(share* M[], int m, int n, share* N[], int mm, int nn,
                         share* res[], BooleanCircuit* c) {
  (void)mm;
  assert(n == mm);
  for (int i = 0; i < m; i++) {     // row
    for (int j = 0; j < nn; j++) {  // col
      float z = 0;
      res[i * nn + j] = c->PutCONSGate((uint32_t*)&z, 32);
      for (int k = 0; k < n; k++) {
        share* temp = c->PutFPGate(M[i * n + k], N[k * nn + j], MUL, bitlen, 1,
                                   no_status);
        share* temp2 = res[i * nn + j];
        res[i * nn + j] =
            c->PutFPGate(res[i * nn + j], temp, ADD, bitlen, 1, no_status);
        delete temp;
        delete temp2;
      }
    }
  }
}

// M mxn (M is secret, m n are public)
// N mmxnn
// res mxnn
void BuildMatmult2DwTransposeCircuit(share** M[], int m, int n,
                                     const bool transposeM, share** N[], int mm,
                                     int nn, const bool transposeN,
                                     share** res[], BooleanCircuit* c) {
  if (transposeM)
    std::swap(m, n);
  if (transposeN)
    std::swap(mm, nn);
  assert(n == mm);

  float z = 0;
  for (int i = 0; i < m; i++) {     // row
    for (int j = 0; j < nn; j++) {  // col
      res[i][j] = c->PutCONSGate((uint32_t*)&z, 32);

      for (int k = 0; k < n; k++) {
        share* prod;
        if (!transposeM && !transposeN)
          prod = c->PutFPGate(M[i][k], N[k][j], MUL, bitlen, 1, no_status);
        else if (!transposeM && transposeN)
          prod = c->PutFPGate(M[i][k], N[j][k], MUL, bitlen, 1, no_status);
        else if (transposeM && !transposeN)
          prod = c->PutFPGate(M[k][i], N[k][j], MUL, bitlen, 1, no_status);
        else if (transposeM && transposeN)
          prod = c->PutFPGate(M[k][i], N[j][k], MUL, bitlen, 1, no_status);

        share* oldres = res[i][j];
        res[i][j] = c->PutFPGate(res[i][j], prod, ADD, bitlen, 1, no_status);
        delete prod;
        delete oldres;
      }
    }
  }
}
