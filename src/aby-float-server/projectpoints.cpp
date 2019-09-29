#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

#include <matmult.h>
#include <rodrigues.h>

using namespace std;

// Note: All 2D matrix inputs passed as single dimension array
// P=4xnumPoints [x1, x2, ...;  y1, y2, ...;  z1, z2, ...; 1, 1, ... ]
// x=6x1 [ rotation angles; translation] (column)
// K=3x3 camera intrinsic
// result=3xnumPoints preallocated, only first two rows are x,y.
//      last row needed for intermediate calculation (homogeneous)
//      [ u1, u2, ...;  v1, v2, ...;  1, 1, ... ]
void BuildProjectPointsCircuit(share* _P[], share* _x[], share* _K[],
                               share* _res[], int numPoints, BooleanCircuit* c,
                               bool skipOneGate) {

  uint32_t bitlen = 32;

  float one = 1;
  share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);

  share* _R[12];

  BuildRodriguesCircuit(_x, _R, (BooleanCircuit*)c);

  // copy translation into last column of R
  _R[3] = _x[3];
  _R[7] = _x[4];
  _R[11] = _x[5];

  // res = K*R*P;
  share* scratch[3 * 4];
  BuildMatmultCircuit(_K, 3, 3, _R, 3, 4, scratch, c);
  BuildMatmultCircuit(scratch, 3, 4, _P, 4, numPoints, _res, c);

  // homogeneous coords to x,y
  for (int i = 0; i < numPoints; i++) {
    _res[i] = c->PutFPGate(_res[i], _res[2 * numPoints + i], DIV, bitlen, 1,
                           no_status);
    _res[numPoints + i] =
        c->PutFPGate(_res[numPoints + i], _res[2 * numPoints + i], DIV, bitlen,
                     1, no_status);
    if (!skipOneGate) {
      _res[2 * numPoints + i] = one_gate;  // can be excluded for efficiency
    }
  }

  for (int i = 0; i < 12; i++) {
    if (i % 4 != 3)
      delete _R[i];  // skip 3,7,11 (don't delete _x's)
  }
}
