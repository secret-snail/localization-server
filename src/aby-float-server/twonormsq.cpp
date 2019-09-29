#include <stdio.h>
#include <iostream>
#include <string>
#include "jlog.h"

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

using namespace std;

/*
 An implementation of SVD from Numerical Recipes in C and Mike Erhdmann's
 lectures
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const uint32_t bitlen = 32;

// M mxn (M is secret, m n are public)
// N mmxnn
// res mxnn
share* BuildTwoNormSqCircuit(share* vect[], int sz, BooleanCircuit* c) {
  float zero = 0.0;
  share* sum = c->PutCONSGate((uint32_t*)&zero, bitlen);
  for (int i = 0; i < sz; i++) {
    share* temp = c->PutFPGate(vect[i], vect[i], MUL, bitlen, 1, no_status);
    share* oldsum = sum;
    sum = c->PutFPGate(sum, temp, ADD, bitlen, 1, no_status);
    delete temp;
    delete oldsum;
  }
  return sum;
}
