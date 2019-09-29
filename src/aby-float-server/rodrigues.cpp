#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

#include <rodrigues.h>
#include <trigfuncs.h>

using namespace std;

// r=3x1, R=3x4 (only 3x3 is used here)
// values are secret, sizes are public
void BuildRodriguesCircuit(share* r[], share* R[], BooleanCircuit* c) {

  uint32_t bitlen = 32;
  float one = 1;
  float zero = 0;
  share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);
  share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);

  // calculate theta
  share* rzerosq = c->PutFPGate(r[0], r[0], MUL, bitlen, 1, no_status);
  share* ronesq = c->PutFPGate(r[1], r[1], MUL, bitlen, 1, no_status);
  share* rtwosq = c->PutFPGate(r[2], r[2], MUL, bitlen, 1, no_status);
  share* theta1 = c->PutFPGate(rzerosq, ronesq, ADD, bitlen, 1, no_status);
  share* theta2 = c->PutFPGate(theta1, rtwosq, ADD, bitlen, 1, no_status);
  share* theta = c->PutFPGate(theta2, SQRT);

  share* co = BuildCosCircuit(theta, c);
  share* si = BuildSinCircuit(theta, c);

  share* c1 = c->PutFPGate(one_gate, co, SUB, bitlen, 1, no_status);

  share* itheta = c->PutFPGate(one_gate, theta, DIV, bitlen, 1, no_status);
  // if theta == 0 then itheta = 0
  share* thetanotz = theta->get_wire_ids_as_share(0);
  for (uint32_t tt = 1; tt < bitlen - 1;
       tt++) {  // -1 -> do not include sign bit
    share* tempbit = theta->get_wire_ids_as_share(tt);
    share* temp = thetanotz;
    thetanotz = c->PutORGate(temp, tempbit);
    delete tempbit;
    delete temp;
  }
  share* tempitheta = itheta;
  itheta = c->PutMUXGate(itheta, zero_gate, thetanotz);
  delete tempitheta;

  share* x = c->PutFPGate(r[0], itheta, MUL, bitlen, 1, no_status);
  share* y = c->PutFPGate(r[1], itheta, MUL, bitlen, 1, no_status);
  share* z = c->PutFPGate(r[2], itheta, MUL, bitlen, 1, no_status);

  share* xx = c->PutFPGate(x, x, MUL, bitlen, 1, no_status);
  share* yy = c->PutFPGate(y, y, MUL, bitlen, 1, no_status);
  share* zz = c->PutFPGate(z, z, MUL, bitlen, 1, no_status);

  share* xy = c->PutFPGate(x, y, MUL, bitlen, 1, no_status);
  share* xz = c->PutFPGate(x, z, MUL, bitlen, 1, no_status);
  share* yz = c->PutFPGate(y, z, MUL, bitlen, 1, no_status);

  share* sx = c->PutFPGate(si, x, MUL, bitlen, 1, no_status);
  share* sy = c->PutFPGate(si, y, MUL, bitlen, 1, no_status);
  share* sz = c->PutFPGate(si, z, MUL, bitlen, 1, no_status);

  share* temp;

  temp = c->PutFPGate(c1, xx, MUL, bitlen, 1, no_status);
  R[0] = c->PutFPGate(temp, co, ADD, bitlen, 1, no_status);
  delete temp;

  temp = c->PutFPGate(c1, xy, MUL, bitlen, 1, no_status);
  R[1] = c->PutFPGate(temp, sz, SUB, bitlen, 1, no_status);
  delete temp;

  temp = c->PutFPGate(c1, xz, MUL, bitlen, 1, no_status);
  R[2] = c->PutFPGate(temp, sy, ADD, bitlen, 1, no_status);
  delete temp;

  temp = c->PutFPGate(c1, xy, MUL, bitlen, 1, no_status);
  R[4] = c->PutFPGate(temp, sz, ADD, bitlen, 1, no_status);
  delete temp;

  temp = c->PutFPGate(c1, yy, MUL, bitlen, 1, no_status);
  R[5] = c->PutFPGate(temp, co, ADD, bitlen, 1, no_status);
  delete temp;

  temp = c->PutFPGate(c1, yz, MUL, bitlen, 1, no_status);
  R[6] = c->PutFPGate(temp, sx, SUB, bitlen, 1, no_status);
  delete temp;

  temp = c->PutFPGate(c1, xz, MUL, bitlen, 1, no_status);
  R[8] = c->PutFPGate(temp, sy, SUB, bitlen, 1, no_status);
  delete temp;

  temp = c->PutFPGate(c1, yz, MUL, bitlen, 1, no_status);
  R[9] = c->PutFPGate(temp, sx, ADD, bitlen, 1, no_status);
  delete temp;

  temp = c->PutFPGate(c1, zz, MUL, bitlen, 1, no_status);
  R[10] = c->PutFPGate(temp, co, ADD, bitlen, 1, no_status);
  delete temp;

  delete rzerosq;
  delete ronesq;
  delete rtwosq;

  delete theta1;
  delete theta2;
  delete theta;
  delete itheta;

  delete c1;
  delete co;
  delete si;

  delete x;
  delete y;
  delete z;

  delete xx;
  delete yy;
  delete zz;

  delete xy;
  delete xz;
  delete yz;

  delete sx;
  delete sy;
  delete sz;

  delete one_gate;
  delete zero_gate;
}
