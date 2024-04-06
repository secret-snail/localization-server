#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

#include <rodrigues.h>

using namespace std;

const uint32_t bitlen = 32;

share* BuildSinCircuit(share* theta, BooleanCircuit* c) {
  float pi = M_PI;
  share* pi_gate = c->PutCONSGate((uint32_t*)&pi, bitlen);
  share* adjusted = c->PutFPGate(theta, pi_gate, DIV, bitlen, 1, no_status);
  share* si = c->PutFPGate(
      adjusted, SIN);  // sets bitlength to 40 (?), other circuit is 33 (?)
  // cout << "sine bitlength " <<si->get_bitlength()<<endl;
  si->set_bitlength(32);
  delete pi_gate;
  delete adjusted;
  return si;

  // Below uses taylor series approx
  ////return x - (pow(x,3)/(3*2)) + (pow(x,5)/(5*4*3*2));// -
  ///(pow(x,7)/(7*6*5*4*3*2));
  // share *t1 = c->PutFPGate(theta, theta, MUL, bitlen, 1, no_status);
  // t1 = c->PutFPGate(t1, theta, MUL, bitlen, 1, no_status);
  // float threeFact = 3*2;
  // share* tFGate = c->PutCONSGate((uint32_t*)&threeFact, bitlen);
  // t1 = c->PutFPGate(t1, tFGate, DIV, bitlen, 1, no_status);

  // share *t2 = c->PutFPGate(theta, theta, MUL, bitlen, 1, no_status);
  // t2 = c->PutFPGate(t2, theta, MUL, bitlen, 1, no_status);
  // t2 = c->PutFPGate(t2, theta, MUL, bitlen, 1, no_status);
  // t2 = c->PutFPGate(t2, theta, MUL, bitlen, 1, no_status);
  // float fiveFact = 5*4*3*2;
  // share* fFGate = c->PutCONSGate((uint32_t*)&fiveFact, bitlen);
  // t2 = c->PutFPGate(t2, fFGate, DIV, bitlen, 1, no_status);

  // share *res = c->PutFPGate(theta, t1, SUB, bitlen, 1, no_status);
  // res = c->PutFPGate(res, t2, ADD, bitlen, 1, no_status);

  // delete t1;
  // delete tFGate;
  // delete t2;
  // delete fFGate;

  // return res;
}

share* BuildCosCircuit(share* theta, BooleanCircuit* c) {
  float pi = M_PI;
  share* pi_gate = c->PutCONSGate((uint32_t*)&pi, bitlen);
  share* adjusted = c->PutFPGate(theta, pi_gate, DIV, bitlen, 1, no_status);
  share* co = c->PutFPGate(
      adjusted, COS);  // sets bitlength to 40 (?), other circuit is 33 (?)
  // cout << "sine bitlength " <<si->get_bitlength()<<endl;
  co->set_bitlength(32);
  delete pi_gate;
  delete adjusted;
  return co;

  // Below uses taylor series approx

  ////return 1 - (pow(x,2)/2) + (pow(x,4)/(4*3*2));// - (pow(x,6)/(6*5*4*3*2));
  // share *t1 = c->PutFPGate(theta, theta, MUL, bitlen, 1, no_status);
  // float twoFact = 2;
  // share* tFGate = c->PutCONSGate((uint32_t*)&twoFact, bitlen);
  // t1 = c->PutFPGate(t1, tFGate, DIV, bitlen, 1, no_status);

  // share *t2 = c->PutFPGate(theta, theta, MUL, bitlen, 1, no_status);
  // t2 = c->PutFPGate(t2, theta, MUL, bitlen, 1, no_status);
  // t2 = c->PutFPGate(t2, theta, MUL, bitlen, 1, no_status);
  // float fourFact = 4*3*2;
  // share* fFGate = c->PutCONSGate((uint32_t*)&fourFact, bitlen);
  // t2 = c->PutFPGate(t2, fFGate, DIV, bitlen, 1, no_status);

  // float one = 1;
  // share* oneGate = c->PutCONSGate((uint32_t*)&one, bitlen);
  // share *res = c->PutFPGate(oneGate, t1, SUB, bitlen, 1, no_status);
  // res = c->PutFPGate(res, t2, ADD, bitlen, 1, no_status);

  // delete t1;
  // delete tFGate;
  // delete t2;
  // delete fFGate;

  // return res;
}
