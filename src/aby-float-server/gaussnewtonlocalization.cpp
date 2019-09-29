#include <jlog.h>
#include <privacyconf.h>
#include <stdio.h>
#include <iostream>

#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>

#include <abycore/aby/abyparty.h>
#include <abycore/circuit/arithmeticcircuits.h>
#include <abycore/circuit/booleancircuits.h>
#include <abycore/circuit/circuit.h>
#include <abycore/sharing/sharing.h>

#include <invert.h>
#include <matmult.h>
#include <projectpoints.h>
#include <svd.h>
#include <twonormsq.h>
#include <util.h>

using namespace std;

const uint32_t bitlen = 32;

share* BuildGaussNewtonIteration(share* threeDPts[], share* y0[], int numPts,
                                 share* f, share* cx, share* cy, share* x[],
                                 BooleanCircuit* c, ABYParty* party,
                                 e_role role) {

  std::vector<Sharing*>& sharings = party->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  share* temp;

  float zero = 0.0;
  float one = 1.0;
  share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);

  // add constant ones to threeDPts
  for (int p = 0; p < numPts; p++) {
    threeDPts[3 * numPts + p] = one_gate;
  }

  // Camera params (3x3)
  share* K[9] = {f,
                 c->PutCONSGate((uint32_t*)&zero, bitlen),
                 cx,
                 c->PutCONSGate((uint32_t*)&zero, bitlen),
                 f,
                 cy,
                 c->PutCONSGate((uint32_t*)&zero, bitlen),
                 c->PutCONSGate((uint32_t*)&zero, bitlen),
                 c->PutCONSGate((uint32_t*)&one, bitlen)};

  float jacob_epsilon =
      JACOB_EPSILON;  // this is needed, can't &JACOB_EPSILON directly
  share* jacob_epsilon_gate = c->PutCONSGate((uint32_t*)&jacob_epsilon, bitlen);

  // project points using x guess
  share** yHomog = new share*[3 * numPts];
  BuildProjectPointsCircuit(threeDPts, x, K, yHomog, numPts, c, true);

  // throw away last "row" of result (constant ones),
  // interleave, and transpose y into 2nx1 vector
  // e.g. [x1; y1; x2; y2 ...]
  share** y = new share*[2 * numPts];
  for (int p = 0; p < numPts; ++p) {
    y[p * 2] = yHomog[p];
    y[(p * 2) + 1] = yHomog[p + numPts];
  }
  delete[] yHomog;

  // calculate jacobian (real 2D matrix not linearized)
  //  _                                            _
  // | du1/dr1 du1/dr2 du1/dr3 du1/t1 du1/t2 du1/t3 |
  // | dv1/dr1 dv1/dr2 dv1/dr3 dv1/t1 dv1/t2 dv1/t3 |
  // | du2/dr1 du2/dr2 du2/dr3 du2/t1 du2/t2 du2/t3 |
  // |                      .                       |
  // |                      .                       |
  // |_                     . (2n)                 _|
  share*** jacob = new share**[numPts * 2];
  for (int p = 0; p < numPts * 2; p++) {
    jacob[p] = new share*[6];
  }
  for (int j = 0; j < 6; j++) {  // for each DOF
    // perturb x
    temp = x[j];
    x[j] = c->PutFPGate(x[j], jacob_epsilon_gate, ADD, bitlen, 1, no_status);

    // project with epsilon
    share** ytempHomog = new share*[3 * numPts];
    BuildProjectPointsCircuit(threeDPts, x, K, ytempHomog, numPts, c, true);

    // throw away last "row" of result,
    // reshape, and transpose y into 2nx1 vector
    // e.g. [x1; y1; x2; y2 ...]
    share** ytemp = new share*[2 * numPts];
    for (int p = 0; p < numPts; ++p) {
      ytemp[p * 2] = ytempHomog[p];
      ytemp[(p * 2) + 1] = ytempHomog[p + numPts];
    }
    delete[] ytempHomog;

    // put into jacob
    for (int p = 0; p < 2 * numPts; ++p) {
      share* ytmy = c->PutFPGate(ytemp[p], y[p], SUB, bitlen, 1, no_status);
      jacob[p][j] =
          c->PutFPGate(ytmy, jacob_epsilon_gate, DIV, bitlen, 1, no_status);
      delete ytmy;
    }

    // un-perturb x
    delete x[j];
    x[j] = temp;
  }

  // calculate error
  // dy = y0 - y;
  share** dy = new share*[2 * numPts];
  for (int p = 0; p < 2 * numPts; p++) {
    dy[p] = c->PutFPGate(y0[p], y[p], SUB, bitlen, 1, no_status);
  }

  // compute pseudo inverse to solve:
  // dy = J dx   ->   dx = Jdagger dy
  // Note: invert requires real
  // 2D arrays not linearized versions
  share*** jacobI = new share**[6];  // 6x2n
  for (int p = 0; p < 6; p++) {
    jacobI[p] = new share*[numPts * 2];
  }

#if PPL_FLOW == PPL_FLOW_DO || PPL_FLOW == PPL_FLOW_SiSL
  BuildInvertCircuit(jacob, 2 * numPts, 6, jacobI, c, party, role, nullptr,
                     nullptr);

#elif PPL_FLOW == PPL_FLOW_LOOP_LEAK
  // Share Carryover.
  // Invert circuit calls SVD which runs multiple
  // circuit executions within. We need to store
  // dy and x variables across these executions as
  // raw shares.
  // First, define a place to store dy and x raw shares
  uint32_t* raw_dy = new uint32_t[2 * numPts];
  uint32_t* raw_x = new uint32_t[6];
  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    for (int i = 0; i < 2 * numPts; ++i) {
      temp = dy[i];
      dy[i] = bc->PutY2BGate(temp);
      delete temp;
    }
    for (int i = 0; i < 6; ++i) {
      temp = x[i];
      x[i] = bc->PutY2BGate(temp);
      delete temp;
    }
  }
  // Next, put output gates on the cirucuit
  for (int i = 0; i < 2 * numPts; ++i) {
    share* temp = dy[i];
    dy[i] = bc->PutSharedOUTGate(dy[i]);
    delete temp;
  }
  for (int i = 0; i < 6; ++i) {
    share* temp = x[i];
    x[i] = bc->PutSharedOUTGate(x[i]);
    delete temp;
  }
  // Next, build lambda function to convert share object to raw share
  std::function<void()> toRawShares = [dy, raw_dy, x, raw_x, numPts]() {
    for (int i = 0; i < 2 * numPts; ++i) {
      raw_dy[i] = dy[i]->get_clear_value<uint32_t>();
    }
    for (int i = 0; i < 6; ++i) {
      raw_x[i] = x[i]->get_clear_value<uint32_t>();
    }
  };
  // Lastly, build lambda function to convert back to share object
  std::function<void()> toShareObjects = [dy, raw_dy, x, raw_x, numPts, bc]() {
    for (int i = 0; i < 2 * numPts; ++i) {
      dy[i] = bc->PutSharedINGate(raw_dy[i], bitlen);
    }
    for (int i = 0; i < 6; ++i) {
      x[i] = bc->PutSharedINGate(raw_x[i], bitlen);
    }
  };

  BuildInvertCircuit(jacob, 2 * numPts, 6, jacobI, (BooleanCircuit*)c, party,
                     role, toRawShares, toShareObjects);

  // if yao was used, convert stored shares from bool back to yao circuit
  if (c->GetContext() == S_YAO) {
    for (int i = 0; i < 2 * numPts; ++i) {
      temp = dy[i];
      dy[i] = c->PutB2YGate(temp);
      delete temp;
    }
    for (int i = 0; i < 6; ++i) {
      temp = x[i];
      x[i] = c->PutB2YGate(temp);
      delete temp;
    }
  }

  // cleanup share carryover
  delete[] raw_dy;
  delete[] raw_x;
#endif

  // linearize for matmult
  share** jacobILinear = new share*[6 * 2 * numPts];
  for (int p = 0; p < 6; p++) {
    for (int pp = 0; pp < 2 * numPts; pp++) {
      jacobILinear[(p * 2 * numPts) + pp] = jacobI[p][pp];
    }
  }
  share* dx[6];
  BuildMatmultCircuit(jacobILinear, 6, 2 * numPts, dy, 2 * numPts,
                      1,  // column vector
                      dx, c);

  // cleanup jacobian
  for (int p = 0; p < numPts * 2; p++) {
    // individual jacob shares are deleted in invert->svd()
    delete[] jacob[p];
  }
  delete[] jacob;

  for (int p = 0; p < 6; p++) {
    delete[] jacobI[p];
  }
  delete[] jacobILinear;

  // return if error under threshold
  //  abs(norm(dy, cv::NORM_L2SQR)) < MIN_ER
  share* norm = BuildTwoNormSqCircuit(dx, 6, c);
  BuildFabsCircuit(norm, c);
  float min_er = MIN_ER;
  share* minErGate = c->PutCONSGate((uint32_t*)&min_er, bitlen);
  share* cmp = c->PutFPGate(minErGate, norm, CMP, bitlen, 1, no_status);
  delete norm;
  delete minErGate;

  // update pose
  for (int p = 0; p < 6; p++) {
#if PPL_FLOW == PPL_FLOW_DO  // set dx to zero if no error
    share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);
    share* olddx = dx[p];
    dx[p] = c->PutMUXGate(dx[p], zero_gate, cmp);
    delete olddx;
    delete zero_gate;
#endif
    temp = x[p];
    x[p] = c->PutFPGate(x[p], dx[p], ADD, bitlen, 1, no_status);
    delete temp;
    delete dx[p];
  }

  return cmp;
}

// runs circuit using pre-shared inputs
bool RunGaussNewtonIteration(uint32_t* threeDPts, uint32_t* y0, int numPts,
                             uint32_t f, uint32_t cx, uint32_t cy, uint32_t* x,
                             BooleanCircuit* c, ABYParty* p, e_role role) {

  std::vector<Sharing*>& sharings = p->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  // Allocate space for shares
  share** s_threeDPts = new share*[4 * numPts];  // leave room for contant 1's
  share** s_y0 = new share*[2 * numPts];
  share** s_x = new share*[6];
  share* s_f;
  share* s_cx;
  share* s_cy;
  share* temp;

  // Prepare inputs
  for (int p = 0; p < 3 * numPts; p++) {
    s_threeDPts[p] = bc->PutSharedINGate(threeDPts[p], bitlen);
  }
  for (int p = 0; p < 2 * numPts; p++) {
    s_y0[p] = bc->PutSharedINGate(y0[p], bitlen);
  }
  for (int p = 0; p < 6; p++) {
    s_x[p] = bc->PutSharedINGate(x[p], bitlen);
  }
  s_f = bc->PutSharedINGate(f, bitlen);
  s_cx = bc->PutSharedINGate(cx, bitlen);
  s_cy = bc->PutSharedINGate(cy, bitlen);
  // if yao, convert bool inputs to yao
  if (c->GetContext() == S_YAO) {
    for (int p = 0; p < 3 * numPts; p++) {
      temp = s_threeDPts[p];
      s_threeDPts[p] = c->PutB2YGate(temp);
      delete temp;
    }
    for (int p = 0; p < 2 * numPts; p++) {
      temp = s_y0[p];
      s_y0[p] = c->PutB2YGate(temp);
      delete temp;
    }
    for (int p = 0; p < 6; p++) {
      temp = s_x[p];
      s_x[p] = c->PutB2YGate(temp);
      delete temp;
    }
    temp = s_f;
    s_f = c->PutB2YGate(temp);
    delete temp;
    temp = s_cx;
    s_cx = c->PutB2YGate(temp);
    delete temp;
    temp = s_cy;
    s_cy = c->PutB2YGate(temp);
    delete temp;
  }
  // prepare circuit
  share* s_done = BuildGaussNewtonIteration(s_threeDPts, s_y0, numPts, s_f,
                                            s_cx, s_cy, s_x, c, p, role);

  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    for (int p = 0; p < 6; p++) {
      temp = s_x[p];
      s_x[p] = bc->PutY2BGate(temp);
      delete temp;
    }
  }

  // prepare output (only need x and done flag)
  for (int p = 0; p < 6; p++) {
    temp = s_x[p];
    s_x[p] = bc->PutSharedOUTGate(temp);  // raw share
    delete temp;
  }
  temp = s_done;
  s_done = c->PutOUTGate(temp, ALL);  // need cleartext
  delete temp;

  p->ExecCircuit();

  // get shared output
  for (int p = 0; p < 6; p++) {
    x[p] = s_x[p]->get_clear_value<uint32_t>();
    delete s_x[p];
  }
  delete[] s_x;
  uint32_t done = s_done->get_clear_value<uint32_t>();
  for (int p = 0; p < 3 * numPts; p++) {  // ignore constant ones?
    delete s_threeDPts[p];
  }
  delete[] s_threeDPts;
  for (int p = 0; p < 2 * numPts; p++) {
    delete s_y0[p];
  }
  delete[] s_y0;
  delete s_f;
  delete s_cx;
  delete s_cy;

  collectTiming();
  collectCommunication();
  p->Reset();
  return done != 0;
}

// uint32_t arguments must be from PutSharedOUTGate()
// then calling get_clear_value() and circuit->Reset();
// This function builds/executes multiple circuits.
void RunGaussNewtonCircuit(uint32_t* threeDPts, uint32_t* y0, int numPts,
                           uint32_t f, uint32_t cx, uint32_t cy, uint32_t* x,
                           BooleanCircuit* c, ABYParty* party, e_role role) {
  // GN Iteration
  for (int i = 0; i < GN_MAX_ITR; i++) {
    if (_verbosity & DBG_FLOW) {
      cout << RED << "\n\nGN Iteration " << i << RESET << endl;
    }

    // break if error under threshold
    if (RunGaussNewtonIteration(threeDPts, y0, numPts, f, cx, cy, x, c, party,
                                role))
      break;
  }
}

// Wrapper function around RunGaussNewtonCircuit which
// takes shares instead of secret shares.
// This (and all inner function calls) creates multiple
// circuits such that whenever an intermediate plaintext value
// is required, the circuit ends and a new circuit begins.
//
// Doing so has the overhead of converting intermediate values
// to secret shares to be used in the next circuit execution.
// It does not however, require any of the "debug" functionality
// from the ABY library to retreive intermediate values.
//
// This wrapper creates dummy circuits to build secret shares
// to be passed to the top level GN function.
//
// Note: All 2D matrix inputs passed as single dimension array
// Note: the constant 1's can be ignored (but space must be allocated)
//
// 3D World Points (4xnumPts) which look like:
//  _             _
// |  x1 x2 x3     |
// |  y1 y2 y3 ... |
// |  z1 z2 z3     |
// |_ 1  1  1     _|
// must be pass as single dimension:
// [x1, x2, ...;  y1, y2, ...;  z1, z2, ...; 1, 1, ... ]
//
// 2D Image Points (3xnumPts) which look like:
//  _             _
// |  u1 u2 u3     |
// |  v1 v2 v3 ... |
// |_ 1  1  1     _|
// must be passed as single dimension:
// [ u1, u2, ...;  v1, v2, ...;  1, 1, ... ]
//
// x (initial guess) : 6x1 [ rotation angles, translations ]
void BuildAndRunGaussNewtonLoopLeak(share* s_threeDPts[], share* s_twoDPts[],
                                    int numPts, share* s_f, share* s_cx,
                                    share* s_cy, share* s_x[],
                                    BooleanCircuit* c, ABYParty* party,
                                    e_role role) {

  std::vector<Sharing*>& sharings = party->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  share* temp;

  // First transform input shares
  //
  // throw away last "row" of 2D points,
  // reshape, and transpose into 2nx1 vector
  // e.g. [x1; y1; x2; y2 ...]
  share** s_y0 = new share*[2 * numPts];
  for (int p = 0; p < numPts; ++p) {
    s_y0[2 * p] = s_twoDPts[p];
    s_y0[(2 * p) + 1] = s_twoDPts[p + numPts];
  }

  // Next, run dummy circuit to create raw shares from share objects
  //   stores secret-shared output
  uint32_t* threeDPts = new uint32_t[4 * numPts];
  uint32_t* y0 = new uint32_t[2 * numPts];
  uint32_t* x = new uint32_t[6];
  uint32_t f;
  uint32_t cx;
  uint32_t cy;

  // if yao is prefered circuit and yao shares passed in, must
  // convert to bool to get raw shares.
  // Tests must pass in bool shares even if yao is preferred.
  if (c->GetContext() == S_YAO && s_threeDPts[0]->get_share_type() == S_YAO) {
    DEBUG_MSG("Converting yao shares to bool for shared output\n");
    for (int p = 0; p < 3 * numPts; p++) {
      temp = s_threeDPts[p];
      s_threeDPts[p] = bc->PutY2BGate(temp);
      delete temp;
    }
    for (int p = 0; p < 2 * numPts; p++) {
      temp = s_y0[p];
      s_y0[p] = bc->PutY2BGate(temp);
      delete temp;
    }
    for (int p = 0; p < 6; p++) {
      temp = s_x[p];
      s_x[p] = bc->PutY2BGate(temp);
      delete temp;
    }
    temp = s_f;
    s_f = bc->PutY2BGate(temp);
    delete temp;
    temp = s_cx;
    s_cx = bc->PutY2BGate(temp);
    delete temp;
    temp = s_cy;
    s_cy = bc->PutY2BGate(temp);
    delete temp;
  }

  //   build shared output objects to get share, ignore constant 1s
  for (int p = 0; p < 3 * numPts; p++) {
    temp = s_threeDPts[p];
    s_threeDPts[p] = bc->PutSharedOUTGate(temp);
    delete temp;
  }
  for (int p = 0; p < 2 * numPts; p++) {
    temp = s_y0[p];
    s_y0[p] = bc->PutSharedOUTGate(temp);
    delete temp;
  }
  for (int p = 0; p < 6; p++) {
    temp = s_x[p];
    s_x[p] = bc->PutSharedOUTGate(temp);
    delete temp;
  }
  temp = s_f;
  s_f = bc->PutSharedOUTGate(temp);
  delete temp;
  temp = s_cx;
  s_cx = bc->PutSharedOUTGate(temp);
  delete temp;
  temp = s_cy;
  s_cy = bc->PutSharedOUTGate(temp);
  delete temp;
  //   run the dummy circuit
  CLOCK(ShareInputs);
  TIC(ShareInputs);
  party->ExecCircuit();
  TOC(ShareInputs);
  //   get shared outputs, ignore constant 1s
  for (int p = 0; p < 3 * numPts; p++) {
    threeDPts[p] = s_threeDPts[p]->get_clear_value<uint32_t>();
    delete s_threeDPts[p];
  }
  for (int p = 0; p < 2 * numPts; p++) {
    y0[p] = s_y0[p]->get_clear_value<uint32_t>();
    delete s_y0[p];
  }
  for (int p = 0; p < 6; p++) {
    x[p] = s_x[p]->get_clear_value<uint32_t>();
    delete s_x[p];
  }
  f = s_f->get_clear_value<uint32_t>();
  delete s_f;
  cx = s_cx->get_clear_value<uint32_t>();
  delete s_cx;
  cy = s_cy->get_clear_value<uint32_t>();
  delete s_cy;

  collectTiming();
  collectCommunication();
  party->Reset();
  // arrays now contain secret shared values of plaintext from each party

  CLOCK(GaussNewtonGD);
  TIC(GaussNewtonGD);
  RunGaussNewtonCircuit(threeDPts, y0, numPts, f, cx, cy, x, c, party, role);
  TOC(GaussNewtonGD);

  // next start another dummy circuit to return share
  // objects which are built from shares
  for (int p = 0; p < 6; p++) {
    s_x[p] = bc->PutSharedINGate(&x[p], bitlen);
  }
}

// Data Oblivious version of GN. Does not require executing
// multiple circuits to retreive intermediate values.
//
// Note: All 2D matrix inputs passed as single dimension array
// Note: the constant 1's can be ignored (but space must be allocated)
//
// 3D World Points (4xnumPts) which look like:
//  _             _
// |  x1 x2 x3     |
// |  y1 y2 y3 ... |
// |  z1 z2 z3     |
// |_ 1  1  1     _|
// must be pass as single dimension:
// [x1, x2, ...;  y1, y2, ...;  z1, z2, ...; 1, 1, ... ]
//
// 2D Image Points (3xnumPts) which look like:
//  _             _
// |  u1 u2 u3     |
// |  v1 v2 v3 ... |
// |_ 1  1  1     _|
// must be passed as single dimension:
// [ u1, u2, ...;  v1, v2, ...;  1, 1, ... ]
//
// x (initial guess) : 6x1 [ rotation angles, translations ]
void BuildAndRunGaussNewtonDO(share* s_threeDPts[], share* s_twoDPts[],
                              int numPts, share* s_f, share* s_cx, share* s_cy,
                              share* s_x[], BooleanCircuit* c, ABYParty* party,
                              e_role role) {
  // First transform input shares
  //
  // throw away last "row" of 2D points,
  // reshape, and transpose into 2nx1 vector
  // e.g. [x1; y1; x2; y2 ...]
  share** s_y0 = new share*[2 * numPts];
  for (int p = 0; p < numPts; ++p) {
    s_y0[2 * p] = s_twoDPts[p];
    s_y0[(2 * p) + 1] = s_twoDPts[p + numPts];
  }

  // GN Iteration
  for (int i = 0; i < GN_MAX_ITR; i++) {
    if (_verbosity & DBG_FLOW) {
      cout << RED << "\n\nGN Iteration " << i << RESET << endl;
    }

    share* s_done = BuildGaussNewtonIteration(s_threeDPts, s_y0, numPts, s_f,
                                              s_cx, s_cy, s_x, c, party, role);
    delete s_done;  // throw it away
  }
}

void BuildAndRunGaussNewton(share* s_threeDPts[], share* s_twoDPts[],
                            int numPts, share* s_f, share* s_cx, share* s_cy,
                            share* s_x[], BooleanCircuit* c, ABYParty* party,
                            e_role role) {
#if PPL_FLOW == PPL_FLOW_DO
  BuildAndRunGaussNewtonDO(s_threeDPts, s_twoDPts, numPts, s_f, s_cx, s_cy, s_x,
                           c, party, role);
#elif PPL_FLOW == PPL_FLOW_LOOP_LEAK || PPL_FLOW == PPL_FLOW_SiSL
  BuildAndRunGaussNewtonLoopLeak(s_threeDPts, s_twoDPts, numPts, s_f, s_cx,
                                 s_cy, s_x, c, party, role);
#endif
}
