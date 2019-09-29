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

// returns error norm (not if lower than min), and this lambda
std::pair<share*, share*> BuildLMIteration(share* threeDPts[], share* y0[],
                                           int numPts, share* f, share* cx,
                                           share* cy, share* x[],
                                           share* prevErrNorm, share* lambda,
                                           BooleanCircuit* c, ABYParty* party,
                                           e_role role) {

  if (_verbosity & DBG_ARGS) {
    for (int p = 0; p < 3 * numPts; p++) {
      char s[50];
      sprintf(s, "threeDPts[%d]", p);
      c->PutPrintValueGate(threeDPts[p], s);
    }
    for (int p = 0; p < 2 * numPts; p++) {
      char s[50];
      sprintf(s, "y0[%d]", p);
      c->PutPrintValueGate(y0[p], s);
    }
    c->PutPrintValueGate(f, "s_f");
    c->PutPrintValueGate(cx, "s_cx");
    c->PutPrintValueGate(cx, "s_cy");
  }

  std::vector<Sharing*>& sharings = party->GetSharings();

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

  float jacob_epsilon = JACOB_EPSILON;
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
  if (_verbosity & DBG_PROJECT) {
    for (int p = 0; p < 2 * numPts; p++) {
      char s[50];
      sprintf(s, "y[%d]", p);
      c->PutPrintValueGate(y[p], s);
    }
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
  if (_verbosity & DBG_JACOB) {
    for (int p = 0; p < numPts * 2; p++) {
      for (int pp = 0; pp < 6; pp++) {
        char s[50];
        sprintf(s, "jacob[%d][%d]", p, pp);
        c->PutPrintValueGate(jacob[p][pp], s);
      }
    }
  }

  // calculate error
  // dy = y0 - y;
  share** dy = new share*[2 * numPts];
  for (int p = 0; p < 2 * numPts; p++) {
    dy[p] = c->PutFPGate(y0[p], y[p], SUB, bitlen, 1, no_status);
  }

  // LM with fletcher improvement
  // inv(Jt*J + lambda*diag(Jt*J)) * Jt * dy
  share*** JtJ = new share**[6];  // 6x6
  for (int p = 0; p < 6; p++) {
    JtJ[p] = new share*[6];
  }
  BuildMatmult2DwTransposeCircuit(jacob, numPts * 2, 6, true, jacob, numPts * 2,
                                  6, false, JtJ,
                                  c);  // must use 2D arrays, not 1D
  // add lambda*diag(Jt*J) in-place
  for (int p = 0; p < 6; p++) {
    // add fletcher column-wise to JtJ
    share* fletcher =
        c->PutFPGate(lambda, JtJ[p][p], MUL, bitlen, 1, no_status);

    for (int pp = 0; pp < 6; pp++) {
      share* temp = JtJ[pp][p];
      JtJ[pp][p] =
          c->PutFPGate(JtJ[pp][p], fletcher, ADD, bitlen, 1, no_status);
      delete temp;
    }
    delete fletcher;
  }

  // compute pseudo inverse
  // Note: invert requires real
  // 2D arrays not linearized versions
  share*** JtJ_i = new share**[6];  // 6x6
  for (int p = 0; p < 6; p++) {
    JtJ_i[p] = new share*[6];
  }

#if PPL_FLOW == PPL_FLOW_DO || PPL_FLOW == PPL_FLOW_SiSL
  BuildInvertCircuit(JtJ, 6, 6, JtJ_i, c, party, role, nullptr, nullptr);

#elif PPL_FLOW == PPL_FLOW_LOOP_LEAK
  // Share Carryover.
  // Invert circuit calls SVD which runs multiple
  // circuit executions within. We need to store
  // dy and jacob variables across these executions as
  // raw shares.
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  // First, define a place to store dy and jacob raw shares
  uint32_t* raw_dy = new uint32_t[2 * numPts];
  uint32_t* raw_jacob = new uint32_t[6 * 2 * numPts];
  uint32_t* raw_x = new uint32_t[6];
  uint32_t* raw_lambda = new uint32_t;
  uint32_t* raw_prevErrNorm = new uint32_t;
  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    DEBUG_MSG("converting to yao\n");
    for (int i = 0; i < 2 * numPts; ++i) {
      temp = dy[i];
      dy[i] = bc->PutY2BGate(temp);
      delete temp;
    }
    for (int p = 0; p < 2 * numPts; ++p) {
      for (int pp = 0; pp < 6; ++pp) {
        temp = jacob[p][pp];
        jacob[p][pp] = bc->PutY2BGate(temp);
        delete temp;
      }
    }
    for (int i = 0; i < 6; ++i) {
      temp = x[i];
      x[i] = bc->PutY2BGate(temp);
      delete temp;
    }
    temp = lambda;
    lambda = bc->PutY2BGate(temp);
    delete temp;
    temp = prevErrNorm;
    prevErrNorm = bc->PutY2BGate(temp);
    delete temp;
  }
  // Next, put output gates on the cirucuit
  for (int i = 0; i < 2 * numPts; ++i) {
    temp = dy[i];
    dy[i] = bc->PutSharedOUTGate(dy[i]);
    delete temp;
  }
  for (int p = 0; p < 2 * numPts; ++p) {
    for (int pp = 0; pp < 6; ++pp) {
      temp = jacob[p][pp];
      jacob[p][pp] = bc->PutSharedOUTGate(jacob[p][pp]);
      delete temp;
    }
  }
  for (int i = 0; i < 6; ++i) {
    temp = x[i];
    x[i] = bc->PutSharedOUTGate(x[i]);
    delete temp;
  }
  temp = lambda;
  lambda = bc->PutSharedOUTGate(temp);
  delete temp;
  temp = prevErrNorm;
  prevErrNorm = bc->PutSharedOUTGate(temp);
  delete temp;
  // Next, build lambda function to convert share object to raw share
  std::function<void()> toRawShares = [dy, raw_dy, jacob, raw_jacob, x, raw_x,
                                       lambda, raw_lambda, prevErrNorm,
                                       raw_prevErrNorm, numPts]() {
    for (int i = 0; i < 2 * numPts; ++i) {
      raw_dy[i] = dy[i]->get_clear_value<uint32_t>();
    }
    for (int p = 0; p < 2 * numPts; ++p) {
      for (int pp = 0; pp < 6; ++pp) {
        raw_jacob[p * 6 + pp] = jacob[p][pp]->get_clear_value<uint32_t>();
      }
    }
    for (int i = 0; i < 6; ++i) {
      raw_x[i] = x[i]->get_clear_value<uint32_t>();
    }
    *raw_lambda = lambda->get_clear_value<uint32_t>();
    *raw_prevErrNorm = prevErrNorm->get_clear_value<uint32_t>();
  };
  // Lastly, build lambda function to convert back to share object
  share** ref_lambda = &lambda;
  share** ref_prevErrNorm = &prevErrNorm;
  std::function<void()> toShareObjects =
      [dy, raw_dy, jacob, raw_jacob, x, raw_x, ref_lambda, raw_lambda,
       ref_prevErrNorm, raw_prevErrNorm, numPts, bc]() {
        for (int i = 0; i < 2 * numPts; ++i) {
          dy[i] = bc->PutSharedINGate(raw_dy[i], bitlen);
        }
        for (int p = 0; p < 2 * numPts; ++p) {
          for (int pp = 0; pp < 6; ++pp) {
            jacob[p][pp] = bc->PutSharedINGate(raw_jacob[p * 6 + pp], bitlen);
          }
        }
        for (int i = 0; i < 6; ++i) {
          x[i] = bc->PutSharedINGate(raw_x[i], bitlen);
        }
        *ref_lambda = bc->PutSharedINGate(*raw_lambda, bitlen);
        *ref_prevErrNorm = bc->PutSharedINGate(*raw_prevErrNorm, bitlen);
      };

  BuildInvertCircuit(JtJ, 6, 6, JtJ_i, (BooleanCircuit*)c, party, role,
                     toRawShares, toShareObjects);

  // if yao was used, convert stored shares from bool back to yao circuit
  if (c->GetContext() == S_YAO) {
    for (int i = 0; i < 2 * numPts; ++i) {
      temp = dy[i];
      dy[i] = c->PutB2YGate(temp);
      delete temp;
    }
    for (int p = 0; p < 2 * numPts; ++p) {
      for (int pp = 0; pp < 6; ++pp) {
        temp = jacob[p][pp];
        jacob[p][pp] = c->PutB2YGate(temp);
        delete temp;
      }
    }
    for (int i = 0; i < 6; ++i) {
      temp = x[i];
      x[i] = c->PutB2YGate(temp);
      delete temp;
    }
    temp = lambda;
    lambda = c->PutB2YGate(temp);
    delete temp;
    temp = prevErrNorm;
    prevErrNorm = c->PutB2YGate(temp);
    delete temp;
  }

  // cleanup share carryover
  delete[] raw_dy;
  delete[] raw_jacob;
  delete[] raw_x;
  delete raw_lambda;
  delete raw_prevErrNorm;
#endif

  share*** JtJ_i_Jt = new share**[6];  // 6x2n
  for (int p = 0; p < 6; p++) {
    JtJ_i_Jt[p] = new share*[numPts * 2];
  }
  BuildMatmult2DwTransposeCircuit(JtJ_i, 6, 6, false, jacob, numPts * 2, 6,
                                  true, JtJ_i_Jt, c);

  // linearize for matmult
  share** JtJ_i_Jt_linear = new share*[6 * 2 * numPts];
  for (int p = 0; p < 6; p++) {
    for (int pp = 0; pp < 2 * numPts; pp++) {
      JtJ_i_Jt_linear[(p * 2 * numPts) + pp] = JtJ_i_Jt[p][pp];
    }
  }
  share* dx[6];
  BuildMatmultCircuit(JtJ_i_Jt_linear, 6, 2 * numPts, dy, 2 * numPts,
                      1,  // column vector
                      dx, c);

  // cleanup
  delete[] JtJ_i_Jt_linear;
  for (int p = 0; p < numPts * 2; p++) {
    delete[] jacob[p];
  }
  delete[] jacob;

  for (int p = 0; p < 6; p++) {
    delete[] JtJ[p];
    delete[] JtJ_i[p];
    delete[] JtJ_i_Jt[p];
  }
  delete[] JtJ;
  delete[] JtJ_i;
  delete[] JtJ_i_Jt;

  // Constants used after invert (invert may break circuit, define after)
  float min_er = MIN_ER;
  share* min_er_gate = c->PutCONSGate((uint32_t*)&min_er, bitlen);

  float lambda_max = LM_LAMBDA_MAX;
  share* lambda_max_gate = c->PutCONSGate((uint32_t*)&lambda_max, bitlen);

  float lambda_min = LM_LAMBDA_MIN;
  share* lambda_min_gate = c->PutCONSGate((uint32_t*)&lambda_min, bitlen);

  float ten = 10;
  share* ten_gate = c->PutCONSGate((uint32_t*)&ten, bitlen);

  // check error under threshold
  //  abs(norm(dy, cv::NORM_L2SQR)) < MIN_ER
  share* absnorm = BuildTwoNormSqCircuit(dx, 6, c);
  BuildFabsCircuit(absnorm, c);
  share* er_flag =
      c->PutFPGate(min_er_gate, absnorm, CMP, bitlen, 1, no_status);
  if (_verbosity & DBG_ER) {
    c->PutPrintValueGate(absnorm, "absnorm");
  }

  share* biggerLambda =
      c->PutFPGate(lambda, ten_gate, MUL, bitlen, 1, no_status);
  share* smallerLambda =
      c->PutFPGate(lambda, ten_gate, DIV, bitlen, 1, no_status);

  share* errWentDown =
      c->PutFPGate(prevErrNorm, absnorm, CMP, bitlen, 1, no_status);
  share* oldlam = lambda;
  lambda = c->PutMUXGate(smallerLambda, biggerLambda, errWentDown);
  delete oldlam;
  delete smallerLambda;
  delete biggerLambda;
  delete errWentDown;

  share* lambdaTooBig =
      c->PutFPGate(lambda, lambda_max_gate, CMP, bitlen, 1, no_status);
  oldlam = lambda;
  lambda = c->PutMUXGate(lambda_max_gate, lambda, lambdaTooBig);
  delete oldlam;
  delete lambdaTooBig;
  share* lambdaTooSmall =
      c->PutFPGate(lambda_min_gate, lambda, CMP, bitlen, 1, no_status);
  oldlam = lambda;
  share* newlambda = c->PutMUXGate(lambda_min_gate, lambda, lambdaTooSmall);
  delete oldlam;
  delete lambdaTooSmall;
  if (_verbosity & DBG_LAMBDA) {
    c->PutPrintValueGate(newlambda, "lambda");
  }

  // update pose - note: this is usually done after checking error
  for (int p = 0; p < 6; p++) {
#if PPL_FLOW == PPL_FLOW_DO  // set dx to zero if no error
    share* olddx = dx[p];
    share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);
    dx[p] = c->PutMUXGate(dx[p], zero_gate, er_flag);
    delete olddx;
#endif
    temp = x[p];
    x[p] = c->PutFPGate(x[p], dx[p], ADD, bitlen, 1, no_status);
    delete temp;
    delete dx[p];
  }

  delete er_flag;
  return {absnorm, newlambda};
}

// runs circuit using pre-shared inputs
bool RunLMIteration(uint32_t* threeDPts, uint32_t* y0, int numPts, uint32_t f,
                    uint32_t cx, uint32_t cy, uint32_t* x,
                    uint32_t* prevErrNorm, uint32_t* lambda, BooleanCircuit* c,
                    ABYParty* p, e_role role) {

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
  share* s_prevErrNorm;
  share* s_lambda;
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
  s_prevErrNorm = bc->PutSharedINGate(*prevErrNorm, bitlen);
  s_lambda = bc->PutSharedINGate(*lambda, bitlen);
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
    temp = s_prevErrNorm;
    s_prevErrNorm = c->PutB2YGate(temp);
    delete temp;
    temp = s_lambda;
    s_lambda = c->PutB2YGate(temp);
    delete temp;
  }
  // build circuit
  auto s_res = BuildLMIteration(s_threeDPts, s_y0, numPts, s_f, s_cx, s_cy, s_x,
                                s_prevErrNorm, s_lambda, c, p, role);

  // manage the outputs
  share* s_absnorm = s_res.first;
  share* s_newlambda = s_res.second;
  float min_er = MIN_ER;
  share* s_minErGate = c->PutCONSGate((uint32_t*)&min_er, bitlen);
  share* s_done =
      c->PutFPGate(s_minErGate, s_absnorm, CMP, bitlen, 1, no_status);
  delete s_minErGate;

  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    for (int p = 0; p < 6; p++) {
      temp = s_x[p];
      s_x[p] = bc->PutY2BGate(temp);
      delete temp;
    }
    temp = s_absnorm;
    s_absnorm = bc->PutY2BGate(temp);
    delete temp;
    temp = s_newlambda;
    s_newlambda = bc->PutY2BGate(temp);
    delete temp;
    // do not convert done, it is not raw share
  }

  // prepare output (only need x and done flag)
  for (int p = 0; p < 6; p++) {
    temp = s_x[p];
    s_x[p] = bc->PutSharedOUTGate(temp);  // raw share
    delete temp;
  }
  temp = s_absnorm;
  s_absnorm = bc->PutSharedOUTGate(temp);  // raw share
  delete temp;
  temp = s_newlambda;
  s_newlambda = bc->PutSharedOUTGate(temp);  // raw share
  delete temp;
  temp = s_done;
  s_done = c->PutOUTGate(temp, ALL);  // loop leak needs cleartext
  delete temp;

  p->ExecCircuit();

  // get shared output
  for (int p = 0; p < 6; p++) {
    x[p] = s_x[p]->get_clear_value<uint32_t>();
    delete s_x[p];
  }
  delete[] s_x;
  *prevErrNorm = s_absnorm->get_clear_value<uint32_t>();
  delete s_absnorm;
  *lambda = s_newlambda->get_clear_value<uint32_t>();
  delete s_newlambda;
  uint32_t done = s_done->get_clear_value<uint32_t>();
  delete s_done;

  // cleanup
  for (int p = 0; p < 3 * numPts; p++) {  // ignore constant ones?
    delete s_threeDPts[p];
  }
  delete[] s_threeDPts;
  for (int p = 0; p < 2 * numPts; p++) {
    delete s_y0[p];
  }
  delete[] s_y0;

  collectTiming();
  collectCommunication();
  p->Reset();
  return done != 0;
}

// uint32_t arguments must be from PutSharedOUTGate()
// then calling get_clear_value() and circuit->Reset();
// This function builds/executes multiple circuits.
void RunLMCircuit(uint32_t* threeDPts, uint32_t* y0, int numPts, uint32_t f,
                  uint32_t cx, uint32_t cy, uint32_t* x, uint32_t* prevErrNorm,
                  uint32_t* lambda, BooleanCircuit* c, ABYParty* party,
                  e_role role) {
  // LM Iteration
  for (int i = 0; i < LM_MAX_ITR; i++) {
    if (_verbosity & DBG_PROJECT) {
      cout << RED << "\n\nLM Iteration " << i << RESET << endl;
    }

    bool done = RunLMIteration(threeDPts, y0, numPts, f, cx, cy, x, prevErrNorm,
                               lambda, c, party, role);

    // break if error under threshold
    if (done)
      break;
  }
}

// Wrapper function around RunLMCircuit which
// takes shares instead of secret shares.
// This (and inner function calls) creates multiple
// circuits such that whenever an intermediate plaintext value
// is required, the circuit ends and a new circuit begins.
//
// Doing so has the overhead of converting intermediate values
// to secret shares to be used in the next circuit execution.
// It does not, however, require any of the "debug" functionality
// from the ABY library to retreive intermediate values.
//
// This wrapper creates dummy circuits to build secret shares
// to be passed to the top level LM function.
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
void BuildAndRunLMLoopLeak(share* s_threeDPts[], share* s_twoDPts[], int numPts,
                           share* s_f, share* s_cx, share* s_cy, share* s_x[],
                           BooleanCircuit* c, ABYParty* party, e_role role) {
  std::vector<Sharing*>& sharings = party->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  share* temp;

  // initialize state required between iterations
  float prevErrNorm_init = std::numeric_limits<float>::max();
  share* s_prevErrNorm = bc->PutCONSGate((uint32_t*)&prevErrNorm_init, bitlen);
  float lambda_init = LM_LAMBDA_INIT;
  share* s_lambda = bc->PutCONSGate((uint32_t*)&lambda_init, bitlen);

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
  uint32_t prevErrNorm;
  uint32_t lambda;

  // if yao is prefered circuit and yao shares passed in, must
  // convert to bool to get raw shares.
  // Tests must pass in bool shares even if yao is preferred.
  if (c->GetContext() == S_YAO && s_threeDPts[0]->get_share_type() == S_YAO) {
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
    // dont need prevErrNorm or lambda, they're always created with bool circuit
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
  temp = s_prevErrNorm;
  s_prevErrNorm = bc->PutSharedOUTGate(temp);
  delete temp;
  temp = s_lambda;
  s_lambda = bc->PutSharedOUTGate(temp);
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
  prevErrNorm = s_prevErrNorm->get_clear_value<uint32_t>();
  delete s_prevErrNorm;
  lambda = s_lambda->get_clear_value<uint32_t>();
  delete s_lambda;

  collectTiming();
  collectCommunication();
  party->Reset();
  // arrays now contain secret shared values of plaintext from each party

  CLOCK(LM);
  TIC(LM);
  RunLMCircuit(threeDPts, y0, numPts, f, cx, cy, x, &prevErrNorm, &lambda, c,
               party, role);
  TOC(LM);

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
void BuildAndRunLMDO(share* s_threeDPts[], share* s_twoDPts[], int numPts,
                     share* s_f, share* s_cx, share* s_cy, share* s_x[],
                     BooleanCircuit* c, ABYParty* party, e_role role) {
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

  // LM Iteration
  float lambda = LM_LAMBDA_INIT;
  share* lambda_gate = c->PutCONSGate((uint32_t*)&lambda, bitlen);

  float f_max = std::numeric_limits<float>::max();
  share* prev_err_norm_gate = c->PutCONSGate((uint32_t*)&f_max, bitlen);

  for (int i = 0; i < LM_MAX_ITR; i++) {
    if (_verbosity & DBG_PROJECT) {
      cout << RED << "\n\nLM Iteration " << i << RESET << endl;
    }

    auto res =
        BuildLMIteration(s_threeDPts, s_y0, numPts, s_f, s_cx, s_cy, s_x,
                         prev_err_norm_gate, lambda_gate, c, party, role);
    prev_err_norm_gate = res.first;
    lambda_gate = res.second;
  }
}

void BuildAndRunLM(share* s_threeDPts[], share* s_twoDPts[], int numPts,
                   share* s_f, share* s_cx, share* s_cy, share* s_x[],
                   BooleanCircuit* c, ABYParty* party, e_role role) {
#if PPL_FLOW == PPL_FLOW_DO
  BuildAndRunLMDO(s_threeDPts, s_twoDPts, numPts, s_f, s_cx, s_cy, s_x, c,
                  party, role);
#elif PPL_FLOW == PPL_FLOW_LOOP_LEAK || PPL_FLOW == PPL_FLOW_SiSL
  BuildAndRunLMLoopLeak(s_threeDPts, s_twoDPts, numPts, s_f, s_cx, s_cy, s_x, c,
                        party, role);
#endif
}

uint32_t test_lm_circuit(
    e_role role, const std::string& address, uint16_t port, seclvl seclvl,
    uint32_t nthreads, e_mt_gen_alg mt_alg, e_sharing sharing,
    vector<cv::Point3f> threeDPts, vector<cv::Point2f> twoDPts, float f,
    float cx, float cy,
    float* x /* initial guess for { r1, r2, r3, t1, t2, t3 } */) {

  uint32_t reservegates = 65536;
  const std::string& abycircdir = "../../extern/ABY/bin/circ";
  ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                                 mt_alg, reservegates, abycircdir);
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();
#if PPL_FLOW == PPL_FLOW_LOOP_LEAK
  //  always use boolean shares (not yao) so BuildAndRunLM can get raw shares.
  //  this is because you cant use Y2B shares on Yao input gates.
  Circuit* bc = sharings[S_BOOL]->GetCircuitBuildRoutine();
#else
  //  if data oblivious, either yao or bool can be used
  Circuit* bc = circ;
#endif

  // sharings[S_BOOL]->SetPreCompPhaseValue(ePreCompRAMWrite);

  int numPts = threeDPts.size();

  // Allocate space for shares
  share** s_threeDPts = new share*[4 * numPts];
  share** s_twoDPts = new share*[3 * numPts];
  share* s_f;
  share* s_cx;
  share* s_cy;
  share** s_x = new share*[6];

  // Prepare inputs
  if (role == SERVER) {
    // float one=1.0;
    for (int p = 0; p < numPts; p++) {
      s_threeDPts[p] = bc->PutINGate((uint32_t*)&threeDPts[p].x, bitlen, role);
    }
    for (int p = 0; p < numPts; p++) {
      s_threeDPts[numPts + p] =
          bc->PutINGate((uint32_t*)&threeDPts[p].y, bitlen, role);
    }
    for (int p = 0; p < numPts; p++) {
      s_threeDPts[(2 * numPts) + p] =
          bc->PutINGate((uint32_t*)&threeDPts[p].z, bitlen, role);
    }
    // for(int p=0; p<numPts; p++) {
    //     s_threeDPts[3*numPts+p] = bc->PutINGate((uint32_t*) &one, bitlen,
    //     role);
    // }
    for (int p = 0; p < numPts; p++) {
      s_twoDPts[p] = bc->PutINGate((uint32_t*)&twoDPts[p].x, bitlen, role);
    }
    for (int p = 0; p < numPts; p++) {
      s_twoDPts[numPts + p] =
          bc->PutINGate((uint32_t*)&twoDPts[p].y, bitlen, role);
    }
    // for(int p=0; p<numPts; p++) {
    //     s_twoDPts[2*numPts+p] = bc->PutINGate((uint32_t*) &one, bitlen,
    //     role);
    // }
    s_f = bc->PutINGate((uint32_t*)&f, bitlen, role);
    s_cx = bc->PutINGate((uint32_t*)&cx, bitlen, role);
    s_cy = bc->PutINGate((uint32_t*)&cy, bitlen, role);
    for (int p = 0; p < 6; p++) {
      s_x[p] = bc->PutINGate((uint32_t*)&x[p], bitlen, role);
    }
  } else {
    for (int p = 0; p < 3 * numPts; p++) {  // ignore constant 1s
      s_threeDPts[p] = bc->PutDummyINGate(bitlen);
    }
    for (int p = 0; p < 2 * numPts; p++) {  // ignore constant 1s
      s_twoDPts[p] = bc->PutDummyINGate(bitlen);
    }
    s_f = bc->PutDummyINGate(bitlen);
    s_cx = bc->PutDummyINGate(bitlen);
    s_cy = bc->PutDummyINGate(bitlen);
    for (int p = 0; p < 6; p++) {
      s_x[p] = bc->PutDummyINGate(bitlen);
    }
  }

  BuildAndRunLM(s_threeDPts, s_twoDPts, numPts, s_f, s_cx, s_cy, s_x,
                (BooleanCircuit*)circ, party, role);

  for (int i = 0; i < 6; i++) {
    share* temp = s_x[i];
    s_x[i] = bc->PutOUTGate(s_x[i], ALL);
    delete temp;
  }

  party->ExecCircuit();

  for (int i = 0; i < 6; i++) {
    uint32_t* output;
    uint32_t out_bitlen, out_nvals;

    // This method only works for an output length of maximum 64 bits in
    // general, if the output length is higher you must use get_clear_value_ptr
    s_x[i]->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
    // s_out[i]->get_clear_value_ptr();

    // cout << *(float*)output << " ";

    x[i] = *(float*)output;
  }

  delete party;

  return 0;
}
