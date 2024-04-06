#include <jlog.h>
#include <privacyconf.h>
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

#include <invert.h>
#include <svd.h>

using namespace std;

const uint32_t bitlen = 32;

void BuildInvertCircuitWSubCircuits(share** s_in[], int m, int n,
                                    share** s_res[], BooleanCircuit* c,
                                    ABYParty* party, e_role role,
                                    std::function<void()> toRawShares,
                                    std::function<void()> toShareObjects) {

  assert(m >= n);
  float zero = 0;

  // svd overwrites with s_in with u matrix, must copy shares
  share*** s_u = new share**[m];
  for (int i = 0; i < m; i++) {
    s_u[i] = new share*[n];
    for (int j = 0; j < n; j++) {
      s_u[i][j] = s_in[i][j];
    }
  }

  share** s_w = new share*[m];    // aka sigma, only diag
  share*** s_v = new share**[n];  // nxn
  for (int i = 0; i < n; i++)
    s_v[i] = new share*[n];

  BuildAndRunSvd(s_u, m, n, s_w, s_v, (BooleanCircuit*)c, party, role,
                 toRawShares, toShareObjects);

  // make new zero gate after BuildAndRunSvd because it makes a new circuit
  share* zerogate = c->PutCONSGate((uint32_t*)&zero, bitlen);

  // res = inv(in) = v*inv(w)*uT
  // inv(w)*uT
  for (int j = 0; j < n; j++) {
    // if (s_w[j]) {
    share* ifw = s_w[j]->get_wire_ids_as_share(0);
    for (uint32_t tt = 1; tt < bitlen - 1;
         tt++) {  // -1 -> do not include sign bit
      share* tempbit = s_w[j]->get_wire_ids_as_share(tt);
      share* temp = ifw;
      ifw = c->PutORGate(temp, tempbit);
      delete tempbit;
      delete temp;
    }

    for (int i = 0; i < m; i++) {
      share* temp = c->PutFPGate(s_u[i][j], s_w[j], DIV, bitlen, 1, no_status);
      share* temp2 = s_u[i][j];
      s_u[i][j] = c->PutMUXGate(temp, zerogate, ifw);
      delete temp;
      // delete temp2; DO NOT DELETE s_u, it is still part of _s_in
    }
    delete ifw;
    //}
  }

  // v*(w*uT) (don't use matmult so we can do transpose ourselves)
  for (int j = 0; j < n; j++) {
    for (int jj = 0; jj < m; jj++) {
      // dont delete s_res[j][jj] - share may be elsewhere
      s_res[j][jj] = c->PutCONSGate((uint32_t*)&zero, bitlen);
      for (int k = 0; k < n; k++) {
        // note u indices do transpose
        share* temp =
            c->PutFPGate(s_v[j][k], s_u[jj][k], MUL, bitlen, 1, no_status);
        share* temp2 = s_res[j][jj];
        s_res[j][jj] =
            c->PutFPGate(s_res[j][jj], temp, ADD, bitlen, 1, no_status);
        delete temp;
        delete temp2;
      }
    }
  }
}

void BuildInvertCircuitDO(share** s_in[], int m, int n, share** s_res[],
                          BooleanCircuit* c, ABYParty* party, e_role role) {
  assert(m >= n);
  float zero = 0;

  share** s_w = new share*[m];    // aka sigma, only diag
  share*** s_v = new share**[n];  // nxn
  for (int i = 0; i < n; i++)
    s_v[i] = new share*[n];

  BuildAndRunSvd(s_in, m, n, s_w, s_v, (BooleanCircuit*)c, party, role);

  share*** s_u = s_in;
  share* zerogate = c->PutCONSGate((uint32_t*)&zero, bitlen);

  // res = inv(in) = v*inv(w)*uT
  // inv(w)*uT
  for (int j = 0; j < n; j++) {
    // if (s_w[j]) {
    share* ifw = s_w[j]->get_wire_ids_as_share(0);
    for (uint32_t tt = 1; tt < bitlen - 1;
         tt++) {  // -1 -> do not include sign bit
      share* tempbit = s_w[j]->get_wire_ids_as_share(tt);
      share* temp = ifw;
      ifw = c->PutORGate(temp, tempbit);
      delete tempbit;
      delete temp;
    }

    for (int i = 0; i < m; i++) {
      share* temp = c->PutFPGate(s_u[i][j], s_w[j], DIV, bitlen, 1, no_status);
      share* temp2 = s_u[i][j];
      s_u[i][j] = c->PutMUXGate(temp, zerogate, ifw);
      delete temp;
      // delete temp2; DO NOT DELETE s_u, it is still part of _s_in
    }
    delete ifw;
    //}
  }

  // v*(w*uT) (don't use matmult so we can do transpose ourselves)
  for (int j = 0; j < n; j++) {
    for (int jj = 0; jj < m; jj++) {
      // dont delete s_res[j][jj] - share may be elsewhere
      s_res[j][jj] = c->PutCONSGate((uint32_t*)&zero, bitlen);
      for (int k = 0; k < n; k++) {
        // note u indices do transpose
        share* temp =
            c->PutFPGate(s_v[j][k], s_u[jj][k], MUL, bitlen, 1, no_status);
        share* temp2 = s_res[j][jj];
        s_res[j][jj] =
            c->PutFPGate(s_res[j][jj], temp, ADD, bitlen, 1, no_status);
        delete temp;
        delete temp2;
      }
    }
  }
}

void BuildInvertCircuit(share** s_in[], int m, int n, share** s_res[],
                        BooleanCircuit* c, ABYParty* party, e_role role,
                        std::function<void()> toRawShares,
                        std::function<void()> toShareObjects) {
#if PPL_FLOW == PPL_FLOW_DO || \
    PPL_FLOW == PPL_FLOW_SiSL  // set dx to zero if no error
  (void)toRawShares;
  (void)toShareObjects;
  BuildInvertCircuitDO(s_in, m, n, s_res, c, party, role);
#elif PPL_FLOW == PPL_FLOW_LOOP_LEAK
  BuildInvertCircuitWSubCircuits(s_in, m, n, s_res, c, party, role, toRawShares,
                                 toShareObjects);
#endif
}
