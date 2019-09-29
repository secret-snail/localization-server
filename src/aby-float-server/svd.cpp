#include <jlog.h>
#include <privacyconf.h>
#include <stdio.h>
#include <svd.h>
#include <iostream>
#include <string>

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/circuit/share.h"
#include "abycore/sharing/sharing.h"

#include <util.h>

#include <math.h>

#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a))

static float maxarg1, maxarg2;
#define FMAX(a, b) \
  (maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1, iminarg2;
#define IMIN(a, b)                 \
  (iminarg1 = (a), iminarg2 = (b), \
   (iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

static float sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

using namespace std;

const uint32_t bitlen = 32;

//#define TWODINDEX(mat,ncols,a,b) (mat[(a*ncols)+b])

void BuildFabsCircuit(share* a, BooleanCircuit* c) {
  // sets first bit of a to zero (can be public)
  uint32_t zero = 0;
  share* zero_gate = c->PutCONSGate((uint32_t*)&zero, 32);
  a->set_wire_id(31, zero_gate->get_wire_id(31));
}

void BuildSignCircuit(share* a, share* b, BooleanCircuit* c) {
  // WARNING:
  // cannot simply set msb of a to msb of b because
  // of the case when b is 0.0
  // TODO - fix memory leak
  uint32_t zero = 0;
  share* zero_gate = c->PutCONSGate((uint32_t*)&zero, 32);
  share* cmp = c->PutFPGate(b, zero_gate, CMP);
  share* tmp = c->PutINVGate(cmp);
  a->set_wire_id(31, tmp->get_wire_id(0));
}

share* BuildFMAXCircuit(share* a, share* b, BooleanCircuit* c) {
  // returns which a or b is bigger
  share* cmp = c->PutFPGate(a, b, CMP);
  share* ret = c->PutMUXGate(a, b, cmp);
  delete cmp;
  return ret;
}

void BuildNegativeCircuit(share* a, BooleanCircuit* c) {
  // flip msb
  // TODO - fix memory leak
  share* tmp = a->get_wire_ids_as_share(31);
  tmp = c->PutINVGate(tmp);
  // cout << tmp->get_bitlength() << endl; // tmp is one bit
  a->set_wire_id(31, tmp->get_wire_id(0));
}

share* BuildPythagCircuit(share* a, share* b, BooleanCircuit* c,
                          uint32_t bitlen) {
  share* aa = c->PutFPGate(a, a, MUL, bitlen, 1, no_status);
  share* bb = c->PutFPGate(b, b, MUL, bitlen, 1, no_status);
  share* sum = c->PutFPGate(aa, bb, ADD, bitlen, 1, no_status);
  share* res = c->PutFPGate(sum, SQRT);
  delete aa;
  delete bb;
  delete sum;
  return res;

  // Doing the same pythag optimization as cleartext is
  // actually slower than the long way to compute pythag
  // float zero = 0;
  // share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);
  // float one = 1;
  // share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);

  // boolshare *acopy = new boolshare(a->get_wires(), c);//copy to new share
  // boolshare *bcopy = new boolshare(b->get_wires(), c);//copy to new share
  // BuildFabsCircuit(acopy, c);
  // BuildFabsCircuit(bcopy, c);

  // share *adb = c->PutFPGate(acopy, bcopy, DIV, bitlen, 1, no_status);
  // share *sqadb = c->PutFPGate(adb, adb, MUL, bitlen, 1, no_status);
  // share *addone = c->PutFPGate(sqadb, one_gate, ADD, bitlen, 1, no_status);
  // share *sqrt = c->PutFPGate(addone, SQRT);

  // share *mul = BuildFMAXCircuit(acopy, bcopy, c);

  // share *res = c->PutFPGate(mul, sqrt, MUL, bitlen, 1, no_status);

  // share *cmp = c->PutFPGate(mul, zero_gate, CMP); // mul > 0
  // share *resorzero = c->PutMUXGate(res, zero_gate, cmp);

  // delete zero_gate;
  // delete one_gate;
  // delete acopy;
  // delete bcopy;
  // delete adb;
  // delete sqadb;
  // delete addone;
  // delete sqrt;
  // delete mul;
  // delete cmp;
  // delete res;
  // return resorzero;
}

// r=3x1, R=3x4 (only 3x3 is used here)
// values are secret, sizes are public
void BuildSvdPart1Circuit(share** a[], int nRows, int nCols, share* w[],
                          share** v[], share** anorm, share* rv1[],
                          BooleanCircuit* c) {

  float zero = 0;
  share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);
  float one = 1;
  share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);
  share* temp;

  int i, j, k, l;  // local plaintext
  share *f, *g, *h, *s,
      *scale;  // local secrets (are reset before being used elsewhere)

  // g = scale = anorm = zero_gate;
  g = new boolshare(zero_gate->get_wires(), c);       // copy to new share
  scale = new boolshare(zero_gate->get_wires(), c);   // copy to new share
  *anorm = new boolshare(zero_gate->get_wires(), c);  // copy to new share

  for (i = 0; i < nCols; i++) {
    l = i + 1;
    rv1[i] = c->PutFPGate(scale, g, MUL, bitlen, 1, no_status);
    g = new boolshare(zero_gate->get_wires(), c);      // copy to new share
    s = new boolshare(zero_gate->get_wires(), c);      // copy to new share
    scale = new boolshare(zero_gate->get_wires(), c);  // copy to new share
    if (i < nRows) {
      for (k = i; k < nRows; k++) {
        boolshare* shcopy =
            new boolshare(a[k][i]->get_wires(), c);  // copy to new share
        BuildFabsCircuit(shcopy, c);
        scale = c->PutFPGate(scale, shcopy, ADD, bitlen, 1, no_status);
      }
      {  // if (scale)
        share* ifscale = scale->get_wire_ids_as_share(0);
        for (uint32_t tt = 1; tt < bitlen - 1;
             tt++) {  // -1 -> do not include sign bit
          share* tempbit = scale->get_wire_ids_as_share(tt);
          temp = ifscale;
          ifscale = c->PutORGate(temp, tempbit);
          delete tempbit;
          delete temp;
        }

        for (k = i; k < nRows; k++) {
          temp = c->PutFPGate(a[k][i], scale, DIV, bitlen, 1, no_status);
          a[k][i] = c->PutMUXGate(temp, a[k][i], ifscale);
          temp = c->PutFPGate(a[k][i], SQR, bitlen, 1, no_status);
          s = c->PutFPGate(s, temp, ADD, bitlen, 1, no_status);
        }
        f = new boolshare(a[i][i]->get_wires(), c);  // no mux okay
        temp = c->PutFPGate(s, SQRT);
        BuildSignCircuit(temp, f, c);
        BuildNegativeCircuit(temp, c);
        g = c->PutMUXGate(temp, g, ifscale);
        temp = c->PutFPGate(f, g, MUL, bitlen, 1, no_status);
        h = c->PutFPGate(temp, s, SUB, bitlen, 1, no_status);  // no mux okay
        temp = c->PutFPGate(f, g, SUB, bitlen, 1, no_status);
        a[i][i] = c->PutMUXGate(temp, a[i][i], ifscale);
        for (j = l; j < nCols; j++) {
          s = new boolshare(zero_gate->get_wires(), c);
          for (k = i; k < nRows; k++) {
            temp = c->PutFPGate(a[k][i], a[k][j], MUL, bitlen, 1, no_status);
            s = c->PutFPGate(s, temp, ADD, bitlen, 1, no_status);
          }
          f = c->PutFPGate(s, h, DIV, bitlen, 1, no_status);
          for (k = i; k < nRows; k++) {
            temp = c->PutFPGate(f, a[k][i], MUL, bitlen, 1, no_status);
            temp = c->PutFPGate(a[k][j], temp, ADD, bitlen, 1, no_status);
            a[k][j] = c->PutMUXGate(temp, a[k][j], ifscale);
          }
        }
        for (k = i; k < nRows; k++) {
          temp = c->PutFPGate(a[k][i], scale, MUL, bitlen, 1, no_status);
          a[k][i] = c->PutMUXGate(temp, a[k][i], ifscale);
        }
      }
    }
    w[i] = c->PutFPGate(scale, g, MUL, bitlen, 1, no_status);
    g = new boolshare(zero_gate->get_wires(), c);      // copy to new share
    s = new boolshare(zero_gate->get_wires(), c);      // copy to new share
    scale = new boolshare(zero_gate->get_wires(), c);  // copy to new share
    if (i < nRows && i != nCols - 1) {
      for (k = l; k < nCols; k++) {
        boolshare* shcopy =
            new boolshare(a[i][k]->get_wires(), c);  // copy to new share
        BuildFabsCircuit(shcopy, c);
        scale = c->PutFPGate(scale, shcopy, ADD, bitlen, 1, no_status);
      }
      {  // if (scale)
        share* ifscale = scale->get_wire_ids_as_share(0);
        for (uint32_t tt = 1; tt < bitlen - 1;
             tt++) {  // -1 -> do not include sign bit
          share* tempbit = scale->get_wire_ids_as_share(tt);
          temp = ifscale;
          ifscale = c->PutORGate(temp, tempbit);
          delete tempbit;
          delete temp;
        }

        for (k = l; k < nCols; k++) {
          a[i][k] = c->PutFPGate(a[i][k], scale, DIV, bitlen, 1, no_status);
          temp = c->PutFPGate(a[i][k], SQR);  //, bitlen, 1, no_status);
          s = c->PutFPGate(s, temp, ADD, bitlen, 1, no_status);
        }
        f = new boolshare(a[i][l]->get_wires(),
                          c);  // not used again, no mux okay
        temp = c->PutFPGate(s, SQRT);
        BuildSignCircuit(temp, f, c);
        BuildNegativeCircuit(temp, c);
        g = c->PutMUXGate(temp, g, ifscale);
        temp = c->PutFPGate(f, g, MUL, bitlen, 1, no_status);
        h = c->PutFPGate(temp, s, SUB, bitlen, 1, no_status);
        share* htoone = new boolshare(one_gate->get_wires(), c);
        h = c->PutMUXGate(
            h, htoone,
            ifscale);  // h is used to set rv1 if scale==0, set to 1 for noop
        temp = c->PutFPGate(f, g, SUB, bitlen, 1, no_status);
        a[i][l] = c->PutMUXGate(temp, a[i][l], ifscale);
        for (k = l; k < nCols; k++) {
          rv1[k] = c->PutFPGate(a[i][k], h, DIV, bitlen, 1,
                                no_status);  // no mux because h is one
        }
        for (j = l; j < nRows; j++) {
          s = new boolshare(zero_gate->get_wires(), c);
          for (k = l; k < nCols; k++) {
            temp = c->PutFPGate(a[j][k], a[i][k], MUL, bitlen, 1, no_status);
            s = c->PutFPGate(s, temp, ADD, bitlen, 1, no_status);
          }
          for (k = l; k < nCols; k++) {
            temp = c->PutFPGate(s, rv1[k], MUL, bitlen, 1, no_status);
            temp = c->PutFPGate(a[j][k], temp, ADD, bitlen, 1, no_status);
            a[j][k] = c->PutMUXGate(temp, a[j][k], ifscale);
          }
        }
        for (k = l; k < nCols; k++) {
          temp = c->PutFPGate(a[i][k], scale, MUL, bitlen, 1, no_status);
          a[i][k] = c->PutMUXGate(temp, a[i][k], ifscale);
        }
      }
    }
    boolshare* shcopy1 =
        new boolshare(w[i]->get_wires(), c);  // copy to new share
    boolshare* shcopy2 =
        new boolshare(rv1[i]->get_wires(), c);  // copy to new share

    BuildFabsCircuit(shcopy1, c);
    BuildFabsCircuit(shcopy2, c);

    temp = c->PutFPGate(shcopy1, shcopy2, ADD, bitlen, 1, no_status);
    *anorm = BuildFMAXCircuit(*anorm, temp, c);
    printf(".");
    fflush(stdout);
  }

  for (i = nCols - 1; i >= 0; i--) {
    if (i < nCols - 1) {
      {  // if (g)
        share* ifg = g->get_wire_ids_as_share(0);
        for (uint32_t tt = 1; tt < bitlen - 1;
             tt++) {  // -1 -> do not include sign bit
          share* tempbit = g->get_wire_ids_as_share(tt);
          temp = ifg;
          ifg = c->PutORGate(temp, tempbit);
          delete tempbit;
          delete temp;
        }
        for (j = l; j < nCols; j++) {
          temp = c->PutFPGate(a[i][j], a[i][l], DIV, bitlen, 1, no_status);
          temp = c->PutFPGate(temp, g, DIV, bitlen, 1, no_status);
          v[j][i] = c->PutMUXGate(temp, g, ifg);
        }
        for (j = l; j < nCols; j++) {
          s = new boolshare(zero_gate->get_wires(), c);
          for (k = l; k < nCols; k++) {
            temp = c->PutFPGate(a[i][k], v[k][j], MUL, bitlen, 1, no_status);
            s = c->PutFPGate(s, temp, ADD, bitlen, 1,
                             no_status);  // s only used here so no mux okay
          }
          for (k = l; k < nCols; k++) {
            temp = c->PutFPGate(s, v[k][i], MUL, bitlen, 1, no_status);
            temp = c->PutFPGate(v[k][j], temp, ADD, bitlen, 1, no_status);
            v[k][j] = c->PutMUXGate(temp, v[k][j], ifg);
          }
        }
      }
      for (j = l; j < nCols; j++) {
        v[i][j] = new boolshare(zero_gate->get_wires(), c);
        v[j][i] = new boolshare(zero_gate->get_wires(), c);
      }
    }
    v[i][i] = new boolshare(one_gate->get_wires(), c);
    g = new boolshare(rv1[i]->get_wires(), c);
    l = i;
    printf(":");
    fflush(stdout);
  }

  for (i = IMIN(nRows, nCols) - 1; i >= 0; i--) {
    l = i + 1;
    g = new boolshare(w[i]->get_wires(), c);
    for (j = l; j < nCols; j++)
      a[i][j] = new boolshare(zero_gate->get_wires(), c);
    {  // if (g)
      share* ifg = g->get_wire_ids_as_share(0);
      for (uint32_t tt = 1; tt < bitlen - 1;
           tt++) {  // -1 -> do not include sign bit
        share* tempbit = g->get_wire_ids_as_share(tt);
        temp = ifg;
        ifg = c->PutORGate(temp, tempbit);
        delete tempbit;
        delete temp;
      }
      g = c->PutFPGate(one_gate, g, DIV, bitlen, 1, no_status);
      for (j = l; j < nCols; j++) {
        s = new boolshare(zero_gate->get_wires(), c);
        for (k = l; k < nRows; k++) {
          temp = c->PutFPGate(a[k][i], a[k][j], MUL, bitlen, 1, no_status);
          s = c->PutFPGate(
              s, temp, ADD, bitlen, 1,
              no_status);  // s is only used here so no mux gate okay
        }
        temp = c->PutFPGate(s, a[i][i], DIV, bitlen, 1, no_status);
        f = c->PutFPGate(temp, g, MUL, bitlen, 1,
                         no_status);  // f is only used here so no mux gate okay
        for (k = i; k < nRows; k++) {
          temp = c->PutFPGate(f, a[k][i], MUL, bitlen, 1, no_status);
          temp = c->PutFPGate(a[k][j], temp, ADD, bitlen, 1, no_status);
          a[k][j] =
              c->PutMUXGate(temp, a[k][j], ifg);  // do nothing if(g) is false
        }
      }
      for (j = i; j < nRows; j++) {
        temp = c->PutFPGate(a[j][i], g, MUL, bitlen, 1, no_status);
        share* zcpy = new boolshare(zero_gate->get_wires(), c);
        a[j][i] = c->PutMUXGate(temp, zcpy, ifg);  // set to zero ifg is false
        delete zcpy;
        delete temp;
      }
    }
    //} else  - done in for loop above with mux
    //    for (j = i; j < nRows; j++)
    //        a[j][i] = 0.0;
    a[i][i] = c->PutFPGate(a[i][i], one_gate, ADD, bitlen, 1, no_status);
    printf("|");
    fflush(stdout);
  }
}

void RunSvdPart1Circuit(uint32_t** a, int nRows, int nCols, uint32_t* w,
                        uint32_t** v, uint32_t* anorm, uint32_t* rv1,
                        BooleanCircuit* c, ABYParty* p) {

  std::vector<Sharing*>& sharings = p->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  // Allocate space for shares
  share* temp;
  share*** s_a = new share**[nRows];
  for (int i = 0; i < nRows; i++) {
    s_a[i] = new share*[nCols];
  }
  share** s_w = new share*[nCols];
  share*** s_v = new share**[nCols];
  for (int i = 0; i < nCols; i++) {
    s_v[i] = new share*[nCols];
  }
  share* s_anorm;
  share** s_rv1 = new share*[nCols];

  // Prepare inputs
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      s_a[i][j] = bc->PutSharedINGate(&a[i][j], bitlen);
    }
  }
  // if yao, convert bool inputs to yao
  if (c->GetContext() == S_YAO) {
    for (int i = 0; i < nRows; i++) {
      for (int j = 0; j < nCols; j++) {
        temp = s_a[i][j];
        s_a[i][j] = c->PutB2YGate(temp);
        delete temp;
      }
    }
  }
  // prepare circuit
  BuildSvdPart1Circuit(s_a, nRows, nCols, s_w, s_v, &s_anorm, s_rv1, c);

  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    for (int t = 0; t < nRows; t++) {
      for (int tt = 0; tt < nCols; tt++) {
        temp = s_a[t][tt];
        s_a[t][tt] = bc->PutY2BGate(temp);
        delete temp;
      }
    }
    for (int t = 0; t < nCols; t++) {
      temp = s_w[t];
      s_w[t] = bc->PutY2BGate(temp);
      delete temp;
    }
    for (int t = 0; t < nCols; t++) {
      for (int tt = 0; tt < nCols; tt++) {
        temp = s_v[t][tt];
        s_v[t][tt] = bc->PutY2BGate(temp);
        delete temp;
      }
    }
    temp = s_anorm;
    s_anorm = bc->PutY2BGate(temp);
    delete temp;
    for (int t = 0; t < nCols; t++) {
      temp = s_rv1[t];
      s_rv1[t] = bc->PutY2BGate(temp);
      delete temp;
    }
  }

  // prepare output
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      temp = s_a[i][j];
      s_a[i][j] = bc->PutSharedOUTGate(temp);
      delete temp;
    }
  }
  for (int i = 0; i < nCols; i++) {
    temp = s_w[i];
    s_w[i] = bc->PutSharedOUTGate(temp);
    delete temp;
  }
  for (int i = 0; i < nCols; i++) {
    for (int j = 0; j < nCols; j++) {
      temp = s_v[i][j];
      s_v[i][j] = bc->PutSharedOUTGate(temp);
      delete temp;
    }
  }

  temp = s_anorm;
  s_anorm = bc->PutSharedOUTGate(temp);
  delete temp;
  for (int i = 0; i < nCols; i++) {
    temp = s_rv1[i];
    s_rv1[i] = bc->PutSharedOUTGate(temp);
    delete temp;
  }

  p->ExecCircuit();

  // get shared output
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      a[i][j] = s_a[i][j]->get_clear_value<uint32_t>();
      delete s_a[i][j];
    }
    delete[] s_a[i];
  }
  delete[] s_a;
  for (int i = 0; i < nCols; i++) {
    w[i] = s_w[i]->get_clear_value<uint32_t>();
    delete s_w[i];
  }
  delete[] s_w;
  for (int i = 0; i < nCols; i++) {
    for (int j = 0; j < nCols; j++) {
      v[i][j] = s_v[i][j]->get_clear_value<uint32_t>();
      delete s_v[i][j];
    }
    delete[] s_v[i];
  }
  delete[] s_v;

  *anorm = s_anorm->get_clear_value<uint32_t>();
  delete s_anorm;

  for (int i = 0; i < nCols; i++) {
    rv1[i] = s_rv1[i]->get_clear_value<uint32_t>();
    delete s_rv1[i];
  }
  collectTiming();
  collectCommunication();
  p->Reset();
}

share* BuildSvdCheckZeroCircuit(share* val, share* anorm, BooleanCircuit* c) {
  // This evaluates:
  //  ((fabs(val) + anorm) == anorm)
  BuildFabsCircuit(val,
                   c);  // overwriting is okay because this is separate circuit
  share* vpa = c->PutFPGate(val, anorm, ADD, bitlen, 1, no_status);
  share* gt =
      c->PutFPGate(vpa, anorm, CMP);  // = val+anorm > anorm (val always > 0)
  share* res = c->PutINVGate(gt);
  delete vpa;
  delete gt;
  return res;
}

uint32_t RunSvdCheckZeroCircuit(uint32_t* val, uint32_t* anorm,
                                BooleanCircuit* c, ABYParty* p) {

  std::vector<Sharing*>& sharings = p->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  // Allocate space for shares
  share* s_val;
  share* s_anorm;
  // Prepare inputs
  s_val = bc->PutSharedINGate(val, bitlen);
  s_anorm = bc->PutSharedINGate(anorm, bitlen);
  // if yao, convert bool inputs to yao
  if (c->GetContext() == S_YAO) {
    share* temp = s_val;
    s_val = c->PutB2YGate(temp);
    delete temp;
    temp = s_anorm;
    s_anorm = c->PutB2YGate(temp);
    delete temp;
  }
  // prepare circuit
  share* s_res = BuildSvdCheckZeroCircuit(s_val, s_anorm, c);
  // prepare CLEARTEXT output
  share* temp = s_res;
  s_res = c->PutOUTGate(temp, ALL);
  delete temp;

  p->ExecCircuit();

  // get shared output
  uint32_t res = s_res->get_clear_value<uint32_t>();
  delete s_res;
  delete s_val;
  delete s_anorm;

  collectTiming();
  collectCommunication();
  p->Reset();
  return res;
}

void BuildSvdFlag1Circuit(share** cc, share** f, share** s, share** rv1i,
                          BooleanCircuit* c) {
  *f = c->PutFPGate(*s, *rv1i, MUL, bitlen, 1, no_status);
  *rv1i = c->PutFPGate(*cc, *rv1i, MUL, bitlen, 1, no_status);
}
void RunSvdFlag1Circuit(uint32_t* cc, uint32_t* f, uint32_t* s, uint32_t* rv1i,
                        BooleanCircuit* c, ABYParty* p) {

  std::vector<Sharing*>& sharings = p->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  // Allocate space for shares
  share* temp;
  share* s_cc;
  share* s_f;
  share* s_s;
  share* s_rv1i;
  // Prepare inputs
  s_cc = bc->PutSharedINGate(cc, bitlen);
  s_f = bc->PutSharedINGate(f, bitlen);
  s_s = bc->PutSharedINGate(s, bitlen);
  s_rv1i = bc->PutSharedINGate(rv1i, bitlen);
  // if yao, convert bool inputs to yao
  if (c->GetContext() == S_YAO) {
    temp = s_cc;
    s_cc = c->PutB2YGate(temp);
    delete temp;
    temp = s_f;
    s_f = c->PutB2YGate(temp);
    delete temp;
    temp = s_s;
    s_s = c->PutB2YGate(temp);
    delete temp;
    temp = s_rv1i;
    s_rv1i = c->PutB2YGate(temp);
    delete temp;
  }
  // prepare circuit
  BuildSvdFlag1Circuit(&s_cc, &s_f, &s_s, &s_rv1i, c);
  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    temp = s_f;
    s_f = bc->PutY2BGate(temp);
    delete temp;
    temp = s_rv1i;
    s_rv1i = bc->PutY2BGate(temp);
    delete temp;
  }
  // prepare output
  temp = s_f;
  s_f = bc->PutSharedOUTGate(temp);
  delete temp;
  temp = s_rv1i;
  s_rv1i = bc->PutSharedOUTGate(temp);
  delete temp;

  p->ExecCircuit();

  // get shared output
  *f = s_f->get_clear_value<uint32_t>();
  *rv1i = s_rv1i->get_clear_value<uint32_t>();
  delete s_cc;
  delete s_f;
  delete s_s;
  delete s_rv1i;

  collectTiming();
  collectCommunication();
  p->Reset();
}

void BuildSvdFlag2Circuit(share* aSTARi[], share* aSTARnm[], int nRows,
                          share** wi, share** cc, share** f, share** s,
                          BooleanCircuit* c) {
  // non secret requires: i, nm
  // secret requires: a[*][i], a[*][nm], w[i], cc, f
  share *g, *h, *y, *z;

  float one = 1.0;
  share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);

  g = new boolshare((*wi)->get_wires(), c);
  h = BuildPythagCircuit(*f, g, c, bitlen);
  *wi = new boolshare(h->get_wires(), c);
  h = c->PutFPGate(one_gate, h, DIV, bitlen, 1, no_status);
  *cc = c->PutFPGate(g, h, MUL, bitlen, 1, no_status);
  share* temp = c->PutFPGate(*f, h, MUL, bitlen, 1, no_status);
  BuildNegativeCircuit(temp, c);
  *s = temp;
  for (int j = 0; j < nRows; j++) {
    y = aSTARnm[j];  // y = a[j][nm];
    z = aSTARi[j];   // z = a[j][i];
    share* tmp1 = c->PutFPGate(y, *cc, MUL, bitlen, 1, no_status);
    share* tmp2 = c->PutFPGate(z, *s, MUL, bitlen, 1, no_status);
    aSTARnm[j] = c->PutFPGate(tmp1, tmp2, ADD, bitlen, 1, no_status);
    // a[j][nm]
    delete tmp1;
    delete tmp2;

    tmp1 = c->PutFPGate(z, *cc, MUL, bitlen, 1, no_status);
    tmp2 = c->PutFPGate(y, *s, MUL, bitlen, 1, no_status);
    aSTARi[j] = c->PutFPGate(tmp1, tmp2, SUB, bitlen, 1, no_status);
    // a[j][i]
    delete tmp1;
    delete tmp2;
  }
}
void RunSvdFlag2Circuit(uint32_t* aSTARi, uint32_t* aSTARnm, int nRows,
                        uint32_t* wi, uint32_t* cc, uint32_t* f, uint32_t* s,
                        BooleanCircuit* c, ABYParty* p) {
  std::vector<Sharing*>& sharings = p->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  // Allocate space for shares
  share* temp;
  share** s_aSTARi = new share*[nRows];
  share** s_aSTARnm = new share*[nRows];
  share* s_wi;
  share* s_cc;
  share* s_f;
  share* s_s;
  // Prepare inputs
  for (int i = 0; i < nRows; i++) {
    s_aSTARi[i] = bc->PutSharedINGate(&aSTARi[i], bitlen);
    s_aSTARnm[i] = bc->PutSharedINGate(&aSTARnm[i], bitlen);
  }
  s_wi = bc->PutSharedINGate(wi, bitlen);
  s_cc = bc->PutSharedINGate(cc, bitlen);
  s_f = bc->PutSharedINGate(f, bitlen);
  s_s = bc->PutSharedINGate(s, bitlen);
  // if yao, convert bool inputs to yao
  if (c->GetContext() == S_YAO) {
    for (int i = 0; i < nRows; i++) {
      temp = s_aSTARi[i];
      s_aSTARi[i] = c->PutB2YGate(temp);
      delete temp;
      temp = s_aSTARnm[i];
      s_aSTARnm[i] = c->PutB2YGate(temp);
      delete temp;
    }
    temp = s_wi;
    s_wi = c->PutB2YGate(temp);
    delete temp;
    temp = s_cc;
    s_cc = c->PutB2YGate(temp);
    delete temp;
    temp = s_f;
    s_f = c->PutB2YGate(temp);
    delete temp;
    temp = s_s;
    s_s = c->PutB2YGate(temp);
    delete temp;
  }
  // prepare circuit
  BuildSvdFlag2Circuit(s_aSTARi, s_aSTARnm, nRows, &s_wi, &s_cc, &s_f, &s_s, c);
  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    for (int i = 0; i < nRows; i++) {
      temp = s_aSTARi[i];
      s_aSTARi[i] = bc->PutY2BGate(temp);
      delete temp;
      temp = s_aSTARnm[i];
      s_aSTARnm[i] = bc->PutY2BGate(temp);
      delete temp;
    }
    temp = s_wi;
    s_wi = bc->PutY2BGate(temp);
    delete temp;
    temp = s_cc;
    s_cc = bc->PutY2BGate(temp);
    delete temp;
    temp = s_f;
    s_f = bc->PutY2BGate(temp);
    delete temp;
    temp = s_s;
    s_s = bc->PutY2BGate(temp);
    delete temp;
  }
  // prepare output (everything except f which is not an output)
  for (int star = 0; star < nRows; star++) {
    temp = s_aSTARi[star];
    s_aSTARi[star] = bc->PutSharedOUTGate(temp);
    delete temp;
    temp = s_aSTARnm[star];
    s_aSTARnm[star] = bc->PutSharedOUTGate(temp);
    delete temp;
  }
  temp = s_wi;
  s_wi = bc->PutSharedOUTGate(temp);
  delete temp;
  temp = s_cc;
  s_cc = bc->PutSharedOUTGate(temp);
  delete temp;
  temp = s_s;
  s_s = bc->PutSharedOUTGate(temp);
  delete temp;

  p->ExecCircuit();

  // get shared output
  for (int star = 0; star < nRows; star++) {
    aSTARi[star] = s_aSTARi[star]->get_clear_value<uint32_t>();
    delete s_aSTARi[star];
    aSTARnm[star] = s_aSTARnm[star]->get_clear_value<uint32_t>();
    delete s_aSTARnm[star];
  }
  delete[] s_aSTARi;
  delete[] s_aSTARnm;
  *wi = s_wi->get_clear_value<uint32_t>();
  delete s_wi;
  *cc = s_cc->get_clear_value<uint32_t>();
  delete s_cc;
  delete s_f;
  *s = s_s->get_clear_value<uint32_t>();
  delete s_s;

  collectTiming();
  collectCommunication();
  p->Reset();
}

void BuildSvdZCircuit(int nCols, share** z, share** wk, share* vSTARk[],
                      BooleanCircuit* c) {
  // secret requires: z, w[k], v[*][k]

  float zero = 0;
  share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);

  // if (z < 0.0)
  share* isneg = c->PutFPGate(zero_gate, *z, CMP);
  share* temp = new boolshare((*z)->get_wires(), c);
  BuildNegativeCircuit(temp, c);
  *wk = c->PutMUXGate(temp, *wk, isneg);
  for (int j = 0; j < nCols; j++) {
    temp = new boolshare(vSTARk[j]->get_wires(), c);
    BuildNegativeCircuit(temp, c);
    vSTARk[j] = c->PutMUXGate(temp, vSTARk[j], isneg);  // v[j][k]
  }
}
void RunSvdZCircuit(int nCols, uint32_t* z, uint32_t* wk, uint32_t* vSTARk,
                    BooleanCircuit* c, ABYParty* p) {

  std::vector<Sharing*>& sharings = p->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  // Allocate space for shares
  share* temp;
  share* s_z;
  share* s_wk;
  share** s_vSTARk = new share*[nCols];
  // Prepare inputs
  s_z = bc->PutSharedINGate(z, bitlen);
  s_wk = bc->PutSharedINGate(wk, bitlen);
  for (int i = 0; i < nCols; i++) {
    s_vSTARk[i] = bc->PutSharedINGate(&vSTARk[i], bitlen);
  }
  // if yao, convert bool inputs to yao
  if (c->GetContext() == S_YAO) {
    temp = s_z;
    s_z = c->PutB2YGate(temp);
    delete temp;
    temp = s_wk;
    s_wk = c->PutB2YGate(temp);
    delete temp;
    for (int i = 0; i < nCols; i++) {
      temp = s_vSTARk[i];
      s_vSTARk[i] = c->PutB2YGate(temp);
      delete temp;
    }
  }
  // prepare circuit
  BuildSvdZCircuit(nCols, &s_z, &s_wk, s_vSTARk, c);
  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    temp = s_wk;
    s_wk = bc->PutY2BGate(temp);
    delete temp;
    for (int i = 0; i < nCols; i++) {
      temp = s_vSTARk[i];
      s_vSTARk[i] = bc->PutY2BGate(temp);
      delete temp;
    }
  }
  // prepare output (everything except z which is not an output)
  for (int star = 0; star < nCols; star++) {
    temp = s_vSTARk[star];
    s_vSTARk[star] = bc->PutSharedOUTGate(temp);
    delete temp;
  }
  temp = s_wk;
  s_wk = bc->PutSharedOUTGate(temp);
  delete temp;

  p->ExecCircuit();

  // get shared output
  delete s_z;
  *wk = s_wk->get_clear_value<uint32_t>();
  delete s_wk;
  for (int star = 0; star < nCols; star++) {
    vSTARk[star] = s_vSTARk[star]->get_clear_value<uint32_t>();
    delete s_vSTARk[star];
  }
  delete[] s_vSTARk;

  collectTiming();
  collectCommunication();
  p->Reset();
}

void BuildSvdPart2Circuit(share** a[], int nRows, int nCols, share* w[],
                          share** v[], share* rv1[], share** z, int k, int l,
                          int nm, BooleanCircuit* c) {
  // non secret requires: k, l, nm
  // secret requires: a, w, v, z, rv1 (all elements)

  float zero = 0;
  share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);
  float one = 1;
  share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);
  float two = 2;
  share* two_gate = c->PutCONSGate((uint32_t*)&two, bitlen);

  int jj;  // local plaintext (not used anywhere else)
  share *cc, *s, *x, *y, *g, *h, *f;

  x = new boolshare(w[l]->get_wires(), c);
  y = new boolshare(w[nm]->get_wires(), c);
  g = new boolshare(rv1[nm]->get_wires(), c);
  h = new boolshare(rv1[k]->get_wires(), c);
  {  // f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
    share* ymz = c->PutFPGate(y, *z, SUB, bitlen, 1, no_status);
    share* ypz = c->PutFPGate(y, *z, ADD, bitlen, 1, no_status);
    share* temp1 = c->PutFPGate(ymz, ypz, MUL, bitlen, 1, no_status);
    share* gmh = c->PutFPGate(g, h, SUB, bitlen, 1, no_status);
    share* gph = c->PutFPGate(g, h, ADD, bitlen, 1, no_status);
    share* temp2 = c->PutFPGate(gmh, gph, MUL, bitlen, 1, no_status);
    share* num = c->PutFPGate(temp1, temp2, ADD, bitlen, 1, no_status);
    share* den = c->PutFPGate(two_gate, h, MUL, bitlen, 1, no_status);
    den = c->PutFPGate(den, y, MUL, bitlen, 1, no_status);
    f = c->PutFPGate(num, den, DIV, bitlen, 1, no_status);
  }
  g = BuildPythagCircuit(f, one_gate, c, bitlen);
  {  // f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f)))- h)) / x;
    share* xmz = c->PutFPGate(x, *z, SUB, bitlen, 1, no_status);
    share* xpz = c->PutFPGate(x, *z, ADD, bitlen, 1, no_status);
    share* term1 = c->PutFPGate(xmz, xpz, MUL, bitlen, 1, no_status);
    share* sgf = new boolshare(g->get_wires(), c);
    BuildSignCircuit(sgf, f, c);
    share* fpsgf = c->PutFPGate(f, sgf, ADD, bitlen, 1, no_status);
    share* ydfpsgf = c->PutFPGate(y, fpsgf, DIV, bitlen, 1, no_status);
    share* uggg = c->PutFPGate(ydfpsgf, h, SUB, bitlen, 1, no_status);
    share* term2 = c->PutFPGate(h, uggg, MUL, bitlen, 1, no_status);
    share* num2 = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
    f = c->PutFPGate(num2, x, DIV, bitlen, 1, no_status);
  }
  cc = new boolshare(one_gate->get_wires(), c);
  s = new boolshare(one_gate->get_wires(), c);

  for (int j = l; j <= nm; j++) {
    int i = j + 1;
    g = new boolshare(rv1[i]->get_wires(), c);
    y = new boolshare(w[i]->get_wires(), c);
    h = c->PutFPGate(s, g, MUL, bitlen, 1, no_status);
    g = c->PutFPGate(cc, g, MUL, bitlen, 1, no_status);
    *z = BuildPythagCircuit(f, h, c, bitlen);
    rv1[j] = new boolshare((*z)->get_wires(), c);
    cc = c->PutFPGate(f, *z, DIV, bitlen, 1, no_status);
    s = c->PutFPGate(h, *z, DIV, bitlen, 1, no_status);
    {  // f = x * c + g * s;
      share* term1 = c->PutFPGate(x, cc, MUL, bitlen, 1, no_status);
      share* term2 = c->PutFPGate(g, s, MUL, bitlen, 1, no_status);
      f = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
    }
    {  // g = g * c - x * s;
      share* term1 = c->PutFPGate(g, cc, MUL, bitlen, 1, no_status);
      share* term2 = c->PutFPGate(x, s, MUL, bitlen, 1, no_status);
      g = c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
    }
    h = c->PutFPGate(y, s, MUL, bitlen, 1, no_status);
    y = c->PutFPGate(y, cc, MUL, bitlen, 1, no_status);
    for (jj = 0; jj < nCols; jj++) {
      x = new boolshare(v[jj][j]->get_wires(), c);
      *z = new boolshare(v[jj][i]->get_wires(), c);
      {  // v[jj][j] = x * c + z * s;
        share* term1 = c->PutFPGate(x, cc, MUL, bitlen, 1, no_status);
        share* term2 = c->PutFPGate(*z, s, MUL, bitlen, 1, no_status);
        v[jj][j] = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
      }
      {  // v[jj][i] = z * c - x * s;
        share* term1 = c->PutFPGate(*z, cc, MUL, bitlen, 1, no_status);
        share* term2 = c->PutFPGate(x, s, MUL, bitlen, 1, no_status);
        v[jj][i] = c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
      }
    }
    *z = BuildPythagCircuit(f, h, c, bitlen);
    w[j] = new boolshare((*z)->get_wires(), c);
    {  // if (z) {
      share* temp;
      share* ifz = (*z)->get_wire_ids_as_share(0);
      for (uint32_t tt = 1; tt < bitlen - 1;
           tt++) {  // -1 -> do not include sign bit
        share* tempbit = (*z)->get_wire_ids_as_share(tt);
        temp = ifz;
        ifz = c->PutORGate(temp, tempbit);
        delete tempbit;
        delete temp;
      }

      temp = c->PutFPGate(one_gate, *z, DIV, bitlen, 1, no_status);
      *z = c->PutMUXGate(temp, *z, ifz);
      delete temp;

      temp = c->PutFPGate(f, *z, MUL, bitlen, 1, no_status);
      cc = c->PutMUXGate(temp, cc, ifz);
      delete temp;

      temp = c->PutFPGate(h, *z, MUL, bitlen, 1, no_status);
      s = c->PutMUXGate(temp, s, ifz);
      delete temp;
    }
    {  // f = c * g + s * y;
      share* term1 = c->PutFPGate(cc, g, MUL, bitlen, 1, no_status);
      share* term2 = c->PutFPGate(s, y, MUL, bitlen, 1, no_status);
      f = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
    }
    {  // x = c * y - s * g;
      share* term1 = c->PutFPGate(cc, y, MUL, bitlen, 1, no_status);
      share* term2 = c->PutFPGate(s, g, MUL, bitlen, 1, no_status);
      x = c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
    }
    for (jj = 0; jj < nRows; jj++) {
      y = new boolshare(a[jj][j]->get_wires(), c);
      *z = new boolshare(a[jj][i]->get_wires(), c);
      {  // a[jj][j] = y * c + z * s;
        share* term1 = c->PutFPGate(y, cc, MUL, bitlen, 1, no_status);
        share* term2 = c->PutFPGate(*z, s, MUL, bitlen, 1, no_status);
        a[jj][j] = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
      }
      {  // a[jj][i] = z * c - y * s;
        share* term1 = c->PutFPGate(*z, cc, MUL, bitlen, 1, no_status);
        share* term2 = c->PutFPGate(y, s, MUL, bitlen, 1, no_status);
        a[jj][i] = c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
      }
    }
  }
  rv1[l] = new boolshare(zero_gate->get_wires(), c);
  rv1[k] = new boolshare(f->get_wires(), c);
  w[k] = new boolshare(x->get_wires(), c);
}
void RunSvdPart2Circuit(uint32_t** a, int nRows, int nCols, uint32_t* w,
                        uint32_t** v, uint32_t* rv1, uint32_t* z, int k, int l,
                        int nm, BooleanCircuit* c, ABYParty* p) {

  std::vector<Sharing*>& sharings = p->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  // Allocate space for shares
  share* temp;
  share*** s_a = new share**[nRows];
  for (int t = 0; t < nRows; t++) {
    s_a[t] = new share*[nCols];
  }
  share** s_w = new share*[nCols];
  share*** s_v = new share**[nCols];
  for (int t = 0; t < nCols; t++) {
    s_v[t] = new share*[nCols];
  }
  share** s_rv1 = new share*[nCols];
  share* s_z;

  // Prepare inputs
  for (int t = 0; t < nRows; t++) {
    for (int tt = 0; tt < nCols; tt++) {
      s_a[t][tt] = bc->PutSharedINGate(&a[t][tt], bitlen);
    }
  }
  for (int t = 0; t < nCols; t++) {
    s_w[t] = bc->PutSharedINGate(&w[t], bitlen);
  }
  for (int t = 0; t < nCols; t++) {
    for (int tt = 0; tt < nCols; tt++) {
      s_v[t][tt] = bc->PutSharedINGate(&v[t][tt], bitlen);
    }
  }
  for (int t = 0; t < nCols; t++) {
    s_rv1[t] = bc->PutSharedINGate(&rv1[t], bitlen);
  }
  s_z = bc->PutSharedINGate(z, bitlen);

  // if yao, convert bool inputs to yao
  if (c->GetContext() == S_YAO) {
    for (int t = 0; t < nRows; t++) {
      for (int tt = 0; tt < nCols; tt++) {
        share* temp = s_a[t][tt];
        s_a[t][tt] = c->PutB2YGate(temp);
        delete temp;
      }
    }
    for (int t = 0; t < nCols; t++) {
      share* temp = s_w[t];
      s_w[t] = c->PutB2YGate(temp);
      delete temp;
    }
    for (int t = 0; t < nCols; t++) {
      for (int tt = 0; tt < nCols; tt++) {
        share* temp = s_v[t][tt];
        s_v[t][tt] = c->PutB2YGate(temp);
        delete temp;
      }
    }
    for (int t = 0; t < nCols; t++) {
      temp = s_rv1[t];
      s_rv1[t] = c->PutB2YGate(temp);
      delete temp;
    }
    temp = s_z;
    s_z = c->PutB2YGate(temp);
    delete temp;
  }

  // prepare circuit
  BuildSvdPart2Circuit(s_a, nRows, nCols, s_w, s_v, s_rv1, &s_z, k, l, nm, c);

  // if yao was used, convert from yao back to boolean circuit
  if (c->GetContext() == S_YAO) {
    for (int t = 0; t < nRows; t++) {
      for (int tt = 0; tt < nCols; tt++) {
        share* temp = s_a[t][tt];
        s_a[t][tt] = bc->PutY2BGate(temp);
        delete temp;
      }
    }
    for (int t = 0; t < nCols; t++) {
      share* temp = s_w[t];
      s_w[t] = bc->PutY2BGate(temp);
      delete temp;
    }
    for (int t = 0; t < nCols; t++) {
      for (int tt = 0; tt < nCols; tt++) {
        share* temp = s_v[t][tt];
        s_v[t][tt] = bc->PutY2BGate(temp);
        delete temp;
      }
    }
    for (int t = 0; t < nCols; t++) {
      temp = s_rv1[t];
      s_rv1[t] = bc->PutY2BGate(temp);
      delete temp;
    }
    temp = s_z;
    s_z = bc->PutY2BGate(temp);
    delete temp;
  }
  // prepare output (everything except z which is not an output)
  for (int t = 0; t < nRows; t++) {
    for (int tt = 0; tt < nCols; tt++) {
      share* temp = s_a[t][tt];
      s_a[t][tt] = bc->PutSharedOUTGate(temp);
      delete temp;
    }
  }
  for (int t = 0; t < nCols; t++) {
    share* temp = s_w[t];
    s_w[t] = bc->PutSharedOUTGate(temp);
    delete temp;
  }
  for (int t = 0; t < nCols; t++) {
    for (int tt = 0; tt < nCols; tt++) {
      share* temp = s_v[t][tt];
      s_v[t][tt] = bc->PutSharedOUTGate(temp);
      delete temp;
    }
  }
  for (int t = 0; t < nCols; t++) {
    temp = s_rv1[t];
    s_rv1[t] = bc->PutSharedOUTGate(temp);
    delete temp;
  }
  temp = s_z;
  s_z = bc->PutSharedOUTGate(temp);
  delete temp;

  p->ExecCircuit();

  // get shared output
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      a[i][j] = s_a[i][j]->get_clear_value<uint32_t>();
      delete s_a[i][j];
    }
    delete[] s_a[i];
  }
  delete[] s_a;
  for (int i = 0; i < nCols; i++) {
    w[i] = s_w[i]->get_clear_value<uint32_t>();
    delete s_w[i];
  }
  delete[] s_w;
  for (int i = 0; i < nCols; i++) {
    for (int j = 0; j < nCols; j++) {
      v[i][j] = s_v[i][j]->get_clear_value<uint32_t>();
      delete s_v[i][j];
    }
    delete[] s_v[i];
  }
  delete[] s_v;
  for (int i = 0; i < nCols; i++) {
    rv1[i] = s_rv1[i]->get_clear_value<uint32_t>();
    delete s_rv1[i];
  }
  delete[] s_rv1;
  *z = s_z->get_clear_value<uint32_t>();
  delete s_z;

  collectTiming();
  collectCommunication();
  p->Reset();
}

// uint32_t arguments must be from PutSharedOUTGate()
// then calling get_clear_value() and circuit->Reset();
// This function builds/executes multiple circuits.
void RunSvdLoopLeak(uint32_t** a, int nRows, int nCols, uint32_t* w,
                    uint32_t** v, BooleanCircuit* c, ABYParty* p, e_role role) {

  // Following values are shared between subcircuits
  int flag, i, its, /*j, jj,*/ k, l, nm;
  uint32_t anorm, cc, f, /*g, h,*/ s, /*scale, x, y,*/ z;
  uint32_t* rv1 = new uint32_t[nCols];

  RunSvdPart1Circuit(a, nRows, nCols, w, v, &anorm, rv1, c, p);

  for (k = nCols - 1; k >= 0; k--) {
    for (its = 0; its < 30; its++) {
      cout << k;
      flag = 1;
      for (l = k; l >= 0; l--) {
        nm = l - 1;
        if (RunSvdCheckZeroCircuit(&rv1[l], &anorm, c, p)) {
          flag = 0;
          break;
        }
        if (RunSvdCheckZeroCircuit(&w[nm], &anorm, c, p)) {
          break;
        }
      }

      if (flag) {
        // setting both shares to zero = constant 0
        cc = 0;
        // setting party0 to 0 and party1 to 1 = constant 1
        float frole = (float)role;
        uint32_t* frolep = (uint32_t*)&frole;
        s = *frolep;

        for (i = l; i <= k; i++) {
          RunSvdFlag1Circuit(&cc, &f, &s, &rv1[i], c, p);

          if (RunSvdCheckZeroCircuit(&f, &anorm, c, p)) {
            break;
          }

          // Slice out column i of a (a[*][i]) and column nm (a[*][nn])
          uint32_t* aSTARi = new uint32_t[nRows];
          uint32_t* aSTARnm = new uint32_t[nRows];
          for (int star = 0; star < nRows; star++) {
            aSTARi[star] = a[star][i];
            aSTARnm[star] = a[star][nm];
          }
          RunSvdFlag2Circuit(aSTARi, aSTARnm, nRows, &w[i], &cc, &f, &s, c, p);
          // Copy column back into matrix
          for (int star = 0; star < nRows; star++) {
            a[star][i] = aSTARi[star];
            a[star][nm] = aSTARnm[star];
          }
          delete[] aSTARi;
          delete[] aSTARnm;
        }
      }
      z = w[k];
      if (l == k) {
        // Slice out column k of v (v[*][k])
        uint32_t* vSTARk = new uint32_t[nCols];
        for (int star = 0; star < nCols; star++) {
          vSTARk[star] = v[star][k];
        }
        RunSvdZCircuit(nCols, &z, &w[k], vSTARk, c, p);
        // Copy column back into matrix
        for (int star = 0; star < nCols; star++) {
          v[star][k] = vSTARk[star];
        }
        delete[] vSTARk;
        break;
      }
      if (its == 29)
        printf("no convergence in 30 svdcmp iterations\n");
      nm = k - 1;

      RunSvdPart2Circuit(a, nRows, nCols, w, v, rv1, &z, k, l, nm, c, p);
    }
    printf("/");
    fflush(stdout);
  }
}

// Wrapper function around RunSvdLoopLeak which
// takes shares instead of secret shares.
// Creates dummy circuits to build secret shares.
void BuildAndRunSvdLoopLeak(share*** s_a, int nRows, int nCols, share** s_w,
                            share*** s_v, BooleanCircuit* c, ABYParty* party,
                            e_role role, std::function<void()> toRawShares,
                            std::function<void()> toShareObjects) {

  std::vector<Sharing*>& sharings = party->GetSharings();
  BooleanCircuit* bc =
      (BooleanCircuit*)sharings[S_BOOL]->GetCircuitBuildRoutine();

  share* temp;

  // stores secret-shared output
  uint32_t** a = new uint32_t*[nRows];
  for (int i = 0; i < nRows; i++) {
    a[i] = new uint32_t[nCols];
  }
  uint32_t* w = new uint32_t[nCols];
  uint32_t** v = new uint32_t*[nCols];
  for (int i = 0; i < nCols; i++) {
    v[i] = new uint32_t[nCols];
  }

  // if yao is prefered circuit and yao shares passed in, must
  // convert to bool to get raw shares.
  // Tests must pass in bool shares even if yao is preferred.
  if (c->GetContext() == S_YAO && s_a[0][0]->get_share_type() == S_YAO) {
    for (int i = 0; i < nRows; i++) {
      for (int j = 0; j < nCols; j++) {
        temp = s_a[i][j];
        s_a[i][j] = bc->PutY2BGate(temp);
        delete temp;
      }
    }
  }
  // build shared output objects
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      temp = s_a[i][j];
      s_a[i][j] = bc->PutSharedOUTGate(temp);
      delete temp;
    }
  }

  // run the dummy circuit
  CLOCK(ShareSvdInputs);
  TIC(ShareSvdInputs);
  party->ExecCircuit();
  TOC(ShareSvdInputs);
  // Get raw output shares
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      a[i][j] = s_a[i][j]->get_clear_value<uint32_t>();
      delete s_a[i][j];
    }
  }

  toRawShares();  // allow caller to store shares across sub-circuit executions

  collectTiming();
  collectCommunication();
  party->Reset();
  // a[][] now contains secret shared values of plaintext from each party

  // Run the svd on the secret shared data
  CLOCK(SVD);
  TIC(SVD);
  RunSvdLoopLeak(a, nRows, nCols, w, v, (BooleanCircuit*)c, party, role);
  TOC(SVD);

  //// next run another dummy circuit to get the cleartext output from shared
  /// output
  // shared inputs to circuit
  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nCols; j++) {
      s_a[i][j] = bc->PutSharedINGate(&a[i][j], bitlen);
    }
  }
  for (int i = 0; i < nCols; i++) {
    s_w[i] = bc->PutSharedINGate(&w[i], bitlen);
  }
  for (int i = 0; i < nCols; i++) {
    for (int j = 0; j < nCols; j++) {
      s_v[i][j] = bc->PutSharedINGate(&v[i][j], bitlen);
    }
  }
  // if yao, convert bool back to yao for future executions
  if (c->GetContext() == S_YAO) {
    for (int i = 0; i < nRows; i++) {
      for (int j = 0; j < nCols; j++) {
        temp = s_a[i][j];
        s_a[i][j] = c->PutB2YGate(temp);
        delete temp;
      }
    }
    for (int i = 0; i < nCols; i++) {
      temp = s_w[i];
      s_w[i] = c->PutB2YGate(temp);
      delete temp;
    }
    for (int i = 0; i < nCols; i++) {
      for (int j = 0; j < nCols; j++) {
        temp = s_v[i][j];
        s_v[i][j] = c->PutB2YGate(temp);
        delete temp;
      }
    }
  }

  toShareObjects();  // allow caller to recover shares across sub-circuit
                     // executions
}

// Build SVD circuit. Since it does not need
// intermediate values, it does not use subcircuits
// thus only builds the circuit but does not execute.
void BuildSvdDO(share*** a, int nRows, int nCols, share** w, share*** v,
                BooleanCircuit* c, ABYParty* party, e_role role) {

  // Following values are shared between subcircuits
  int /*flag,*/ i, its, j, jj, k, l, nm;
  share *anorm = nullptr, *cc = nullptr, *f = nullptr, *g = nullptr,
        *h = nullptr, *s = nullptr, *scale = nullptr, *x = nullptr,
        *y = nullptr, *z = nullptr;
  share** rv1 = new share*[nCols];

  float zero = 0;
  share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);
  int int_zero = 0;
  share* int_zero_gate = c->PutCONSGate((uint32_t*)&int_zero, bitlen);
  float one = 1;
  share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);
  int int_one = 1;
  share* int_one_gate = c->PutCONSGate((uint32_t*)&int_one, bitlen);
  float two = 2;
  share* two_gate = c->PutCONSGate((uint32_t*)&two, bitlen);
  uint8_t true_val = 1;
  share* true_gate = c->PutCONSGate((uint8_t*)&true_val, 1);
  uint8_t false_val = 0;
  share* false_gate = c->PutCONSGate((uint8_t*)&false_val, 1);
  share* temp;

  x = new boolshare(zero_gate->get_wires(), c);  // x must be initialized

  BuildSvdPart1Circuit(a, nRows, nCols, w, v, &anorm, rv1, c);

  for (k = nCols - 1; k >= 0; k--) {
    share* converged = new boolshare(false_gate->get_wires(), c);
    for (its = 0; its < 30; its++) {
      cout << k << flush;
      share* flag = new boolshare(true_gate->get_wires(), c);
      share* l_found = new boolshare(false_gate->get_wires(), c);
      share* priv_l = new boolshare(int_zero_gate->get_wires(), c);

      for (l = k; l >= 0; l--) {
        nm = l - 1;

        share* this_l = c->PutCONSGate((uint32_t*)&l, bitlen);

        share* rv1lcopy = new boolshare(rv1[l]->get_wires(),
                                        c);  // check zero overwrites with Fabs
        share* rv1isanorm = BuildSvdCheckZeroCircuit(rv1lcopy, anorm, c);
        delete rv1lcopy;
        // flag = flag & !(rv1isanorm & !l_found);
        share* nl_found = c->PutINVGate(l_found);
        share* r_a_nl_found = c->PutANDGate(rv1isanorm, nl_found);
        share* nr_a_nl_found = c->PutINVGate(r_a_nl_found);
        temp = flag;
        flag = c->PutANDGate(flag, nr_a_nl_found);
        delete temp;
        temp = priv_l;
        priv_l = c->PutMUXGate(this_l, priv_l, r_a_nl_found);
        delete temp;
        temp = l_found;
        l_found = c->PutORGate(l_found, rv1isanorm);
        delete temp;
        delete nr_a_nl_found;
        delete r_a_nl_found;
        delete nl_found;
        delete rv1isanorm;

        if (nm >= 0) {
          share* wnmcopy = new boolshare(w[nm]->get_wires(),
                                         c);  // check zero overwrites with Fabs
          share* wnmisanorm = BuildSvdCheckZeroCircuit(wnmcopy, anorm, c);
          delete wnmcopy;
          share* nl_found = c->PutINVGate(l_found);
          share* wnm_a_nl_found = c->PutANDGate(wnmisanorm, nl_found);
          temp = priv_l;
          priv_l = c->PutMUXGate(this_l, priv_l, wnm_a_nl_found);
          delete temp;
          temp = l_found;
          l_found = c->PutORGate(l_found, wnmisanorm);
          delete temp;
          delete wnm_a_nl_found;
          delete nl_found;
          delete wnmisanorm;
        }
      }

      share* priv_nm =
          c->PutSUBGate(priv_l, int_one_gate);  // integer subtraction
      {                                         // if (flag)
        cc = new boolshare(zero_gate->get_wires(), c);
        s = new boolshare(one_gate->get_wires(), c);
        for (i = 1; i <= k; i++) {
          share* priv_i = c->PutCONSGate((uint32_t*)&i, bitlen);
          temp = c->PutGTGate(priv_l, priv_i);
          share* igel = c->PutINVGate(temp);
          delete temp;
          delete priv_i;
          temp = c->PutANDGate(flag, igel);
          share* nconverged = c->PutINVGate(converged);
          share* modify = c->PutANDGate(temp, nconverged);
          delete temp;
          delete nconverged;

          if (f)
            delete f;
          f = c->PutFPGate(s, rv1[i], MUL, bitlen, 1, no_status);
          temp = rv1[i];
          share* tmp2 = c->PutFPGate(cc, rv1[i], MUL, bitlen, 1, no_status);
          rv1[i] = c->PutMUXGate(tmp2, rv1[i], modify);
          delete temp;
          delete tmp2;

          share* fcopy = new boolshare(f->get_wires(),
                                       c);  // check zero overwrites with Fabs
          share* fisanorm = BuildSvdCheckZeroCircuit(fcopy, anorm, c);
          delete fcopy;
          share* nfisanorm = c->PutINVGate(fisanorm);
          temp = modify;
          modify = c->PutANDGate(modify, nfisanorm);
          delete temp;
          delete nfisanorm;
          delete fisanorm;

          if (g)
            delete g;
          g = new boolshare(w[i]->get_wires(), c);
          if (h)
            delete h;
          h = BuildPythagCircuit(f, g, c, bitlen);
          temp = w[i];
          w[i] = c->PutMUXGate(h, w[i], modify);
          delete temp;
          share* oldh = h;
          h = c->PutFPGate(one_gate, h, DIV, bitlen, 1, no_status);
          delete oldh;
          temp = c->PutFPGate(g, h, MUL, bitlen, 1, no_status);
          share* oldcc = cc;
          cc = c->PutMUXGate(temp, cc, modify);
          delete temp;
          delete oldcc;
          temp = c->PutFPGate(f, h, MUL, bitlen, 1, no_status);
          BuildNegativeCircuit(temp, c);
          share* olds = s;
          s = c->PutMUXGate(temp, s, modify);
          delete temp;
          delete olds;
          for (int j = 0; j < nRows; j++) {
            for (int nm_finder = 0; nm_finder < k; nm_finder++) {
              share* priv_nm_finder =
                  c->PutCONSGate((uint32_t*)&nm_finder, bitlen);
              share* found_nm = c->PutEQGate(priv_nm, priv_nm_finder);
              delete priv_nm_finder;
              share* mafnm = c->PutANDGate(modify, found_nm);
              delete found_nm;

              y = new boolshare(a[j][nm_finder]->get_wires(), c);
              z = new boolshare(a[j][i]->get_wires(), c);
              share* tmp1 = c->PutFPGate(y, cc, MUL, bitlen, 1, no_status);
              share* tmp2 = c->PutFPGate(z, s, MUL, bitlen, 1, no_status);
              share* tmp3 = c->PutFPGate(tmp1, tmp2, ADD, bitlen, 1, no_status);
              temp = a[j][nm_finder];
              a[j][nm_finder] = c->PutMUXGate(tmp3, a[j][nm_finder], mafnm);
              delete temp;
              delete tmp1;
              delete tmp2;
              delete tmp3;

              tmp1 = c->PutFPGate(z, cc, MUL, bitlen, 1, no_status);
              tmp2 = c->PutFPGate(y, s, MUL, bitlen, 1,
                                  no_status);  // y is broken?
              tmp3 = c->PutFPGate(tmp1, tmp2, SUB, bitlen, 1, no_status);
              temp = a[j][i];
              a[j][i] = c->PutMUXGate(tmp3, a[j][i], mafnm);
              delete temp;
              delete tmp1;
              delete tmp2;
              delete tmp3;

              delete y;
              delete z;
              delete mafnm;
            }
          }
        }
      }
      z = w[k];
      share* tempk = c->PutCONSGate((uint32_t*)&k, 32);
      share* newconverged = c->PutEQGate(priv_l, tempk);
      delete tempk;
      {    // if (l == k)
        {  // if (z < 0.0) singular value is made nonnegative
          share* isneg = c->PutFPGate(zero_gate, z, CMP);
          share* modify = c->PutANDGate(isneg, newconverged);
          share* nconverged = c->PutINVGate(converged);
          share* oldmodify = modify;
          modify = c->PutANDGate(modify, nconverged);
          delete oldmodify;

          share* temp = new boolshare(z->get_wires(), c);
          BuildNegativeCircuit(temp, c);
          w[k] = c->PutMUXGate(temp, w[k], modify);
          for (int j = 0; j < nCols; j++) {
            share* negv = new boolshare(v[j][k]->get_wires(), c);
            BuildNegativeCircuit(negv, c);
            share* oldvjk = v[j][k];
            v[j][k] = c->PutMUXGate(negv, v[j][k], modify);
            delete negv;
            delete oldvjk;
          }
        }
      }
      share* oldconverged = converged;
      converged = c->PutORGate(converged, newconverged);
      delete oldconverged;
      if (k == 0)
        break;  // if k == 0 we are done, code below segfaults from nm

      // shift from bottom 2-by-2 minor
      for (int l_finder = 0; l_finder < k; l_finder++) {
        share* pl_finder = c->PutCONSGate((uint32_t*)&l_finder, bitlen);
        share* found_l = c->PutEQGate(priv_l, pl_finder);
        delete pl_finder;
        share* flanc = c->PutINVGate(converged);
        share* oldflanc = flanc;
        flanc = c->PutANDGate(found_l, flanc);
        delete oldflanc;
        delete found_l;
        share* oldx = x;
        x = c->PutMUXGate(w[l_finder], x, flanc);
        delete oldx;
        delete flanc;
      }

      nm = k - 1;
      y = new boolshare(w[nm]->get_wires(),
                        c);  // do not delete before assigning
      if (g)
        delete g;
      g = new boolshare(rv1[nm]->get_wires(), c);
      if (h)
        delete h;
      h = new boolshare(rv1[k]->get_wires(), c);
      {  // f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        share* ymz = c->PutFPGate(y, z, SUB, bitlen, 1, no_status);
        share* ypz = c->PutFPGate(y, z, ADD, bitlen, 1, no_status);
        share* temp1 = c->PutFPGate(ymz, ypz, MUL, bitlen, 1, no_status);
        share* gmh = c->PutFPGate(g, h, SUB, bitlen, 1, no_status);
        share* gph = c->PutFPGate(g, h, ADD, bitlen, 1, no_status);
        share* temp2 = c->PutFPGate(gmh, gph, MUL, bitlen, 1, no_status);
        share* num = c->PutFPGate(temp1, temp2, ADD, bitlen, 1, no_status);
        share* den1 = c->PutFPGate(two_gate, h, MUL, bitlen, 1, no_status);
        share* den2 = c->PutFPGate(den1, y, MUL, bitlen, 1, no_status);
        f = c->PutFPGate(num, den2, DIV, bitlen, 1, no_status);
        delete ymz;
        delete ypz;
        delete temp1;
        delete gmh;
        delete gph;
        delete temp2;
        delete num;
        delete den1;
        delete den2;
      }
      g = BuildPythagCircuit(f, one_gate, c, bitlen);
      {  // f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f)))- h)) / x;
        share* xmz = c->PutFPGate(x, z, SUB, bitlen, 1, no_status);
        share* xpz = c->PutFPGate(x, z, ADD, bitlen, 1, no_status);
        share* term1 = c->PutFPGate(xmz, xpz, MUL, bitlen, 1, no_status);
        share* sgf = new boolshare(g->get_wires(), c);
        BuildSignCircuit(sgf, f, c);
        share* fpsgf = c->PutFPGate(f, sgf, ADD, bitlen, 1, no_status);
        share* ydfpsgf = c->PutFPGate(y, fpsgf, DIV, bitlen, 1, no_status);
        share* uggg = c->PutFPGate(ydfpsgf, h, SUB, bitlen, 1, no_status);
        share* term2 = c->PutFPGate(h, uggg, MUL, bitlen, 1, no_status);
        share* num2 = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
        delete f;
        f = c->PutFPGate(num2, x, DIV, bitlen, 1, no_status);
        delete xmz;
        delete xpz;
        delete term1;
        delete sgf;
        delete fpsgf;
        delete ydfpsgf;
        delete uggg;
        delete term2;
        delete num2;
      }
      cc = new boolshare(one_gate->get_wires(), c);
      s = new boolshare(one_gate->get_wires(), c);

      for (int j = 0; j <= nm; j++) {
        share* priv_j = c->PutCONSGate((uint32_t*)&j, bitlen);
        share* jgtl = c->PutGTGate(priv_l, priv_j);  // have to invert for >=
        share* oldjgtl = jgtl;
        jgtl = c->PutINVGate(jgtl);
        delete oldjgtl;
        share* nconverged = c->PutINVGate(converged);
        share* modify = c->PutANDGate(jgtl, nconverged);
        delete nconverged;
        delete jgtl;
        delete priv_j;

        int i = j + 1;
        delete g;
        g = new boolshare(rv1[i]->get_wires(), c);
        delete y;
        y = new boolshare(w[i]->get_wires(), c);
        delete h;
        h = c->PutFPGate(s, g, MUL, bitlen, 1, no_status);
        share* oldg = g;
        g = c->PutFPGate(cc, g, MUL, bitlen, 1, no_status);
        delete oldg;
        delete z;
        z = BuildPythagCircuit(f, h, c, bitlen);
        share* oldr = rv1[j];
        rv1[j] = c->PutMUXGate(z, rv1[j], modify);
        delete oldr;
        temp = c->PutFPGate(f, z, DIV, bitlen, 1, no_status);
        share* oldcc = cc;
        cc = c->PutMUXGate(temp, cc, modify);
        delete oldcc;
        delete temp;
        share* olds = s;
        temp = c->PutFPGate(h, z, DIV, bitlen, 1, no_status);
        s = c->PutMUXGate(temp, s, modify);
        delete olds;
        delete temp;
        {  // f = x * c + g * s;
          share* term1 = c->PutFPGate(x, cc, MUL, bitlen, 1, no_status);
          share* term2 = c->PutFPGate(g, s, MUL, bitlen, 1, no_status);
          share* added = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
          share* oldf = f;
          f = c->PutMUXGate(added, f, modify);
          delete oldf;
          delete added;
          delete term1;
          delete term2;
        }
        {  // g = g * c - x * s;
          share* term1 = c->PutFPGate(g, cc, MUL, bitlen, 1, no_status);
          share* term2 = c->PutFPGate(x, s, MUL, bitlen, 1, no_status);
          share* oldg = g;
          g = c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
          delete oldg;
          delete term1;
          delete term2;
        }
        h = c->PutFPGate(y, s, MUL, bitlen, 1, no_status);
        y = c->PutFPGate(y, cc, MUL, bitlen, 1, no_status);
        for (jj = 0; jj < nCols; jj++) {
          share* oldx = x;
          x = c->PutMUXGate(v[jj][j], x, modify);
          delete oldx;
          z = new boolshare(v[jj][i]->get_wires(), c);
          {  // v[jj][j] = x * c + z * s;
            share* term1 = c->PutFPGate(x, cc, MUL, bitlen, 1, no_status);
            share* term2 = c->PutFPGate(z, s, MUL, bitlen, 1, no_status);
            share* added =
                c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
            share* oldv = v[jj][j];
            v[jj][j] = c->PutMUXGate(added, v[jj][j], modify);
            delete oldv;
            delete added;
            delete term1;
            delete term2;
          }
          {  // v[jj][i] = z * c - x * s;
            share* term1 = c->PutFPGate(z, cc, MUL, bitlen, 1, no_status);
            share* term2 = c->PutFPGate(x, s, MUL, bitlen, 1, no_status);
            share* subbed =
                c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
            share* oldv = v[jj][i];
            v[jj][i] = c->PutMUXGate(subbed, v[jj][i], modify);
            delete oldv;
            delete subbed;
            delete term1;
            delete term2;
          }
        }
        z = BuildPythagCircuit(f, h, c, bitlen);
        share* oldw = w[j];
        w[j] = c->PutMUXGate(z, w[j], modify);
        delete oldw;
        {  // if (z) {
          share* temp;
          share* ifz = z->get_wire_ids_as_share(0);
          for (uint32_t tt = 1; tt < bitlen - 1;
               tt++) {  // -1 -> do not include sign bit
            share* tempbit = z->get_wire_ids_as_share(tt);
            temp = ifz;
            ifz = c->PutORGate(temp, tempbit);
            delete tempbit;
            delete temp;
          }
          share* ifzam = c->PutANDGate(ifz, modify);
          delete ifz;

          temp = c->PutFPGate(one_gate, z, DIV, bitlen, 1, no_status);
          z = c->PutMUXGate(temp, z, ifzam);
          delete temp;

          temp = c->PutFPGate(f, z, MUL, bitlen, 1, no_status);
          cc = c->PutMUXGate(temp, cc, ifzam);
          delete temp;

          temp = c->PutFPGate(h, z, MUL, bitlen, 1, no_status);
          s = c->PutMUXGate(temp, s, ifzam);
          delete temp;
        }
        {  // f = c * g + s * y;
          share* term1 = c->PutFPGate(cc, g, MUL, bitlen, 1, no_status);
          share* term2 = c->PutFPGate(s, y, MUL, bitlen, 1, no_status);
          share* added = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
          share* oldf = f;
          f = c->PutMUXGate(added, f, modify);
          delete oldf;
          delete added;
          delete term1;
          delete term2;
        }
        {  // x = c * y - s * g;
          share* term1 = c->PutFPGate(cc, y, MUL, bitlen, 1, no_status);
          share* term2 = c->PutFPGate(s, g, MUL, bitlen, 1, no_status);
          share* subbed = c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
          share* oldx = x;
          x = c->PutMUXGate(subbed, x, modify);
          delete oldx;
          delete subbed;
          delete term1;
          delete term2;
        }
        for (jj = 0; jj < nRows; jj++) {
          y = new boolshare(a[jj][j]->get_wires(), c);
          z = new boolshare(a[jj][i]->get_wires(), c);
          {  // a[jj][j] = y * c + z * s;
            share* term1 = c->PutFPGate(y, cc, MUL, bitlen, 1, no_status);
            share* term2 = c->PutFPGate(z, s, MUL, bitlen, 1, no_status);
            share* added =
                c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
            share* olda = a[jj][j];
            a[jj][j] = c->PutMUXGate(added, a[jj][j], modify);
            delete olda;
            delete added;
            delete term1;
            delete term2;
          }
          {  // a[jj][i] = z * c - y * s;
            share* term1 = c->PutFPGate(z, cc, MUL, bitlen, 1, no_status);
            share* term2 = c->PutFPGate(y, s, MUL, bitlen, 1, no_status);
            share* subbed =
                c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
            share* olda = a[jj][i];
            a[jj][i] = c->PutMUXGate(subbed, a[jj][i], modify);
            delete olda;
            delete subbed;
            delete term1;
            delete term2;
          }
        }
      }
      for (int l_finder = 0; l_finder < k; l_finder++) {
        share* pl_finder = c->PutCONSGate((uint32_t*)&l_finder, 32);
        share* found_l = c->PutEQGate(priv_l, pl_finder);
        delete pl_finder;
        share* flanc = c->PutINVGate(converged);
        share* oldflanc = flanc;
        flanc = c->PutANDGate(found_l, flanc);
        delete oldflanc;
        delete found_l;
        share* oldr = rv1[l_finder];
        rv1[l_finder] = c->PutMUXGate(zero_gate, rv1[l_finder], flanc);
        delete oldr;
        delete flanc;
      }
      share* nconverged = c->PutINVGate(converged);
      share* oldr = rv1[k];
      rv1[k] = c->PutMUXGate(f, rv1[k], nconverged);
      delete oldr;
      share* oldw = w[k];
      w[k] = c->PutMUXGate(x, w[k], nconverged);
      delete oldw;
      delete nconverged;
    }
    delete converged;
    printf("/");
    fflush(stdout);
  }
}

// Build SVD circuit. Since it does not need
// intermediate values, it does not use subcircuits
// thus only builds the circuit but does not execute.
void BuildSvdOSL(share*** a, int nRows, int nCols, share** w, share*** v,
                 BooleanCircuit* c, ABYParty* party, e_role role) {

  // Following values are shared between subcircuits
  int /*flag,*/ i, its, j, jj, k, l, nm;
  share *anorm = nullptr, *cc = nullptr, *f = nullptr, *g = nullptr,
        *h = nullptr, *s = nullptr, *scale = nullptr, *x = nullptr,
        *y = nullptr, *z = nullptr;
  share** rv1 = new share*[nCols];

  float zero = 0;
  share* zero_gate = c->PutCONSGate((uint32_t*)&zero, bitlen);
  int int_zero = 0;
  share* int_zero_gate = c->PutCONSGate((uint32_t*)&int_zero, bitlen);
  float one = 1;
  share* one_gate = c->PutCONSGate((uint32_t*)&one, bitlen);
  int int_one = 1;
  share* int_one_gate = c->PutCONSGate((uint32_t*)&int_one, bitlen);
  float two = 2;
  share* two_gate = c->PutCONSGate((uint32_t*)&two, bitlen);
  uint8_t true_val = 1;
  share* true_gate = c->PutCONSGate((uint8_t*)&true_val, 1);
  uint8_t false_val = 0;
  share* false_gate = c->PutCONSGate((uint8_t*)&false_val, 1);
  share* temp;

  x = new boolshare(zero_gate->get_wires(), c);  // x must be initialized

  BuildSvdPart1Circuit(a, nRows, nCols, w, v, &anorm, rv1, c);

  for (k = nCols - 1; k >= 0; k--) {
    share* converged = new boolshare(false_gate->get_wires(), c);
    for (its = 0; its < 2; its++) {
      cout << k << flush;
      share* priv_l = new boolshare(int_zero_gate->get_wires(), c);

      z = w[k];
      share* tempk = c->PutCONSGate((uint32_t*)&k, 32);
      share* newconverged = c->PutEQGate(priv_l, tempk);
      delete tempk;
      {    // if (l == k)
        {  // if (z < 0.0) singular value is made nonnegative
          share* isneg = c->PutFPGate(zero_gate, z, CMP);
          share* modify = c->PutANDGate(isneg, newconverged);
          share* nconverged = c->PutINVGate(converged);
          share* oldmodify = modify;
          modify = c->PutANDGate(modify, nconverged);
          delete oldmodify;

          share* temp = new boolshare(z->get_wires(), c);
          BuildNegativeCircuit(temp, c);
          w[k] = c->PutMUXGate(temp, w[k], modify);
          for (int j = 0; j < nCols; j++) {
            share* negv = new boolshare(v[j][k]->get_wires(), c);
            BuildNegativeCircuit(negv, c);
            share* oldvjk = v[j][k];
            v[j][k] = c->PutMUXGate(negv, v[j][k], modify);
            delete negv;
            delete oldvjk;
          }
        }
      }
      share* oldconverged = converged;
      converged = c->PutORGate(converged, newconverged);
      delete oldconverged;
      if (k == 0)
        break;  // if k == 0 we are done, code below segfaults from nm

      // shift from bottom 2-by-2 minor
      for (int l_finder = 0; l_finder < k; l_finder++) {
        share* pl_finder = c->PutCONSGate((uint32_t*)&l_finder, bitlen);
        share* found_l = c->PutEQGate(priv_l, pl_finder);
        delete pl_finder;
        share* flanc = c->PutINVGate(converged);
        share* oldflanc = flanc;
        flanc = c->PutANDGate(found_l, flanc);
        delete oldflanc;
        delete found_l;
        share* oldx = x;
        x = c->PutMUXGate(w[l_finder], x, flanc);
        delete oldx;
        delete flanc;
      }

      nm = k - 1;
      y = new boolshare(w[nm]->get_wires(),
                        c);  // do not delete before assigning
      if (g)
        delete g;
      g = new boolshare(rv1[nm]->get_wires(), c);
      if (h)
        delete h;
      h = new boolshare(rv1[k]->get_wires(), c);
      {  // f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        share* ymz = c->PutFPGate(y, z, SUB, bitlen, 1, no_status);
        share* ypz = c->PutFPGate(y, z, ADD, bitlen, 1, no_status);
        share* temp1 = c->PutFPGate(ymz, ypz, MUL, bitlen, 1, no_status);
        share* gmh = c->PutFPGate(g, h, SUB, bitlen, 1, no_status);
        share* gph = c->PutFPGate(g, h, ADD, bitlen, 1, no_status);
        share* temp2 = c->PutFPGate(gmh, gph, MUL, bitlen, 1, no_status);
        share* num = c->PutFPGate(temp1, temp2, ADD, bitlen, 1, no_status);
        share* den1 = c->PutFPGate(two_gate, h, MUL, bitlen, 1, no_status);
        share* den2 = c->PutFPGate(den1, y, MUL, bitlen, 1, no_status);
        f = c->PutFPGate(num, den2, DIV, bitlen, 1, no_status);
        delete ymz;
        delete ypz;
        delete temp1;
        delete gmh;
        delete gph;
        delete temp2;
        delete num;
        delete den1;
        delete den2;
      }
      g = BuildPythagCircuit(f, one_gate, c, bitlen);
      {  // f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f)))- h)) / x;
        share* xmz = c->PutFPGate(x, z, SUB, bitlen, 1, no_status);
        share* xpz = c->PutFPGate(x, z, ADD, bitlen, 1, no_status);
        share* term1 = c->PutFPGate(xmz, xpz, MUL, bitlen, 1, no_status);
        share* sgf = new boolshare(g->get_wires(), c);
        BuildSignCircuit(sgf, f, c);
        share* fpsgf = c->PutFPGate(f, sgf, ADD, bitlen, 1, no_status);
        share* ydfpsgf = c->PutFPGate(y, fpsgf, DIV, bitlen, 1, no_status);
        share* uggg = c->PutFPGate(ydfpsgf, h, SUB, bitlen, 1, no_status);
        share* term2 = c->PutFPGate(h, uggg, MUL, bitlen, 1, no_status);
        share* num2 = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
        delete f;
        f = c->PutFPGate(num2, x, DIV, bitlen, 1, no_status);
        delete xmz;
        delete xpz;
        delete term1;
        delete sgf;
        delete fpsgf;
        delete ydfpsgf;
        delete uggg;
        delete term2;
        delete num2;
      }
      cc = new boolshare(one_gate->get_wires(), c);
      s = new boolshare(one_gate->get_wires(), c);

      for (int j = 0; j <= nm; j++) {
        share* priv_j = c->PutCONSGate((uint32_t*)&j, bitlen);
        share* jgtl = c->PutGTGate(priv_l, priv_j);  // have to invert for >=
        share* oldjgtl = jgtl;
        jgtl = c->PutINVGate(jgtl);
        delete oldjgtl;
        share* nconverged = c->PutINVGate(converged);
        share* modify = c->PutANDGate(jgtl, nconverged);
        delete nconverged;
        delete jgtl;
        delete priv_j;

        int i = j + 1;
        delete g;
        g = new boolshare(rv1[i]->get_wires(), c);
        delete y;
        y = new boolshare(w[i]->get_wires(), c);
        delete h;
        h = c->PutFPGate(s, g, MUL, bitlen, 1, no_status);
        share* oldg = g;
        g = c->PutFPGate(cc, g, MUL, bitlen, 1, no_status);
        delete oldg;
        delete z;
        z = BuildPythagCircuit(f, h, c, bitlen);
        share* oldr = rv1[j];
        rv1[j] = c->PutMUXGate(z, rv1[j], modify);
        delete oldr;
        temp = c->PutFPGate(f, z, DIV, bitlen, 1, no_status);
        share* oldcc = cc;
        cc = c->PutMUXGate(temp, cc, modify);
        delete oldcc;
        delete temp;
        share* olds = s;
        temp = c->PutFPGate(h, z, DIV, bitlen, 1, no_status);
        s = c->PutMUXGate(temp, s, modify);
        delete olds;
        delete temp;
        {  // f = x * c + g * s;
          share* term1 = c->PutFPGate(x, cc, MUL, bitlen, 1, no_status);
          share* term2 = c->PutFPGate(g, s, MUL, bitlen, 1, no_status);
          share* added = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
          share* oldf = f;
          f = c->PutMUXGate(added, f, modify);
          delete oldf;
          delete added;
          delete term1;
          delete term2;
        }
        {  // g = g * c - x * s;
          share* term1 = c->PutFPGate(g, cc, MUL, bitlen, 1, no_status);
          share* term2 = c->PutFPGate(x, s, MUL, bitlen, 1, no_status);
          share* oldg = g;
          g = c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
          delete oldg;
          delete term1;
          delete term2;
        }
        h = c->PutFPGate(y, s, MUL, bitlen, 1, no_status);
        y = c->PutFPGate(y, cc, MUL, bitlen, 1, no_status);
        for (jj = 0; jj < nCols; jj++) {
          share* oldx = x;
          x = c->PutMUXGate(v[jj][j], x, modify);
          delete oldx;
          z = new boolshare(v[jj][i]->get_wires(), c);
          {  // v[jj][j] = x * c + z * s;
            share* term1 = c->PutFPGate(x, cc, MUL, bitlen, 1, no_status);
            share* term2 = c->PutFPGate(z, s, MUL, bitlen, 1, no_status);
            share* added =
                c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
            share* oldv = v[jj][j];
            v[jj][j] = c->PutMUXGate(added, v[jj][j], modify);
            delete oldv;
            delete added;
            delete term1;
            delete term2;
          }
          {  // v[jj][i] = z * c - x * s;
            share* term1 = c->PutFPGate(z, cc, MUL, bitlen, 1, no_status);
            share* term2 = c->PutFPGate(x, s, MUL, bitlen, 1, no_status);
            share* subbed =
                c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
            share* oldv = v[jj][i];
            v[jj][i] = c->PutMUXGate(subbed, v[jj][i], modify);
            delete oldv;
            delete subbed;
            delete term1;
            delete term2;
          }
        }
        z = BuildPythagCircuit(f, h, c, bitlen);
        share* oldw = w[j];
        w[j] = c->PutMUXGate(z, w[j], modify);
        delete oldw;
        {  // if (z) {
          share* temp;
          share* ifz = z->get_wire_ids_as_share(0);
          for (uint32_t tt = 1; tt < bitlen - 1;
               tt++) {  // -1 -> do not include sign bit
            share* tempbit = z->get_wire_ids_as_share(tt);
            temp = ifz;
            ifz = c->PutORGate(temp, tempbit);
            delete tempbit;
            delete temp;
          }
          share* ifzam = c->PutANDGate(ifz, modify);
          delete ifz;

          temp = c->PutFPGate(one_gate, z, DIV, bitlen, 1, no_status);
          z = c->PutMUXGate(temp, z, ifzam);
          delete temp;

          temp = c->PutFPGate(f, z, MUL, bitlen, 1, no_status);
          cc = c->PutMUXGate(temp, cc, ifzam);
          delete temp;

          temp = c->PutFPGate(h, z, MUL, bitlen, 1, no_status);
          s = c->PutMUXGate(temp, s, ifzam);
          delete temp;
        }
        {  // f = c * g + s * y;
          share* term1 = c->PutFPGate(cc, g, MUL, bitlen, 1, no_status);
          share* term2 = c->PutFPGate(s, y, MUL, bitlen, 1, no_status);
          share* added = c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
          share* oldf = f;
          f = c->PutMUXGate(added, f, modify);
          delete oldf;
          delete added;
          delete term1;
          delete term2;
        }
        {  // x = c * y - s * g;
          share* term1 = c->PutFPGate(cc, y, MUL, bitlen, 1, no_status);
          share* term2 = c->PutFPGate(s, g, MUL, bitlen, 1, no_status);
          share* subbed = c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
          share* oldx = x;
          x = c->PutMUXGate(subbed, x, modify);
          delete oldx;
          delete subbed;
          delete term1;
          delete term2;
        }
        for (jj = 0; jj < nRows; jj++) {
          y = new boolshare(a[jj][j]->get_wires(), c);
          z = new boolshare(a[jj][i]->get_wires(), c);
          {  // a[jj][j] = y * c + z * s;
            share* term1 = c->PutFPGate(y, cc, MUL, bitlen, 1, no_status);
            share* term2 = c->PutFPGate(z, s, MUL, bitlen, 1, no_status);
            share* added =
                c->PutFPGate(term1, term2, ADD, bitlen, 1, no_status);
            share* olda = a[jj][j];
            a[jj][j] = c->PutMUXGate(added, a[jj][j], modify);
            delete olda;
            delete added;
            delete term1;
            delete term2;
          }
          {  // a[jj][i] = z * c - y * s;
            share* term1 = c->PutFPGate(z, cc, MUL, bitlen, 1, no_status);
            share* term2 = c->PutFPGate(y, s, MUL, bitlen, 1, no_status);
            share* subbed =
                c->PutFPGate(term1, term2, SUB, bitlen, 1, no_status);
            share* olda = a[jj][i];
            a[jj][i] = c->PutMUXGate(subbed, a[jj][i], modify);
            delete olda;
            delete subbed;
            delete term1;
            delete term2;
          }
        }
      }
      for (int l_finder = 0; l_finder < k; l_finder++) {
        share* pl_finder = c->PutCONSGate((uint32_t*)&l_finder, 32);
        share* found_l = c->PutEQGate(priv_l, pl_finder);
        delete pl_finder;
        share* flanc = c->PutINVGate(converged);
        share* oldflanc = flanc;
        flanc = c->PutANDGate(found_l, flanc);
        delete oldflanc;
        delete found_l;
        share* oldr = rv1[l_finder];
        rv1[l_finder] = c->PutMUXGate(zero_gate, rv1[l_finder], flanc);
        delete oldr;
        delete flanc;
      }
      share* nconverged = c->PutINVGate(converged);
      share* oldr = rv1[k];
      rv1[k] = c->PutMUXGate(f, rv1[k], nconverged);
      delete oldr;
      share* oldw = w[k];
      w[k] = c->PutMUXGate(x, w[k], nconverged);
      delete oldw;
      delete nconverged;
    }
    delete converged;
    printf("/");
    fflush(stdout);
  }
}

void BuildAndRunSvd(share*** s_a, int nRows, int nCols, share** s_w,
                    share*** s_v, BooleanCircuit* c, ABYParty* party,
                    e_role role, std::function<void()> toRawShares,
                    std::function<void()> toShareObjects) {
#if PPL_FLOW == PPL_FLOW_DO  // set dx to zero if no error
  (void)toRawShares;
  (void)toShareObjects;
  BuildSvdDO(s_a, nRows, nCols, s_w, s_v, c, party, role);
#elif PPL_FLOW == PPL_FLOW_LOOP_LEAK
  BuildAndRunSvdLoopLeak(s_a, nRows, nCols, s_w, s_v, c, party, role,
                         toRawShares, toShareObjects);
#elif PPL_FLOW == PPL_FLOW_SiSL
  BuildSvdOSL(s_a, nRows, nCols, s_w, s_v, c, party, role);
#endif
}

uint32_t test_fabs_circuit(e_role role, const std::string& address,
                           uint16_t port, seclvl seclvl, uint32_t nthreads,
                           e_mt_gen_alg mt_alg, e_sharing sharing, float a) {
  uint32_t bitlen = 32;

  uint32_t reservegates = 65536;
  const std::string& abycircdir = "../../extern/ABY/bin/circ";
  ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                                 mt_alg, reservegates, abycircdir);

  std::vector<Sharing*>& sharings = party->GetSharings();

  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

  share* s_a;

  if (role == SERVER) {
    s_a = circ->PutINGate((uint32_t*)&a, 32, role);
  } else {
    s_a = circ->PutDummyINGate(32);
  }

  share* out;
  BuildFabsCircuit(s_a, (BooleanCircuit*)circ);

  out = circ->PutOUTGate(s_a, ALL);

  CLOCK(ExecCircuit);
  TIC(ExecCircuit);
  party->ExecCircuit();
  TOC(ExecCircuit);

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;

  // This method only works for an output length of maximum 64 bits in general,
  // if the output length is higher you must use get_clear_value_ptr
  out->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
  // R[i]->get_clear_value_ptr();

  // cout << out_bitlen << "jfjf" << endl;
  cout << *(float*)output << endl;

  delete party;
  // delete[] A;
  // delete[] B;
  // delete[] s_A;
  // delete[] s_B;
  // delete[] s_out;
  return 0;
}

uint32_t test_fmax_circuit(e_role role, const std::string& address,
                           uint16_t port, seclvl seclvl, uint32_t nthreads,
                           e_mt_gen_alg mt_alg, e_sharing sharing, float a,
                           float b) {
  uint32_t bitlen = 32;

  uint32_t reservegates = 65536;
  const std::string& abycircdir = "../../extern/ABY/bin/circ";
  ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                                 mt_alg, reservegates, abycircdir);

  std::vector<Sharing*>& sharings = party->GetSharings();

  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

  share* s_a;
  share* s_b;

  if (role == SERVER) {
    s_a = circ->PutINGate((uint32_t*)&a, 32, role);
    s_b = circ->PutINGate((uint32_t*)&b, 32, role);
  } else {
    s_a = circ->PutDummyINGate(32);
    s_b = circ->PutDummyINGate(32);
  }

  share* s_out;
  share* out;
  s_out = BuildFMAXCircuit(s_a, s_b, (BooleanCircuit*)circ);

  out = circ->PutOUTGate(s_out, ALL);

  CLOCK(ExecCircuit);
  TIC(ExecCircuit);
  party->ExecCircuit();
  TOC(ExecCircuit);

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;

  // This method only works for an output length of maximum 64 bits in general,
  // if the output length is higher you must use get_clear_value_ptr
  out->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
  // R[i]->get_clear_value_ptr();

  // cout << out_bitlen << "jfjf" << endl;
  cout << *(float*)output << endl;

  delete party;
  // delete[] A;
  // delete[] B;
  // delete[] s_A;
  // delete[] s_B;
  // delete[] s_out;
  return 0;
}

uint32_t test_negative_circuit(e_role role, const std::string& address,
                               uint16_t port, seclvl seclvl, uint32_t nthreads,
                               e_mt_gen_alg mt_alg, e_sharing sharing,
                               float a) {
  uint32_t bitlen = 32;

  uint32_t reservegates = 65536;
  const std::string& abycircdir = "../../extern/ABY/bin/circ";
  ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                                 mt_alg, reservegates, abycircdir);

  std::vector<Sharing*>& sharings = party->GetSharings();

  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

  share* s_a;

  if (role == SERVER) {
    s_a = circ->PutINGate((uint32_t*)&a, 32, role);
  } else {
    s_a = circ->PutDummyINGate(32);
  }

  // share* s_out;
  share* out;
  BuildNegativeCircuit(s_a, (BooleanCircuit*)circ);

  out = circ->PutOUTGate(s_a, ALL);

  CLOCK(ExecCircuit);
  TIC(ExecCircuit);
  party->ExecCircuit();
  TOC(ExecCircuit);

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;

  // This method only works for an output length of maximum 64 bits in general,
  // if the output length is higher you must use get_clear_value_ptr
  out->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
  // R[i]->get_clear_value_ptr();

  // cout << out_bitlen << "jfjf" << endl;
  cout << *(float*)output << endl;

  delete party;
  // delete[] A;
  // delete[] B;
  // delete[] s_A;
  // delete[] s_B;
  // delete[] s_out;
  return 0;
}
