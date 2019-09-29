#include <opencv2/opencv.hpp>

#include <abycore/aby/abyparty.h>
#include <abycore/circuit/arithmeticcircuits.h>
#include <abycore/circuit/booleancircuits.h>
#include <abycore/circuit/circuit.h>
#include <abycore/sharing/sharing.h>

using std::vector;

const std::string address = "127.0.0.1";
constexpr const int port = 7766;
constexpr const uint32_t nthreads = 1;
constexpr const e_mt_gen_alg mt_alg = MT_OT;
const seclvl slvl = LT;  // get_sec_lvl(secparam);
const vector<e_sharing> ctypes = {/*S_ARITH,*/ S_BOOL,
                                  S_YAO};  // no floats with arith
std::map<e_sharing, std::string> cnames = {
    {S_ARITH, "arith"},
    {S_BOOL, "bool"},
    {S_YAO, "yao"},
};
constexpr const uint32_t reservegates = 65536;
constexpr const uint32_t bitlen = 32;

// TODO(jc): should return pair<bool, int>
typedef void (*aby_localizer)(share* s_threeDPts[], share* s_twoDPts[],
                              int numPts, share* s_f, share* s_cx, share* s_cy,
                              share* s_x[], BooleanCircuit* c, ABYParty* party,
                              e_role role);

static std::string get_circuit_dir() {
  char result[PATH_MAX];
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  const char* path;
  if (count != -1) {
    path = dirname(result);
  } else {
    std::cerr << "cant find circuit path\n";
    assert(false);
  }
  std::string base_path = string(path) + "/../../extern/ABY/bin/circ";
  return base_path;
}

void aby_localize_wrapper(e_role role, e_sharing sharing, cv::Mat rvec,
                          cv::Mat tvec, cv::Mat cameraMatrix,
                          cv::Mat distCoeffs [[maybe_unused]],
                          vector<cv::Point3f> objectPoints,
                          vector<cv::Point2f> imagePoints, vector<float>& res,
                          auto secure_localize_func) {
  float f = cameraMatrix.at<float>(0, 0);
  float cx = cameraMatrix.at<float>(0, 2);
  float cy = cameraMatrix.at<float>(1, 2);

  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

#if PPL_FLOW == PPL_FLOW_LOOP_LEAK || PPL_FLOW == PPL_FLOW_SiSL
  //  always use boolean shares (not yao) so BuildAndRunLM can get raw shares.
  //  this is because you cant use Y2B shares on Yao input gates.
  Circuit* bc = sharings[S_BOOL]->GetCircuitBuildRoutine();
#else
  Circuit* bc = circ;
#endif

  assert(objectPoints.size() == imagePoints.size());
  int numPts = objectPoints.size();

  // Allocate space for shares
  share** s_objectPoints = new share*[4 * numPts];
  share** s_imagePoints = new share*[3 * numPts];
  share* s_f;
  share* s_cx;
  share* s_cy;
  share** s_x = new share*[6];

  // Prepare inputs
  if (role == SERVER) {
    // float one=1.0;
    for (int p = 0; p < numPts; p++) {
      s_objectPoints[p] =
          bc->PutINGate((uint32_t*)&objectPoints[p].x, bitlen, role);
    }
    for (int p = 0; p < numPts; p++) {
      s_objectPoints[numPts + p] =
          bc->PutINGate((uint32_t*)&objectPoints[p].y, bitlen, role);
    }
    for (int p = 0; p < numPts; p++) {
      s_objectPoints[(2 * numPts) + p] =
          bc->PutINGate((uint32_t*)&objectPoints[p].z, bitlen, role);
    }
    // for(int p=0; p<numPts; p++) {
    //     s_objectPoints[3*numPts+p] = bc->PutINGate((uint32_t*) &one, bitlen,
    //     role);
    // }
    for (int p = 0; p < numPts; p++) {
      s_imagePoints[p] =
          bc->PutINGate((uint32_t*)&imagePoints[p].x, bitlen, role);
    }
    for (int p = 0; p < numPts; p++) {
      s_imagePoints[numPts + p] =
          bc->PutINGate((uint32_t*)&imagePoints[p].y, bitlen, role);
    }
    // for(int p=0; p<numPts; p++) {
    //     s_imagePoints[2*numPts+p] = bc->PutINGate((uint32_t*) &one, bitlen,
    //     role);
    // }
    s_f = bc->PutINGate((uint32_t*)&f, bitlen, role);
    s_cx = bc->PutINGate((uint32_t*)&cx, bitlen, role);
    s_cy = bc->PutINGate((uint32_t*)&cy, bitlen, role);
    for (int p = 0; p < 3; p++) {
      s_x[p] = bc->PutINGate((uint32_t*)&rvec.at<float>(p), bitlen, role);
    }
    for (int p = 0; p < 3; p++) {
      s_x[p + 3] = bc->PutINGate((uint32_t*)&tvec.at<float>(p), bitlen, role);
    }
  } else {
    for (int p = 0; p < 3 * numPts; p++) {  // ignore constant 1s
      s_objectPoints[p] = bc->PutDummyINGate(bitlen);
    }
    for (int p = 0; p < 2 * numPts; p++) {  // ignore constant 1s
      s_imagePoints[p] = bc->PutDummyINGate(bitlen);
    }
    s_f = bc->PutDummyINGate(bitlen);
    s_cx = bc->PutDummyINGate(bitlen);
    s_cy = bc->PutDummyINGate(bitlen);
    for (int p = 0; p < 6; p++) {
      s_x[p] = bc->PutDummyINGate(bitlen);
    }
  }

  secure_localize_func(s_objectPoints, s_imagePoints, numPts, s_f, s_cx, s_cy,
                       s_x, (BooleanCircuit*)circ, party, role);

  for (int i = 0; i < 6; i++) {
    share* temp = s_x[i];
    s_x[i] = bc->PutOUTGate(s_x[i], ALL);
    delete temp;
  }

  party->ExecCircuit();

  if (role == SERVER) {
    cout << "secure result: ";
  }
  for (int i = 0; i < 6; i++) {
    uint32_t* output;
    uint32_t out_bitlen, out_nvals;
    s_x[i]->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
    if (role == SERVER) {
      cout << *(float*)output << " ";
    }
    res[i] = *(float*)output;
  }
  if (role == SERVER) {
    cout << '\n';
  }

  delete party;
}