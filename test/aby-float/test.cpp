#include <opencv2/opencv.hpp>

#include <abycore/aby/abyparty.h>
#include <abycore/circuit/arithmeticcircuits.h>
#include <abycore/circuit/booleancircuits.h>
#include <abycore/circuit/circuit.h>
#include <abycore/sharing/sharing.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <libgen.h>        // dirname
#include <linux/limits.h>  // PATH_MAX
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <test_params.h>
#include <unistd.h>  // readlink
#include <iostream>
#include <string>
#include <thread>

#include <gaussnewtonlocalization.h>
#include <invert.h>
#include <jlog.h>
#include <lmlocalization.h>
#include <matmult.h>
#include <projectpoints.h>
#include <rodrigues.h>
#include <svd.h>
#include <trigfuncs.h>
#include <twonormsq.h>
#include <util.h>
#include <cleartext-ref/gaussnewtonlocalization.hpp>
#include <cleartext-ref/invert.hpp>
#include <cleartext-ref/lmlocalization.hpp>
#include <cleartext-ref/matmult.hpp>
#include <cleartext-ref/projectpoints.hpp>
#include <cleartext-ref/rodrigues.hpp>
#include <cleartext-ref/svd.hpp>
#include <cleartext-ref/trigfuncs.hpp>
#include <cleartext-ref/twonormsq.hpp>
#include <hoff_features.hpp>
#include <localize_wrapper.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using std::cout;
using std::map;
using std::vector;

enum class TRIG_FUNC : uint8_t {
  SIN,
  COS,
};

void aby_trig(e_role role, e_sharing sharing, float angle, TRIG_FUNC tf) {
  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

  share* s_in;
  if (role == SERVER) {
    s_in = circ->PutINGate((uint32_t*)&angle, 32, role);
  } else {
    s_in = circ->PutDummyINGate(32);
  }

  share* s_res = NULL;
  share* res = NULL;  // stores plaintext
  float cleartext_res;
  if (tf == TRIG_FUNC::SIN) {
    s_res = BuildSinCircuit(s_in, (BooleanCircuit*)circ);
    cleartext_res = sin(angle);
  } else if (tf == TRIG_FUNC::COS) {
    s_res = BuildCosCircuit(s_in, (BooleanCircuit*)circ);
    cleartext_res = cos(angle);
  } else {
    assert(false);
  }

  res = circ->PutOUTGate(s_res, ALL);
  party->ExecCircuit();
  assert(res != NULL);

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;
  res->get_clear_value_vec(&output, &out_bitlen, &out_nvals);

  REQUIRE(out_bitlen == 32);
  REQUIRE(out_nvals == 1);
  REQUIRE_THAT(*(float*)output, WithinAbs(cleartext_res, float_test_epsilon));

  delete s_in;
  delete s_res;
  delete res;
  delete party;
}

TEST_CASE("ABY trig functions are computed", "[aby_trig]") {
  float angle = .2;

  std::thread bob([angle]() {
    for (auto tfunc : {TRIG_FUNC::SIN, TRIG_FUNC::COS}) {
      for (auto ctype : ctypes) {
        aby_trig(CLIENT, ctype, angle, tfunc);
        aby_trig(CLIENT, ctype, angle, tfunc);
      }
    }
  });

  for (auto tfunc : {TRIG_FUNC::SIN, TRIG_FUNC::COS}) {
    for (auto ctype : ctypes) {
      aby_trig(SERVER, ctype, angle, tfunc);
      aby_trig(SERVER, ctype, angle, tfunc);
    }
  }
  bob.join();
}

void aby_matmult(e_role role, e_sharing sharing, cv::Mat A, cv::Mat B) {
  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

  int m = A.rows, n = A.cols, mm = B.rows, nn = B.cols;
  share *s_A[m * n], *s_B[mm * nn];

  for (int i = 0; i < m * n; ++i) {
    if (role == SERVER)
      s_A[i] = circ->PutINGate((uint32_t*)&A.at<uint32_t>(i), 32, role);
    else
      s_A[i] = circ->PutDummyINGate(32);
  }

  for (int i = 0; i < mm * nn; ++i) {
    if (role == SERVER)
      s_B[i] = circ->PutINGate((uint32_t*)&B.at<uint32_t>(i), 32, role);
    else
      s_B[i] = circ->PutDummyINGate(32);
  }

  share* res[m * nn];
  share* s_out[m * nn];  // stores plaintext
  BuildMatmultCircuit(s_A, m, n, s_B, mm, nn, res, (BooleanCircuit*)circ);

  for (int i = 0; i < m * nn; ++i) {
    s_out[i] = circ->PutOUTGate(res[i], ALL);
  }

  party->ExecCircuit();

  cv::Mat C = cv::Mat::zeros(m, nn, cv::DataType<float>::type);
  matmult(A.ptr<float>(), A.rows, A.cols, B.ptr<float>(), B.rows, B.cols,
          C.ptr<float>());

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;
  for (int i = 0; i < m * nn; i++) {
    s_out[i]->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
    REQUIRE(out_bitlen == 32);
    REQUIRE(out_nvals == 1);
    REQUIRE_THAT(*(float*)output,
                 WithinAbs(C.ptr<float>()[i], float_test_epsilon));
  }
  delete party;
}

TEST_CASE("ABY matrix multiply is computed", "[aby_mat]") {
  cv::Mat A = cv::Mat::ones(3, 5, cv::DataType<float>::type);
  cv::Mat B = cv::Mat::ones(5, 4, cv::DataType<float>::type) * 2;

  std::thread bob([A, B]() {
    for (auto ctype : ctypes) {
      aby_matmult(CLIENT, ctype, A, B);
    }
  });

  for (auto ctype : ctypes) {
    aby_matmult(SERVER, ctype, A, B);
  }
  bob.join();
}

void aby_rod(e_role role, e_sharing sharing, float* r) {
  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

  share* s_r[3];

  if (role == SERVER) {
    s_r[0] = circ->PutINGate((uint32_t*)&r[0], 32, role);
    s_r[1] = circ->PutINGate((uint32_t*)&r[1], 32, role);
    s_r[2] = circ->PutINGate((uint32_t*)&r[2], 32, role);
  } else {
    s_r[0] = circ->PutDummyINGate(32);
    s_r[1] = circ->PutDummyINGate(32);
    s_r[2] = circ->PutDummyINGate(32);
  }

  share* s_R[12];
  share* R[12];  // stores plaintext
  BuildRodriguesCircuit(s_r, s_R, (BooleanCircuit*)circ);

  for (int i = 0; i < 12; ++i) {
    if (i == 3 || i == 7 || i == 11)
      continue;
    R[i] = circ->PutOUTGate(s_R[i], ALL);
  }

  party->ExecCircuit();

  cv::Mat clear_res = cv::Mat::zeros(3, 4, cv::DataType<float>::type);
  rodrigues(r, clear_res.ptr<float>());

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;
  for (int i = 0; i < 12; i++) {
    if (i == 3 || i == 7 || i == 11)
      continue;
    assert(R[i] != NULL);
    R[i]->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
    REQUIRE(out_bitlen == 32);
    REQUIRE(out_nvals == 1);
    REQUIRE_THAT(*(float*)output,
                 WithinAbs(clear_res.ptr<float>()[i], float_test_epsilon));
  }

  delete party;
}

TEST_CASE("ABY rodriguez transformation is computed", "[aby_rod]") {
  float r[3] = {.2, .4, .2};  // 3x1

  std::thread bob([&r]() {
    for (auto ctype : ctypes) {
      aby_rod(CLIENT, ctype, r);
    }
  });

  for (auto ctype : ctypes) {
    aby_rod(SERVER, ctype, r);
  }
  bob.join();
}

// From points, first estimate the pose using opencv.
// Then reproject using the estimated pose
// and compare to the origional points.
void aby_proj(e_role role, e_sharing sharing, cv::Mat cameraMatrix,
              cv::Mat distCoeffs, vector<cv::Point2f> imagePoints,
              vector<cv::Point3f> objectPoints) {
  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();
  int numPoints = imagePoints.size();

  // Find rotation and translation
  cv::Mat rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::Mat tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs, rvec, tvec,
               false, cv::SOLVEPNP_ITERATIVE);

  cv::Mat x = cv::Mat(6, 1, cv::DataType<float>::type);
  rvec.copyTo(x(cv::Range(0, 3), cv::Range(0, 1)));
  tvec.copyTo(x(cv::Range(3, 6), cv::Range(0, 1)));

  // cleartext projection
  cv::Mat matOP = cv::Mat(objectPoints, false);
  matOP = matOP.reshape(1, matOP.total());
  matOP.convertTo(matOP, cv::DataType<float>::type);
  cv::Mat P = cv::Mat::ones(4, numPoints, cv::DataType<float>::type);
  P(cv::Range(0, 3), cv::Range(0, numPoints)) =
      matOP.t();  // copy points into new matrix
  cv::Mat projected = cv::Mat(3, numPoints, cv::DataType<float>::type);
  projectPoints(P.ptr<float>(), x.ptr<float>(), cameraMatrix.ptr<float>(),
                projected.ptr<float>(), numPoints);
  for (int i = 0; i < numPoints; ++i) {
    REQUIRE_THAT(imagePoints[i].x,
                 WithinAbs(projected.at<float>(0, i), reprojection_tol));
    REQUIRE_THAT(imagePoints[i].y,
                 WithinAbs(projected.at<float>(1, i), reprojection_tol));
  }

  share* s_P[4 * numPoints];
  share* s_x[6];
  share* s_K[9];

  float one = 1;
  share* one_gate = circ->PutCONSGate((uint32_t*)&one, bitlen);

  if (role == SERVER) {
    for (int i = 0; i < 9; i++) {
      s_K[i] = circ->PutINGate((uint32_t*)&cameraMatrix.at<float>(i), 32, role);
    }
    for (int i = 0; i < 6; i++) {
      s_x[i] = circ->PutINGate((uint32_t*)&x.at<float>(i), 32, role);
    }
    for (int i = 0; i < 3 * numPoints; i++) {
      s_P[i] = circ->PutINGate((uint32_t*)&P.at<float>(i), 32, role);
    }
    // in gate for contant 1s
    for (int i = 3 * numPoints; i < 4 * numPoints; i++) {
      s_P[i] = circ->PutCONSGate((uint32_t*)&one, bitlen);
    }
  } else {
    for (int i = 0; i < 9; i++) {
      s_K[i] = circ->PutDummyINGate(32);
    }
    for (int i = 0; i < 6; i++) {
      s_x[i] = circ->PutDummyINGate(32);
    }
    for (int i = 0; i < 3 * numPoints; i++) {
      s_P[i] = circ->PutDummyINGate(32);
    }
    // in gate for contant 1s
    for (int i = 3 * numPoints; i < 4 * numPoints; i++) {
      s_P[i] = circ->PutCONSGate((uint32_t*)&one, bitlen);
    }
  }

  share*
      s_res[3 * numPoints];   // dont forget extra space for homog constant 1's
  share* res[3 * numPoints];  // stores plaintext
  BuildProjectPointsCircuit(s_P, s_x, s_K, s_res, numPoints,
                            (BooleanCircuit*)circ);

  for (int i = 0; i < 2 * numPoints; ++i) {
    res[i] = circ->PutOUTGate(s_res[i], ALL);
  }

  party->ExecCircuit();

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;
  for (int i = 0; i < numPoints; ++i) {
    assert(res[i] != NULL);
    res[i]->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
    REQUIRE(out_bitlen == 32);
    REQUIRE(out_nvals == 1);
    REQUIRE_THAT(*(float*)output,
                 WithinAbs(imagePoints[i].x, reprojection_tol));

    assert(res[i + numPoints] != NULL);
    res[i + numPoints]->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
    REQUIRE(out_bitlen == 32);
    REQUIRE(out_nvals == 1);
    REQUIRE_THAT(*(float*)output,
                 WithinAbs(imagePoints[i].y, reprojection_tol));
  }

  delete party;
}

TEST_CASE("ABY point projection is computed", "[aby_proj]") {
  float f = 715;
  float cx = 354;
  float cy = 245;
  float _cM[] = {f, 0, cx, 0, f, cy, 0, 0, 1};
  cv::Mat cameraMatrix = cv::Mat(3, 3, cv::DataType<float>::type, _cM);
  cv::Mat distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
  cv::Mat rvec(3, 1, cv::DataType<float>::type);
  cv::Mat tvec(3, 1, cv::DataType<float>::type);
  vector<cv::Point2f> imagePoints;
  vector<cv::Point3f> objectPoints;
  Hoffs2DPoints(imagePoints);
  Hoffs3DPoints(objectPoints);

  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_proj(CLIENT, ctype, cameraMatrix, distCoeffs, imagePoints,
               objectPoints);
    }
  });

  for (auto ctype : ctypes) {
    aby_proj(SERVER, ctype, cameraMatrix, distCoeffs, imagePoints,
             objectPoints);
  }
  bob.join();
}

void aby_svd_sign(e_role role, e_sharing sharing, float a, float b) {
  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
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

  share *out_a, *out_b;
  BuildSignCircuit(s_a, s_b, (BooleanCircuit*)circ);
  out_a = circ->PutOUTGate(s_a, ALL);
  out_b = circ->PutOUTGate(s_b, ALL);
  party->ExecCircuit();

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;

  float cleartext_res = mysign(a, b);

  out_a->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
  REQUIRE(out_bitlen == 32);
  REQUIRE(out_nvals == 1);
  REQUIRE_THAT(*(float*)output, WithinAbs(cleartext_res, float_test_epsilon));

  out_b->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
  REQUIRE(out_bitlen == 32);
  REQUIRE(out_nvals == 1);
  REQUIRE_THAT(*(float*)output, WithinAbs(b, float_test_epsilon));

  delete party;
}

TEST_CASE("ABY SVD sign function is computed", "[aby_svd]") {
  float a = 2.40, b = 0;
  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_svd_sign(CLIENT, ctype, a, b);
    }
  });

  for (auto ctype : ctypes) {
    aby_svd_sign(SERVER, ctype, a, b);
  }
  bob.join();
}

void aby_svd_pythag(e_role role, e_sharing sharing, float a, float b) {
  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
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
  s_out = BuildPythagCircuit(s_a, s_b, (BooleanCircuit*)circ, bitlen);

  out = circ->PutOUTGate(s_out, ALL);

  party->ExecCircuit();

  float cleartext_res = mypythag(a, b);

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;
  out->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
  REQUIRE(out_bitlen == 32);
  REQUIRE(out_nvals == 1);
  REQUIRE_THAT(*(float*)output, WithinAbs(cleartext_res, float_test_epsilon));

  delete party;
}

TEST_CASE("ABY SVD pythag function is computed", "[aby_svd]") {
  float a = -2.40, b = 96.0;
  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_svd_pythag(CLIENT, ctype, a, b);
    }
  });

  for (auto ctype : ctypes) {
    aby_svd_pythag(SERVER, ctype, a, b);
  }
  bob.join();
}

void aby_svd(e_role role, e_sharing sharing, float** in, int m, int n) {
  // OpenCV SVD - requires linearize
  float* cvin = new float[m * n];
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      cvin[i * n + j] = in[i][j];
    }
  }
  cv::Mat cvinM = cv::Mat(m, n, cv::DataType<float>::type, in);
  cv::SVD cvsvd(cvinM, cv::SVD::FULL_UV);  // constructor
  if (role == SERVER) {
    cout << "opencv result\n"
         << cvsvd.u << '\n'
         << cvsvd.w << '\n'
         << cvsvd.vt << '\n';
  }
  delete[] cvin;

  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

#if PPL_FLOW == PPL_FLOW_LOOP_LEAK
  //  always use boolean shares (not yao) so BuildAndRunSvd can get raw shares.
  //  this is because you cant use Y2B shares on Yao input gates.
  Circuit* bc = sharings[S_BOOL]->GetCircuitBuildRoutine();
#else
  Circuit* bc = circ;
#endif

  // cleartext SVD - svd overwrites with u matrix, so copy
  float** a = new float*[m];
  for (int i = 0; i < m; i++) {
    a[i] = new float[n];
    for (int j = 0; j < n; j++) {
      a[i][j] = in[i][j];
    }
  }
  float* w = new float[n]();  // aka sigma, only diag
  float** v = new float*[n];  // nxn
  for (int i = 0; i < n; i++)
    v[i] = new float[n]();

  svdcmp(a, m, n, w, v);
  if (role == SERVER) {
    cout << "cleartext" << '\n';
    printMatrix("a", a, m, n);
    printVector("w", w, n);
    printMatrix("v", v, n, n);
  }

  share*** s_a = new share**[m];
  for (int i = 0; i < m; i++) {
    s_a[i] = new share*[n];
  }
  share** s_w = new share*[n];
  share*** s_v = new share**[n];
  for (int i = 0; i < n; i++) {
    s_v[i] = new share*[n];
  }

  // create shared output from plaintext input
  //  initialize input shares
  //  always use boolean shares (not yao) so BuildAndRunSvd can get raw shares.
  //  this is because you cant use Y2B shares on Yao input gates.
  if (role == SERVER) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        s_a[i][j] = bc->PutINGate((uint32_t*)&in[i][j], 32, role);
      }
    }
  } else {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        s_a[i][j] = bc->PutDummyINGate(32);
      }
    }
  }

  BuildAndRunSvd(s_a, m, n, s_w, s_v, (BooleanCircuit*)circ, party, role);

  // shared outputs from circuit
  //  again always use boolean shares (not yao) because output is coming
  //  directly from a boolean circuit ->PutSharedInGate()
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      share* temp = s_a[i][j];
      s_a[i][j] = circ->PutOUTGate(temp, ALL);
      delete temp;
    }
  }
  for (int i = 0; i < n; i++) {
    share* temp = s_w[i];
    s_w[i] = circ->PutOUTGate(temp, ALL);
    delete temp;
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      share* temp = s_v[i][j];
      s_v[i][j] = circ->PutOUTGate(temp, ALL);
      delete temp;
    }
  }
  // run the circuit
  party->ExecCircuit();

  // get cleartext output
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      // compare cleartext to opencv
      // REQUIRE_THAT(a[i][j], WithinAbs(cvsvd.u.at<float>(i,j), svd_tol));

      // compare secure to cleartext
      float* valptr = (float*)s_a[i][j]->get_clear_value_ptr();
      REQUIRE_THAT(*(float*)valptr, WithinAbs(a[i][j], svd_tol));
      // in[i][j] = *valptr;
      free(valptr);
      delete s_a[i][j];
    }
    delete[] s_a[i];
  }
  delete[] s_a;
  for (int i = 0; i < n; i++) {
    // compare cleartext to opencv
    // REQUIRE_THAT(w[i], WithinAbs(cvsvd.w.at<float>(i), svd_tol));

    float* valptr = (float*)s_w[i]->get_clear_value_ptr();
    REQUIRE_THAT(*(float*)valptr, WithinAbs(w[i], svd_tol));
    // w[i] = *valptr;
    free(valptr);
    delete s_w[i];
  }
  delete[] s_w;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // compare cleartext to opencv
      // REQUIRE_THAT(v[i][j], WithinAbs(cvsvd.vt.at<float>(i,j), svd_tol));

      float* valptr = (float*)s_v[i][j]->get_clear_value_ptr();
      REQUIRE_THAT(*(float*)valptr, WithinAbs(v[i][j], svd_tol));
      // v[i][j] = *valptr;
      free(valptr);
      delete s_v[i][j];
    }
    delete[] s_v[i];
  }
  delete[] s_v;
  delete[] w;
  delete[] v;
  delete party;
}

TEST_CASE("ABY large SVD function is computed", "[aby_svd]") {
  int m = 12, n = 6;
  float one_d_M[] = {
      -9.1552734, 149.53613,  -59.509277, 18.310547, 0,         3.0517578,
      -164.79492, -64.086914, -47.302246, 0,         19.836426, 1.5258789,
      27.46582,   77.819824,  50.354004,  22.888184, 0,         1.5258789,
      -56.45752,  -42.724609, -53.405762, 0,         22.888184, 3.0517578,
      45.776367,  79.345703,  79.345703,  21.362305, 0,         0,
      -32.806396, -13.73291,  -20.599365, 0,         24.414062, 4.5776367,
      -21.362305, 108.3374,   -93.078613, 16.784668, 0,         3.0517578,
      -152.58789, -45.776367, -10.681152, 0,         19.836426, 3.0517578,
      6.1035156,  39.672852,  3.0517578,  21.362305, 0,         0,
      -39.672852, -21.362305, -16.784668, 0,         22.888184, 1.5258789,
      24.414062,  48.828125,  27.46582,   24.414062, 0,         0,
      -13.73291,  6.1035156,  12.207031,  0,         22.888184, 1.5258789};
  float** in = new float*[m];
  for (int i = 0; i < m; i++) {
    in[i] = new float[n]();
    for (int j = 0; j < n; j++) {
      in[i][j] = one_d_M[i * n + j];
    }
  }

  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_svd(CLIENT, ctype, in, m, n);
    }
  });

  for (auto ctype : ctypes) {
    aby_svd(SERVER, ctype, in, m, n);
  }
  bob.join();
}

void aby_invert(e_role role, e_sharing sharing, float* in, int m, int n) {
  cv::Mat cvM = cv::Mat(m, n, cv::DataType<float>::type, in);
  cv::Mat cvres;
  invert(cvM, cvres, cv::DECOMP_SVD);
  cout << "opencv result:\n" << cvres << endl;

  float** M = new float*[m];
  for (int i = 0; i < m; i++) {
    M[i] = new float[n];
    for (int j = 0; j < n; j++) {
      M[i][j] = in[i * n + j];
    }
  }

  // res is nxm (not mxn)
  float** res = new float*[n];
  for (int i = 0; i < n; i++) {
    res[i] = new float[m];
  }

  myinvert(M, m, n, res);
  printMatrix("cleartext result:", res, n, m);

  // check cleartext vs opencv
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      REQUIRE_THAT(cvres.at<float>(i, j),
                   WithinAbs(res[i][j], float_test_epsilon));
    }
  }

  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

#if PPL_FLOW == PPL_FLOW_LOOP_LEAK
  //  always use boolean shares (not yao) so BuildAndRunSvd can get raw shares.
  //  this is because you cant use Y2B shares on Yao input gates.
  Circuit* bc = sharings[S_BOOL]->GetCircuitBuildRoutine();
#else
  Circuit* bc = circ;
#endif

  // share arrays
  share*** s_in = new share**[m];
  for (int i = 0; i < m; i++) {
    s_in[i] = new share*[n];
  }

  // initialize input shares
  //  always use boolean shares (not yao) so BuildAndRunSvd (called in
  //  BuildInvertCircuit) can get raw shares. This is because you cant use Y2B
  //  shares on Yao input gates.
  if (role == SERVER) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        s_in[i][j] = bc->PutINGate((uint32_t*)&M[i][j], 32, role);
      }
    }
  } else {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        s_in[i][j] = bc->PutDummyINGate(32);
      }
    }
  }

  // res is nxm not mxn
  share*** s_res = new share**[n];
  for (int i = 0; i < n; i++) {
    s_res[i] = new share*[m];
  }

  BuildInvertCircuit(s_in, m, n, s_res, (BooleanCircuit*)circ, party, role);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      share* temp = s_res[i][j];
      s_res[i][j] = circ->PutOUTGate(s_res[i][j], ALL);
      delete temp;
    }
  }

  party->ExecCircuit();

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      uint32_t* output;
      uint32_t out_bitlen, out_nvals;
      s_res[i][j]->get_clear_value_vec(&output, &out_bitlen, &out_nvals);
      REQUIRE(out_bitlen == 32);
      REQUIRE(out_nvals == 1);
      REQUIRE_THAT(*(float*)output, WithinAbs(res[i][j], float_test_epsilon));
    }
  }

  delete party;
}

TEST_CASE("ABY small invert function is computed", "[aby_small_invert]") {
  int m = 7, n = 3;
  float one_d_M[] = {3, 0,
                     1,  // toy example
                     0, 0, 0, 0, 4, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0};
  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_invert(CLIENT, ctype, one_d_M, m, n);
    }
  });

  for (auto ctype : ctypes) {
    aby_invert(SERVER, ctype, one_d_M, m, n);
  }
  bob.join();
}

TEST_CASE("ABY large invert function is computed", "[aby_large_invert]") {
  int m = 12, n = 6;
  float one_d_M[] = {
      -9.1552734, 149.53613,  -59.509277, 18.310547, 0,         3.0517578,
      -164.79492, -64.086914, -47.302246, 0,         19.836426, 1.5258789,
      27.46582,   77.819824,  50.354004,  22.888184, 0,         1.5258789,
      -56.45752,  -42.724609, -53.405762, 0,         22.888184, 3.0517578,
      45.776367,  79.345703,  79.345703,  21.362305, 0,         0,
      -32.806396, -13.73291,  -20.599365, 0,         24.414062, 4.5776367,
      -21.362305, 108.3374,   -93.078613, 16.784668, 0,         3.0517578,
      -152.58789, -45.776367, -10.681152, 0,         19.836426, 3.0517578,
      6.1035156,  39.672852,  3.0517578,  21.362305, 0,         0,
      -39.672852, -21.362305, -16.784668, 0,         22.888184, 1.5258789,
      24.414062,  48.828125,  27.46582,   24.414062, 0,         0,
      -13.73291,  6.1035156,  12.207031,  0,         22.888184, 1.5258789};
  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_invert(CLIENT, ctype, one_d_M, m, n);
    }
  });

  for (auto ctype : ctypes) {
    aby_invert(SERVER, ctype, one_d_M, m, n);
  }
  bob.join();
}

void aby_twonorm(e_role role, e_sharing sharing, float* vect, int sz) {
  ABYParty* party = new ABYParty(role, address, port, slvl, bitlen, nthreads,
                                 mt_alg, reservegates, get_circuit_dir());
  std::vector<Sharing*>& sharings = party->GetSharings();
  Circuit* circ = sharings[sharing]->GetCircuitBuildRoutine();

  share** s_vect = new share*[sz];

  for (int i = 0; i < sz; ++i) {
    if (role == SERVER)
      s_vect[i] = circ->PutINGate((uint32_t*)&vect[i], 32, role);
    else
      s_vect[i] = circ->PutDummyINGate(32);
  }

  share* s_res = BuildTwoNormSqCircuit(s_vect, sz, (BooleanCircuit*)circ);

  share* s_out = circ->PutOUTGate(s_res, ALL);

  party->ExecCircuit();

  float cleartext_res = twonormsq(vect, sz);

  uint32_t* output;
  uint32_t out_bitlen, out_nvals;
  s_out->get_clear_value_vec(&output, &out_bitlen, &out_nvals);

  REQUIRE(out_bitlen == 32);
  REQUIRE(out_nvals == 1);
  REQUIRE_THAT(*(float*)output, WithinAbs(cleartext_res, float_test_epsilon));

  delete party;
  delete[] s_vect;
  delete s_out;
}

TEST_CASE("ABY two norm is computed", "[aby_twonorm]") {
  int sz = 5;
  float vec[] = {3.1, 5, 2, 4, 1};
  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_twonorm(CLIENT, ctype, vec, sz);
    }
  });

  for (auto ctype : ctypes) {
    aby_twonorm(SERVER, ctype, vec, sz);
  }
  bob.join();
}

void aby_localize(e_role role, e_sharing sharing, cv::Mat rvec, cv::Mat tvec,
                  cv::Mat cameraMatrix, cv::Mat distCoeffs,
                  vector<cv::Point3f> objectPoints,
                  vector<cv::Point2f> imagePoints, auto cleartext_localize_func,
                  auto secure_localize_func) {
  cv::Mat cvrvec = rvec.clone();
  cv::Mat cvtvec = tvec.clone();
  cv::solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs, cvrvec,
               cvtvec, false, cv::SOLVEPNP_ITERATIVE);
  if (role == SERVER) {
    cout << "opencv res: " << cvrvec.t() << cvtvec.t() << '\n';
  }

  float rt[6] = {rvec.at<float>(0), rvec.at<float>(1), rvec.at<float>(2),
                 tvec.at<float>(0), tvec.at<float>(1), tvec.at<float>(2)};
  float f = cameraMatrix.at<float>(0, 0);
  float cx = cameraMatrix.at<float>(0, 2);
  float cy = cameraMatrix.at<float>(1, 2);
  cleartext_localize_func(objectPoints, imagePoints, f, cx, cy, rt, "");
  if (role == SERVER) {
    printVector("cleartext res", rt, 6);
  }

  // check cleartext vs opencv
  for (int i = 0; i < 3; i++) {
    REQUIRE_THAT(cvrvec.at<float>(i),
                 WithinRel(rt[i], localization_tol_rel) ||
                     WithinAbs(rt[i], localization_tol_abs));
    REQUIRE_THAT(cvtvec.at<float>(i),
                 WithinRel(rt[i + 3], localization_tol_rel) ||
                     WithinAbs(rt[i + 3], localization_tol_abs));
  }

  vector<float> res(6, 0);
  aby_localize_wrapper(role, sharing, rvec, tvec, cameraMatrix, distCoeffs,
                       objectPoints, imagePoints, res, secure_localize_func);

  if (role == SERVER) {
    cout << "secure result: ";
    for (auto const& f : res)
      cout << f << " ";
    cout << '\n';
  }
  for (int i = 0; i < 6; i++) {
    REQUIRE_THAT(res[i], WithinRel(rt[i], localization_tol_rel) ||
                             WithinAbs(rt[i], localization_tol_abs));
  }
}

//#if PPL_FLOW != PPL_FLOW_SiSL
TEST_CASE("ABY Gauss Newton pose estimation is computed", "[aby_gn]") {
  float f = 715;
  float cx = 354;
  float cy = 245;
  float _cM[] = {f, 0, cx, 0, f, cy, 0, 0, 1};
  cv::Mat cameraMatrix = cv::Mat(3, 3, cv::DataType<float>::type, _cM);
  cv::Mat distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
  cv::Mat rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::Mat tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0
  vector<cv::Point2f> imagePoints;
  vector<cv::Point3f> objectPoints;
  Hoffs2DPoints(imagePoints);
  Hoffs3DPoints(objectPoints);

  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_localize(CLIENT, ctype, rvec, tvec, cameraMatrix, distCoeffs,
                   objectPoints, imagePoints, gaussNewton<float>,
                   BuildAndRunGaussNewton);
    }
  });

  for (auto ctype : ctypes) {
    aby_localize(SERVER, ctype, rvec, tvec, cameraMatrix, distCoeffs,
                 objectPoints, imagePoints, gaussNewton<float>,
                 BuildAndRunGaussNewton);
  }
  bob.join();
}
//#endif

// Note this test can take around an hour
TEST_CASE("ABY Levenburg Marquardt pose estimation is computed", "[aby_lm]") {
  float f = 715;
  float cx = 354;
  float cy = 245;
  float _cM[] = {f, 0, cx, 0, f, cy, 0, 0, 1};
  cv::Mat cameraMatrix = cv::Mat(3, 3, cv::DataType<float>::type, _cM);
  cv::Mat distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
  cv::Mat rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::Mat tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0
  vector<cv::Point2f> imagePoints;
  vector<cv::Point3f> objectPoints;
  Hoffs2DPoints(imagePoints);
  Hoffs3DPoints(objectPoints);

  std::thread bob([&]() {
    for (auto ctype : ctypes) {
      aby_localize(CLIENT, ctype, rvec, tvec, cameraMatrix, distCoeffs,
                   objectPoints, imagePoints, lm<float>, BuildAndRunLM);
    }
  });

  for (auto ctype : ctypes) {
    aby_localize(SERVER, ctype, rvec, tvec, cameraMatrix, distCoeffs,
                 objectPoints, imagePoints, lm<float>, BuildAndRunLM);
  }
  bob.join();
}