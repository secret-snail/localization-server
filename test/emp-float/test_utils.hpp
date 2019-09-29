//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <math.h>
#include <stdio.h>

#include <iostream>
#include <string>
#include <thread>

#include <jlog.h>
#include <test_params.h>
#include <localize_wrapper.hpp>

#include <gaussnewtonlocalization.h>
#include <invert.h>
#include <lmlocalization.h>
#include <matmult.h>
#include <printutil.h>
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

#include "emp-sh2pc/emp-sh2pc.h"

using namespace emp;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using std::cout;
using std::vector;

void emp_trig(NetIO* io, int party, float angle) {
  setup_semi_honest(io, party);
  Float v = Float(angle, ALICE);
  Float sinres = BuildSinCircuit(v);
  REQUIRE_THAT(sinres.reveal<double>(),
               WithinAbs(sin(angle), float_test_epsilon));

  Float cosres = BuildCosCircuit(v);
  REQUIRE_THAT(cosres.reveal<double>(),
               WithinAbs(cos(angle), float_test_epsilon));
  finalize_semi_honest();
}

void emp_mat(NetIO* io, int party, cv::Mat A, cv::Mat B) {
  setup_semi_honest(io, party);
  cv::Mat clear_res = cv::Mat::zeros(A.rows, B.cols, cv::DataType<float>::type);

  matmult(A.ptr<float>(), A.rows, A.cols, B.ptr<float>(), B.rows, B.cols,
          clear_res.ptr<float>());

  void* raw_memoryA = operator new[](A.rows* A.cols * sizeof(Float));
  void* raw_memoryB = operator new[](B.rows* B.cols * sizeof(Float));
  void* raw_memoryres = operator new[](A.rows* B.cols * sizeof(Float));
  Float* eA = static_cast<Float*>(raw_memoryA);
  Float* eB = static_cast<Float*>(raw_memoryB);
  Float* eres = static_cast<Float*>(raw_memoryres);

  for (int i = 0; i < A.rows * A.cols; ++i) {
    new (&eA[i]) Float(A.ptr<float>()[i], ALICE);
  }

  for (int i = 0; i < B.rows * B.cols; ++i) {
    new (&eB[i]) Float(B.ptr<float>()[i], ALICE);
  }

  BuildMatmultCircuit(eA, A.rows, A.cols, eB, B.rows, B.cols, eres);

  for (int i = 0; i < A.rows * B.cols; ++i) {
    REQUIRE_THAT(eres[i].reveal<double>(),
                 WithinAbs(clear_res.ptr<float>()[i], float_test_epsilon));
  }

  finalize_semi_honest();
  delete[] eA;
  delete[] eB;
  delete[] eres;
}

void emp_rod(NetIO* io, int party, float* r) {
  setup_semi_honest(io, party);

  cv::Mat R = cv::Mat::zeros(3, 4, cv::DataType<float>::type);
  rodrigues(r, R.ptr<float>());

  void* raw_memorysr = operator new[](3 * sizeof(Float));
  Float* sr = static_cast<Float*>(raw_memorysr);
  new (&sr[0]) Float(r[0], ALICE);
  new (&sr[1]) Float(r[1], ALICE);
  new (&sr[2]) Float(r[2], ALICE);

  void* raw_memoryR = operator new[](12 * sizeof(Float));
  Float* sR = static_cast<Float*>(raw_memoryR);

  BuildRodriguesCircuit(sr, sR);

  cout << "rodrigues result:" << endl;
  for (int i = 0; i < 12; ++i) {
    if (i == 3 || i == 7 || i == 11) {
      continue;
    }
    REQUIRE_THAT(sR[i].reveal<double>(),
                 WithinAbs(R.ptr<float>()[i], float_test_epsilon));
  }
  delete[] sr;
  delete[] sR;

  finalize_semi_honest();
}

// From points, first estimate the pose using opencv.
// Then reproject into 3d using the estimated pose
// and compare to the origional points.
void emp_proj(NetIO* io, int party, cv::Mat cameraMatrix, cv::Mat distCoeffs,
              vector<cv::Point3f> objectPoints,
              vector<cv::Point2f> imagePoints) {
  cv::Mat rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::Mat tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);

  // Find rotation and translation
  cv::solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs, rvec, tvec,
               false, cv::SOLVEPNP_ITERATIVE);

  int nOp = objectPoints.size();
  cv::Mat x = cv::Mat::zeros(6, 1, cv::DataType<float>::type);
  float* xptr = x.ptr<float>();
  xptr[0] = rvec.at<float>(0);
  xptr[1] = rvec.at<float>(1);
  xptr[2] = rvec.at<float>(2);
  xptr[3] = tvec.at<float>(0);
  xptr[4] = tvec.at<float>(1);
  xptr[5] = tvec.at<float>(2);

  cv::Mat matOP = cv::Mat(objectPoints, false);
  matOP = matOP.reshape(1, matOP.total());
  matOP.convertTo(matOP, cv::DataType<float>::type);
  cv::Mat P = cv::Mat::ones(4, nOp, cv::DataType<float>::type);
  P(cv::Range(0, 3), cv::Range(0, nOp)) =
      matOP.t();  // copy points into new matrix

  // cleartext projection
  cv::Mat projected = cv::Mat(3, nOp, cv::DataType<float>::type);
  projectPoints(P.ptr<float>(), x.ptr<float>(), cameraMatrix.ptr<float>(),
                projected.ptr<float>(), nOp);
  for (int i = 0; i < nOp; ++i) {
    REQUIRE_THAT(imagePoints[i].x,
                 WithinAbs(projected.at<float>(0, i), reprojection_tol));
    REQUIRE_THAT(imagePoints[i].y,
                 WithinAbs(projected.at<float>(1, i), reprojection_tol));
  }

  // mpc projection
  setup_semi_honest(io, party);
  void* raw_memoryP = operator new[](4 * nOp * sizeof(Float));
  void* raw_memoryx = operator new[](6 * sizeof(Float));
  void* raw_memoryK = operator new[](9 * sizeof(Float));
  void* raw_memoryres = operator new[](3 * nOp * sizeof(Float));
  Float* sP = static_cast<Float*>(raw_memoryP);
  Float* sx = static_cast<Float*>(raw_memoryx);
  Float* sK = static_cast<Float*>(raw_memoryK);
  Float* sres = static_cast<Float*>(raw_memoryres);

  for (int i = 0; i < 3 * nOp; i++) {
    new (&sP[i]) Float(P.ptr<float>()[i], ALICE);
  }
  for (int i = 3 * nOp; i < 4 * nOp; i++) {  // need contant 1s
    new (&sP[i]) Float(1.0, PUBLIC);
  }
  for (int i = 0; i < 6; i++) {
    new (&sx[i]) Float(x.ptr<float>()[i], ALICE);
  }
  for (int i = 0; i < 9; i++) {
    new (&sK[i]) Float(cameraMatrix.ptr<float>()[i], ALICE);
  }

  BuildProjectPointsCircuit(sP, sx, sK, sres, nOp, true);

  for (int i = 0; i < 2 * nOp; ++i) {
    REQUIRE_THAT(projected.ptr<float>()[i],
                 WithinAbs(sres[i].reveal<double>(), reprojection_tol));
  }
  delete[] sP;
  delete[] sx;
  delete[] sK;
  delete[] sres;
  finalize_semi_honest();
}

void emp_svd_sign(NetIO* io, int party, float a, float b) {
  setup_semi_honest(io, party);
  float cleartext_res = mysign(a, b);
  Float sa = Float(a, ALICE);
  Float sb = Float(b, ALICE);
  BuildSignCircuit(&sa, &sb);
  REQUIRE_THAT(sa.reveal<double>(),
               WithinAbs(cleartext_res, float_test_epsilon));
  finalize_semi_honest();
}

void emp_svd_pythag(NetIO* io, int party, float a, float b) {
  setup_semi_honest(io, party);
  float cleartext_res = mypythag(a, b);
  Float sa = Float(a, ALICE);
  Float sb = Float(b, ALICE);
  Float out = BuildPythagCircuit(&sa, &sb);
  REQUIRE_THAT(out.reveal<double>(),
               WithinAbs(cleartext_res, float_test_epsilon));
  finalize_semi_honest();
}

void emp_svd(NetIO* io, int party, float** in, int m, int n) {
  // linearize for cv
  float* cvin = new float[m * n];
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      cvin[i * n + j] = in[i][j];
    }
  }
  cv::Mat cvinM = cv::Mat(m, n, cv::DataType<float>::type, in);
  cv::SVD cvsvd(cvinM, cv::SVD::FULL_UV);  // constructor
  if (party == ALICE) {
    cout << "opencv result\n"
         << cvsvd.u << '\n'
         << cvsvd.w << '\n'
         << cvsvd.vt << '\n';
  }
  delete[] cvin;

  // svd overwrites with u matrix, so copy
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

  setup_semi_honest(io, party);
  Float** sa = new Float*[m];
  for (int i = 0; i < m; i++) {
    sa[i] = static_cast<Float*>(operator new[](n * sizeof(Float)));
    for (int j = 0; j < n; j++) {
      new (&sa[i][j]) Float(in[i][j], ALICE);
    }
  }
  Float* sw = static_cast<Float*>(operator new[](n * sizeof(Float)));
  for (int i = 0; i < n; i++) {
    sw[i] = Float(0.0, PUBLIC);
  }
  Float** sv = new Float*[n];
  for (int i = 0; i < n; i++) {
    sv[i] = static_cast<Float*>(operator new[](n * sizeof(Float)));
    for (int j = 0; j < n; j++) {
      sv[i][j] = Float(0.0, PUBLIC);
    }
  }

  svdcmp(a, m, n, w, v);
  if (party == ALICE) {
    cout << "cleartext" << '\n';
    printMatrix("a", a, m, n);
    printVector("w", w, n);
    printMatrix("v", v, n, n);
  }

  BuildSvdCircuit(sa, m, n, sw, sv);
  // both parties must print if uncommented
  cout << "secure" << '\n';
  printFloatMatrix(sa, m, n, party == BOB);
  printFloatVector(sw, n, party == BOB);
  printFloatMatrix(sv, n, n, party == BOB);

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      // svd is not unique, do not compare with opencv
      // REQUIRE_THAT(a[i][j], WithinAbs(cvsvd.u.at<float>(i,j), svd_tol));
      REQUIRE_THAT(sa[i][j].reveal<double>(), WithinAbs(a[i][j], svd_tol));
    }
    delete[] a[i];
    delete[] sa[i];
  }
  delete[] a;
  delete[] sa;

  for (int i = 0; i < n; i++) {
    // svd is not unique, do not compare with opencv
    // REQUIRE_THAT(w[i], WithinAbs(cvsvd.w.at<float>(i), svd_tol));
    REQUIRE_THAT(sw[i].reveal<double>(), WithinAbs(w[i], svd_tol));
  }
  delete[] w;
  delete[] sw;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // svd is not unique, do not compare with opencv
      // REQUIRE_THAT(v[i][j], WithinAbs(cvsvd.vt.at<float>(i,j), svd_tol));
      REQUIRE_THAT(sv[i][j].reveal<double>(), WithinAbs(v[i][j], svd_tol));
    }
    delete[] v[i];
    delete[] sv[i];
  }
  delete[] v;
  delete[] sv;

  finalize_semi_honest();
}

void emp_invert(NetIO* io, int party, float* in, int m, int n) {
  setup_semi_honest(io, party);

  // opencv computation
  cv::Mat cvM = cv::Mat(m, n, cv::DataType<float>::type, in);
  cv::Mat cvres;
  invert(cvM, cvres, cv::DECOMP_SVD);
  // cout << "opencv result:\n" << ores << endl;

  // convert to 2d matrix
  float** M = new float*[m];
  for (int i = 0; i < m; i++) {
    M[i] = new float[n];
    for (int j = 0; j < n; j++) {
      M[i][j] = in[i * n + j];
    }
  }

  // result is nxm (not mxn)
  float** res = new float*[n];
  for (int i = 0; i < n; i++) {
    res[i] = new float[m];
  }

  // share arrays
  Float** a = new Float*[m];
  for (int i = 0; i < m; i++) {
    a[i] = static_cast<Float*>(operator new[](n * sizeof(Float)));
    for (int j = 0; j < n; j++) {
      new (&a[i][j]) Float(M[i][j], ALICE);
    }
  }

  Float** sres = new Float*[n];
  for (int i = 0; i < n; i++) {
    sres[i] = static_cast<Float*>(operator new[](m * sizeof(Float)));
    // for(int j=0; j<nRows; j++) {
    //   sres[i][j] = Float(0.0, PUBLIC);
    // }
  }

  myinvert(M, m, n, res);
  // printMatrix("cleartext result", res, n, m);

  // check cleartext vs opencv
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      REQUIRE_THAT(cvres.at<float>(i, j),
                   WithinAbs(res[i][j], float_test_epsilon));
    }
  }

  BuildInvertCircuit(a, m, n, sres);

  // check privacy preserving vs cleartext
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      REQUIRE_THAT(sres[i][j].reveal<double>(),
                   WithinAbs(res[i][j], float_test_epsilon));
    }
    delete[] res[i];
    delete[] sres[i];
  }
  delete[] res;
  delete[] sres;

  for (int i = 0; i < m; i++) {
    delete[] M[i];
    delete[] a[i];
  }
  delete[] M;
  delete[] a;

  finalize_semi_honest();
}

void emp_twonorm(NetIO* io, int party, float* vec, int sz) {
  setup_semi_honest(io, party);
  Float* svec = static_cast<Float*>(operator new[](sz * sizeof(Float)));
  for (int i = 0; i < sz; i++) {
    svec[i] = Float(vec[i], PUBLIC);
  }

  Float sres = BuildTwoNormSqCircuit(svec, sz);

  REQUIRE_THAT(sres.reveal<double>(),
               WithinAbs(twonormsq(vec, sz), float_test_epsilon));
  delete[] svec;
  finalize_semi_honest();
}

void emp_localize(NetIO* io, int party, cv::Mat rvec, cv::Mat tvec,
                  cv::Mat cameraMatrix, cv::Mat distCoeffs,
                  vector<cv::Point3f> objectPoints,
                  vector<cv::Point2f> imagePoints, auto cleartext_localize_func,
                  auto secure_localize_func) {
  setup_semi_honest(io, party);

  cv::Mat cvrvec = rvec.clone();
  cv::Mat cvtvec = tvec.clone();
  cv::solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs, cvrvec,
               cvtvec, false, cv::SOLVEPNP_ITERATIVE);
  if (party == ALICE) {
    cout << "opencv\n" << cvrvec.t() << cvtvec.t() << '\n';
  }

  float rt[6] = {rvec.at<float>(0), rvec.at<float>(1), rvec.at<float>(2),
                 tvec.at<float>(0), tvec.at<float>(1), tvec.at<float>(2)};
  float f = cameraMatrix.at<float>(0, 0);
  float cx = cameraMatrix.at<float>(0, 2);
  float cy = cameraMatrix.at<float>(1, 2);
  cleartext_localize_func(objectPoints, imagePoints, f, cx, cy, rt, "");
  if (party == ALICE) {
    cout << "cleartext\n";
    printVector("rt", rt, 6);
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

  vector<float> res(6);
  emp_localize_wrapper(party, rvec, tvec, cameraMatrix, distCoeffs,
                       objectPoints, imagePoints, res, secure_localize_func);

  // check opencv vs secure
  for (int i = 0; i < 6; i++) {
    if (party == BOB) {
      cout << "secure result: ";
      for (auto const& f : res)
        cout << f << ' ';
      cout << '\n';
    }
    REQUIRE_THAT(res[i], WithinRel(rt[i], localization_tol_rel) ||
                             WithinAbs(rt[i], localization_tol_abs));
  }

  finalize_semi_honest();
}
