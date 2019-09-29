//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <stdio.h>
#include <iostream>
#include <string>
#include <thread>

#include <jlog.h>
#include <april_snail_features.hpp>
#include <hoff_features.hpp>
#include <test_utils.hpp>

#include <gaussnewtonlocalization.h>
#include <lmlocalization.h>
#include <cleartext-ref/gaussnewtonlocalization.hpp>
#include <cleartext-ref/lmlocalization.hpp>

#include "emp-sh2pc/emp-sh2pc.h"

#include <math.h>

using namespace emp;
using Catch::Matchers::WithinAbs;
using std::cout;
using std::vector;

constexpr const int port = 8080;

TEST_CASE("EMP trig functions are computed", "[emp_trig]") {
  float angle = .2;

  std::thread bob([angle]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_trig(io, BOB, angle);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_trig(io, ALICE, angle);
  bob.join();
  delete io;
}

TEST_CASE("EMP matrix multiply is computed", "[emp_mat]") {
  cv::Mat A = cv::Mat::ones(3, 5, cv::DataType<float>::type);
  cv::Mat B = cv::Mat::ones(5, 4, cv::DataType<float>::type) * 2;

  std::thread bob([A, B]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_mat(io, BOB, A, B);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_mat(io, ALICE, A, B);
  bob.join();
  delete io;
}

TEST_CASE("EMP rodriguez transformation is computed", "[emp_rod]") {
  float r[3] = {.2, .4, .2};  // 3x1

  std::thread bob([&r]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_rod(io, BOB, r);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_rod(io, ALICE, r);
  bob.join();
  delete io;
}

TEST_CASE("EMP point projection is computed", "[emp_proj]") {
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
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_proj(io, BOB, cameraMatrix, distCoeffs, objectPoints, imagePoints);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_proj(io, ALICE, cameraMatrix, distCoeffs, objectPoints, imagePoints);
  bob.join();
  delete io;
}

TEST_CASE("EMP SVD sign function is computed", "[emp_svd]") {
  float a = 2.40, b = 0;

  std::thread bob([a, b]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_svd_sign(io, BOB, a, b);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_svd_sign(io, ALICE, a, b);
  bob.join();
  delete io;
}

TEST_CASE("EMP SVD pythag function is computed", "[emp_svd]") {
  float a = -2.40, b = 96.0;

  std::thread bob([a, b]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_svd_pythag(io, BOB, a, b);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_svd_pythag(io, ALICE, a, b);
  bob.join();
  delete io;
}

// TEST_CASE("Small SVD function is computed", "[svd]") {
//   int m=5,n=4;
//   float** in = new float*[m];
//   for(int i=0; i<m; i++) {
//     in[i] = new float[n]();
//     for (int j=0; j<n; j++) {
//       in[i][j]=i+j;
//     }
//   }
//   in[0][0] = 0;
//   in[1][0] = 0;
//   in[2][0] = 0;
//   in[3][0] = 0;
//   in[4][0] = 0;
//
//   std::thread bob([in, m, n]() {
//     NetIO *io = new NetIO("127.0.0.1", port);
//     emp_svd(io, BOB, in, m, n);
//     delete io;
//   });
//
//   NetIO *io = new NetIO(nullptr, port);
//   emp_svd(io, ALICE, in, m, n);
//   bob.join();
//   for (int i=0; i<m; i++) {
//     delete[] in[i];
//   }
//   delete[] in;
//   delete io;
// }

TEST_CASE("EMP large SVD function is computed", "[emp_svd]") {
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

  std::thread bob([in, m, n]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_svd(io, BOB, in, m, n);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_svd(io, ALICE, in, m, n);
  bob.join();
  for (int i = 0; i < m; i++) {
    delete[] in[i];
  }
  delete[] in;
  delete io;
}

TEST_CASE("EMP small invert function is computed", "[emp_small_invert]") {
  int m = 7, n = 3;
  float one_d_M[] = {3, 0,
                     1,  // toy example
                     0, 0, 0, 0, 4, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0};

  std::thread bob([&one_d_M, m, n]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_invert(io, BOB, one_d_M, m, n);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_invert(io, ALICE, one_d_M, m, n);
  bob.join();
  delete io;
}

TEST_CASE("EMP large invert function is computed", "[emp_large_invert]") {
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

  std::thread bob([&one_d_M, m, n]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_invert(io, BOB, one_d_M, m, n);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_invert(io, ALICE, one_d_M, m, n);
  bob.join();
  delete io;
}

TEST_CASE("EMP two norm is computed", "[emp_twonorm]") {
  int sz = 5;
  float vec[] = {3.1, 5, 2, 4, 1};

  std::thread bob([&]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_twonorm(io, BOB, vec, sz);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_twonorm(io, ALICE, vec, sz);
  bob.join();
  delete io;
}

// SiSL is lm only
//#if PPL_FLOW != PPL_FLOW_SiSL
TEST_CASE("EMP Gauss Newton pose estimation is computed", "[emp_gn]") {
  float f = 715;
  float cx = 354;
  float cy = 245;
  float _cM[] = {f, 0, cx, 0, f, cy, 0, 0, 1};
  cv::Mat cameraMatrix = cv::Mat(3, 3, cv::DataType<float>::type, _cM);
  cv::Mat distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
  cv::Mat rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::Mat tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  vector<cv::Point2f> imagePoints;
  vector<cv::Point3f> objectPoints;
  Hoffs2DPoints(imagePoints);
  Hoffs3DPoints(objectPoints);

  std::thread bob([&]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_localize(io, BOB, rvec, tvec, cameraMatrix, distCoeffs, objectPoints,
                 imagePoints, gaussNewton<float>, BuildGaussNewton);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_localize(io, ALICE, rvec, tvec, cameraMatrix, distCoeffs, objectPoints,
               imagePoints, gaussNewton<float>, BuildGaussNewton);
  bob.join();
  delete io;
}
//#endif

TEST_CASE("EMP Levenburg Marquardt pose estimation is computed on Hoff points",
          "[emp_lm]") {
  float f = 715;
  float cx = 354;
  float cy = 245;
  float _cM[] = {f, 0, cx, 0, f, cy, 0, 0, 1};
  cv::Mat cameraMatrix = cv::Mat(3, 3, cv::DataType<float>::type, _cM);
  cv::Mat distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
  cv::Mat rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::Mat tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  vector<cv::Point2f> imagePoints;
  vector<cv::Point3f> objectPoints;
  Hoffs2DPoints(imagePoints);
  Hoffs3DPoints(objectPoints);

  std::thread bob([&]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_localize(io, BOB, rvec, tvec, cameraMatrix, distCoeffs, objectPoints,
                 imagePoints, lm<float>, BuildLM);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_localize(io, ALICE, rvec, tvec, cameraMatrix, distCoeffs, objectPoints,
               imagePoints, lm<float>, BuildLM);
  bob.join();
  delete io;
}

TEST_CASE(
    "EMP Levenburg Marquardt pose estimation is computed on April tag points",
    "[emp_lm_april]") {
  float f = 715;
  float cx = 354;
  float cy = 245;
  float _cM[] = {f, 0, cx, 0, f, cy, 0, 0, 1};
  cv::Mat cameraMatrix = cv::Mat(3, 3, cv::DataType<float>::type, _cM);
  cv::Mat distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
  cv::Mat rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::Mat tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  tvec.at<float>(2) = 1.0f;  // cleartext has bug where it fails if z=0
  vector<cv::Point2f> imagePoints;
  vector<cv::Point3f> objectPoints;
  AprilSnail3DPoints(objectPoints);
  AprilSnail2DPoints(imagePoints);

  std::thread bob([&]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    emp_localize(io, BOB, rvec, tvec, cameraMatrix, distCoeffs, objectPoints,
                 imagePoints, lm<float>, BuildLM);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  emp_localize(io, ALICE, rvec, tvec, cameraMatrix, distCoeffs, objectPoints,
               imagePoints, lm<float>, BuildLM);
  bob.join();
  delete io;
}
