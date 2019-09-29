#pragma once
//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <jlog.h>
#include <printutil.h>
#include <privacyconf.h>
#include <stdio.h>
#include <iostream>
#include <typeinfo>

#include <cleartext-ref/invert.hpp>
#include <cleartext-ref/matmult.hpp>
#include <cleartext-ref/projectpoints.hpp>
#include <cleartext-ref/svd.hpp>
#include <cleartext-ref/twonormsq.hpp>

using namespace std;

// See William Hoff's lecure: https://www.youtube.com/watch?v=kq3c6QpcAGc
template <typename T>
bool gaussNewton(vector<cv::Point3_<T>> threeDPts,
                 vector<cv::Point_<T>> twoDPts, T _f, T _cx, T _cy,
                 T* _x, /* initial guess for { r1, r2, r3, t1, t2, t3 } */
                 string logprint = "") {

  if (typeid(T) == typeid(int) || typeid(T) == typeid(int16_t) ||
      typeid(T) == typeid(uint16_t) || typeid(T) == typeid(int32_t) ||
      typeid(T) == typeid(uint32_t) || typeid(T) == typeid(int64_t) ||
      typeid(T) == typeid(uint64_t)) {
    cout << "type not supported\n";
    assert(false);
  }

  if (_verbosity & DBG_FLOW)
    cout << "Gauss-Newton Gradient Decent" << endl;
  int n = threeDPts.size();

  T min_er;
  T jacob_epsilon;

  // Initial Guess
  T x[6];
  // Camera Params
  T K[9] = {_f, 0, _cx, 0, _f, _cy, 0, 0, 1};
  // 3D World Points
  //  _             _
  // |  x1 x2 x3     |
  // |  y1 y2 y3 ... |
  // |  z1 z2 z3     |
  // |_ 1  1  1     _|
  T P[4 * n];
  // 2D Image Points
  //  _             _
  // |  u1 u2 u3     |
  // |  v1 v2 v3 ... |
  // |_ 1  1  1     _|
  T y0Homog[3 * n];

  min_er = MIN_ER;
  jacob_epsilon = JACOB_EPSILON;
  for (int i = 0; i < 6; i++) {
    x[i] = _x[i];
  }
  for (int i = 0; i < n; i++) {
    P[i] = threeDPts[i].x;
    P[n + i] = threeDPts[i].y;
    P[(2 * n) + i] = threeDPts[i].z;
    P[(3 * n) + i] = 1;

    y0Homog[i] = twoDPts[i].x;
    y0Homog[n + i] = twoDPts[i].y;
    y0Homog[(2 * n) + i] = 1;
  }

  // throw away last "row" of y0Homog (constant ones),
  // interleave, and transpose into 2nx1 vector
  // e.g. [x1; y1; x2; y2 ...]
  T y0[2 * n];
  for (int p = 0; p < n; ++p) {
    y0[p * 2] = y0Homog[p];
    y0[(p * 2) + 1] = y0Homog[p + n];
  }

  if (_verbosity & DBG_ARGS) {
    printVector("P", P, 4 * n);
    printVector("y0Homog", y0Homog, 3 * n);
    printVector("y0", y0, 2 * n);
  }

  // GN Iteration
  int i;
  for (i = 0; i < GN_MAX_ITR; i++) {
    if (_verbosity & DBG_FLOW) {
      cout << RED << "\n\nGN Iteration " << i << RESET << endl;
      printVector("x", x, 6);
    }

    // project points using x guess
    T yHomog[3 * n];
    projectPoints<T>(P, x, K, yHomog, n);
    // throw away last "row" of result (constant ones),
    // interleave, and transpose y into 2nx1 vector
    // e.g. [x1; y1; x2; y2 ...]
    T y[2 * n];
    for (int p = 0; p < n; ++p) {
      y[p * 2] = yHomog[p];
      y[(p * 2) + 1] = yHomog[p + n];
    }
    if (_verbosity & DBG_PROJECT) {
      printVector("y", y, 2 * n);
    }

    // calculate jacobian
    //  -                                      -
    // | df/dr1 df/dr2 df/dr3 df/t1 df/t2 df/t3 |
    // |                .                       |
    // |                .                       |
    // |                . (2*n)                 |
    //  -                                      -
    T** jacob = new T*[n * 2];
    for (int p = 0; p < n * 2; p++)
      jacob[p] = new T[6];
    for (int j = 0; j < 6; j++) {  // each dof
      // perturb x
      T oldx = x[j];
      x[j] += jacob_epsilon;

      // project with epsilon
      T ytempHomog[3 * n];
      projectPoints<T>(P, x, K, ytempHomog, n);
      // throw away last "row" of result,
      // reshape, and transpose y into 2nx1 vector
      // e.g. [x1; y1; x2; y2 ...]
      T ytemp[2 * n];
      for (int p = 0; p < n; ++p) {
        ytemp[p * 2] = ytempHomog[p];
        ytemp[(p * 2) + 1] = ytempHomog[p + n];
      }
      if (_verbosity & DBG_PROJECT) {
        printVector("ytemp", ytemp, 2 * n);
      }

      // copy into jacob
      for (int k = 0; k < 2 * n; k++) {
        jacob[k][j] = (ytemp[k] - y[k]) / jacob_epsilon;
      }

      // put x back
      x[j] = oldx;
    }
    if (_verbosity & DBG_JACOB) {
      printMatrix("jacob", jacob, 2 * n, 6);
    }

    // error (nx2)
    T dy[2 * n];
    for (int p = 0; p < 2 * n; p++) {
      dy[p] = y0[p] - y[p];
    }
    if (_verbosity & DBG_ER)
      printVector("dy", dy, 2 * n);

    T** jacobI = new T*[6];  // 6x2n
    for (int p = 0; p < 6; p++)
      jacobI[p] = new T[2 * n];
    // pseudo inverse to solve dy = J dx -> dx = Jdag dy
    // OpenCV's invert function
    // cv::Mat jacobI;
    // invert(jacob, jacobI, cv::DECOMP_SVD);

    //// Compute SVD directly with type T
    // if (myinvert<T>(jacob, 2*n, 6, jacobI)) {
    //   for (int p=0; p<2*n; p++) {
    //       delete[] jacob[p];
    //   }
    //   delete[] jacob;
    //   for (int p=0; p<6; p++) {
    //       delete[] jacobI[p];
    //   }
    //   delete[] jacobI;
    //   return false;
    // }

    // Force compute SVD with floats
    // Convert to floats
    float** f_jacob = new float*[2 * n];  // 2nx6
    for (int p = 0; p < 2 * n; p++) {
      f_jacob[p] = new float[6];
      for (int pp = 0; pp < 6; pp++) {
        f_jacob[p][pp] = jacob[p][pp];
      }
    }
    float** f_jacobI = new float*[6];  // 6x2n
    for (int p = 0; p < 6; p++)
      f_jacobI[p] = new float[2 * n];
    if (myinvert<float>(f_jacob, 2 * n, 6, f_jacobI)) {
      for (int p = 0; p < 2 * n; p++) {
        delete[] jacob[p];
        delete[] f_jacob[p];
      }
      delete[] jacob;
      delete[] f_jacob;
      for (int p = 0; p < 6; p++) {
        delete[] jacobI[p];
        delete[] f_jacobI[p];
      }
      delete[] jacobI;
      delete[] f_jacobI;
      return false;
    }
    // Convert back to T template type
    for (int p = 0; p < 6; p++) {
      for (int pp = 0; pp < 2 * n; pp++) {
        jacobI[p][pp] = f_jacobI[p][pp];
      }
    }
    for (int p = 0; p < 6; p++)
      delete[] f_jacobI[p];
    delete[] f_jacobI;
    for (int p = 0; p < 2 * n; p++)
      delete[] f_jacob[p];
    delete[] f_jacob;

    // cout << "jacobI" << endl << jacobI << endl;

    for (int p = 0; p < 2 * n; p++) {
      delete[] jacob[p];
    }
    delete[] jacob;

    // dx = jacobI * dy;
    // linearize jacobI for matmult
    T jacobILinear[6 * 2 * n];
    for (int p = 0; p < 6; p++) {
      for (int pp = 0; pp < 2 * n; pp++) {
        jacobILinear[(p * 2 * n) + pp] = jacobI[p][pp];
      }
    }
    for (int p = 0; p < 6; p++) {
      delete[] jacobI[p];
    }
    delete[] jacobI;
    T dx[6];
    matmult<T>(jacobILinear, 6, 2 * n, dy, 2 * n, 1, dx);
    if (_verbosity & DBG_POSE_UPDATE)
      printVector("dx", dx, 6);

    // break if error under threshold
    // if (twonormsq(dy, 2 * n) < min_er)
    // break;
    if (twonormsq(dx, 6) < min_er)
      break;

    // update pose
    for (int p = 0; p < 6; p++)
      x[p] = x[p] + dx[p];
  }
  if (_verbosity & DBG_FLOW)
    cout << "Did " << i << " GN iterations\n";
  printf("SeNtInAl,xy,%s,gn_%s_iterations_vs_numpts,%d,%d\n", __FUNCTION__,
         logprint.c_str(), n, i);

  // copy data back to res
  for (int i = 0; i < 6; i++) {
    _x[i] = x[i];
  }

  return i < GN_MAX_ITR;
}
