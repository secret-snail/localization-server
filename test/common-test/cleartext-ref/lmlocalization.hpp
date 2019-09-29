#pragma once
//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <jlog.h>
#include <privacyconf.h>
#include <stdio.h>
#include <iostream>
#include <typeinfo>
#include <vector>

#include <cleartext-ref/invert.hpp>
#include <cleartext-ref/matmult.hpp>
#include <cleartext-ref/projectpoints.hpp>
#include <cleartext-ref/svd.hpp>
#include <cleartext-ref/twonormsq.hpp>

using namespace std;

template <typename T>
bool lm(vector<cv::Point3_<T>> threeDPts, vector<cv::Point_<T>> twoDPts, T _f,
        T _cx, T _cy, T* _x, /* initial guess for { r1, r2, r3, t1, t2, t3 } */
        string logprint = "") {

  if (typeid(T) == typeid(int) || typeid(T) == typeid(int16_t) ||
      typeid(T) == typeid(uint16_t) || typeid(T) == typeid(int32_t) ||
      typeid(T) == typeid(uint32_t) || typeid(T) == typeid(int64_t) ||
      typeid(T) == typeid(uint64_t)) {
    cout << "type not supported\n";
    assert(false);
  }

  if (_verbosity & DBG_FLOW)
    cout << "Levenberg–Marquardt\n";
  int n = threeDPts.size();

  T min_er;
  T jacob_epsilon;
  T lambda_init;
  T lambda_max;
  T lambda_min;

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
  lambda_init = LM_LAMBDA_INIT;
  lambda_max = LM_LAMBDA_MAX;
  lambda_min = LM_LAMBDA_MIN;
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

  T lambda = lambda_init;
  T prevErrNorm = std::numeric_limits<T>::max();

  // LM Iteration
  int i;
  for (i = 0; i < LM_MAX_ITR; i++) {
    if (_verbosity & DBG_FLOW) {
      cout << RED << "\n\nLM Iteration " << i << RESET << endl;
      printVector("x", x, 6);
    }

    // project points using x guess
    T yHomog[3 * n];
    projectPoints(P, x, K, yHomog, n);
    // throw away last "row" of result (constant ones),
    // interleave, and transpose y into 2nx1 vector
    // e.g. [x1; y1; x2; y2 ...]
    T y[2 * n];
    for (int p = 0; p < n; ++p) {
      y[p * 2] = yHomog[p];
      y[(p * 2) + 1] = yHomog[p + n];
    }
    if (_verbosity & DBG_PROJECT) {
      printVector("reshapen y", y, 2 * n);
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
      for (int k = 0; k < 2 * n; k++) {  // calculate ^ by hand
        jacob[k][j] = (ytemp[k] - y[k]) / jacob_epsilon;
      }

      // put x back
      x[j] = oldx;
    }
    if (_verbosity & DBG_JACOB)
      printMatrix("jacob", jacob, 2 * n, 6);

    // error (nx2)
    T dy[2 * n];
    for (int p = 0; p < 2 * n; p++) {
      dy[p] = y0[p] - y[p];
    }
    if (_verbosity & DBG_ER)
      printVector("dy", dy, 2 * n);

    // LM with fletcher improvement
    // inv(Jt*J + lambda*diag(Jt*J)) * Jt * dy
    T** JtJ = new T*[6];  // 6x6
    for (int p = 0; p < 6; p++)
      JtJ[p] = new T[6];
    matmult2DwTranspose<T>(jacob, n * 2, 6, true, jacob, n * 2, 6, false, JtJ);

    if (_verbosity & DBG_JACOB)
      printMatrix("JtJ", JtJ, 6, 6);

    // add lambda * diag(Jt*J) in-place
    for (int li = 0; li < 6; li++) {
      // add fletcher column-wise to JtJ
      T fletcher = lambda * JtJ[li][li];

      for (int lj = 0; lj < 6; lj++) {
        JtJ[lj][li] += fletcher;
      }
    }

    T** JtJ_i = new T*[6];  // 6x6
    for (int p = 0; p < 6; p++)
      JtJ_i[p] = new T[6];

    // OpenCV's invert
    // cv::Mat JtJ_i = cv::Mat(6, 6, cv::DataType<T>::type);
    // invert(JtJ, JtJ_i);

    // Compute SVD with type T
    // if (myinvert<T>(JtJ, 6, 6, JtJ_i)) {
    //  for (int p=0; p<2*n; p++)
    //      delete[] jacob[p];
    //  delete[] jacob;
    //  for (int p=0; p<6; p++) {
    //      delete[] JtJ[p];
    //  }
    //  delete[] JtJ;
    //  for (int p=0; p<6; p++)
    //      delete[] JtJ_i[p];
    //  delete[] JtJ_i;
    //  return -1;
    //}

    // Compute SVD with floats
    // Convert to floats
    float** f_JtJ = new float*[6];  // 6x6
    for (int p = 0; p < 6; p++) {
      f_JtJ[p] = new float[6];
      for (int pp = 0; pp < 6; pp++) {
        f_JtJ[p][pp] = JtJ[p][pp];
      }
    }
    float** f_JtJ_i = new float*[6];  // 6x6
    for (int p = 0; p < 6; p++)
      f_JtJ_i[p] = new float[6];
    if (myinvert<float>(f_JtJ, 6, 6, f_JtJ_i)) {
      for (int p = 0; p < 2 * n; p++) {
        delete[] jacob[p];
      }
      delete[] jacob;
      for (int p = 0; p < 6; p++) {
        delete[] JtJ[p];
        delete[] f_JtJ[p];
        delete[] JtJ_i[p];
        delete[] f_JtJ_i[p];
      }
      delete[] JtJ;
      delete[] f_JtJ;
      delete[] JtJ_i;
      delete[] f_JtJ_i;
      return -1;
    }
    // Convert back to T template type
    for (int p = 0; p < 6; p++) {
      for (int pp = 0; pp < 6; pp++) {
        JtJ_i[p][pp] = f_JtJ_i[p][pp];
      }
    }
    for (int p = 0; p < 6; p++) {
      delete[] f_JtJ[p];
      delete[] f_JtJ_i[p];
    }
    delete[] f_JtJ;
    delete[] f_JtJ_i;

    for (int p = 0; p < 6; p++) {
      delete[] JtJ[p];
    }
    delete[] JtJ;

    if (_verbosity & DBG_JACOB)
      printMatrix("JtJ_i", JtJ_i, 6, 6);

    T** JtJ_i_Jt = new T*[6];  // 6x2n
    for (int p = 0; p < 6; p++)
      JtJ_i_Jt[p] = new T[2 * n];
    matmult2DwTranspose<T>(JtJ_i, 6, 6, false, jacob, n * 2, 6, true, JtJ_i_Jt);
    for (int p = 0; p < 2 * n; p++)
      delete[] jacob[p];
    delete[] jacob;
    for (int p = 0; p < 6; p++)
      delete[] JtJ_i[p];
    delete[] JtJ_i;

    if (_verbosity & DBG_JACOB)
      printMatrix("JtJ_i_Jt", JtJ_i_Jt, 6, 2 * n);

    // linearize for matmult
    T JtJ_i_Jt_linear[6 * 2 * n];
    for (int p = 0; p < 6; p++) {
      for (int pp = 0; pp < n * 2; pp++) {
        JtJ_i_Jt_linear[(p * n * 2) + pp] = JtJ_i_Jt[p][pp];
      }
    }
    for (int p = 0; p < 6; p++)
      delete[] JtJ_i_Jt[p];
    delete[] JtJ_i_Jt;

    T dx[6];
    matmult(JtJ_i_Jt_linear, 6, 2 * n, dy, 2 * n,
            1,  // column vector
            dx);

    if (_verbosity & DBG_POSE_UPDATE)
      printVector("dx", dx, 6);

    // T errNorm = twonormsq(dy, 2*n);
    T errNorm = twonormsq(dx, 6);
    if (_verbosity & DBG_ER) {
      if (printints) {
        cout << "errNorm " << *(int*)&errNorm << endl;
      } else {
        cout << "errNorm " << errNorm << endl;
      }
    }

    if (errNorm > prevErrNorm)
      lambda *= 10;
    else
      lambda /= 10;
    lambda = MIN(lambda, lambda_max);
    lambda = MAX(lambda, lambda_min);
    if (_verbosity & DBG_LAMBDA) {
      if (printints) {
        cout << "lambda " << *(int*)&lambda << endl;
      } else {
        cout << "lambda " << lambda << endl;
      }
    }

    prevErrNorm = errNorm;

    // update pose
    for (int p = 0; p < 6; p++)
      x[p] = x[p] + dx[p];

    // break if error under threshold
    if (errNorm < min_er) {
      if (_verbosity & DBG_ER) {
        cout << "breaking because errNorm " << errNorm << " is smaller than "
             << min_er << endl;
      }
      break;
    }
  }

  if (_verbosity & DBG_FLOW)
    cout << "Did " << i + 1 << " GN iterations\n";
  printf("SeNtInAl,xy,%s,lm_%s_iterations_vs_numpts,%d,%d\n", __FUNCTION__,
         logprint.c_str(), n, i);

  // copy data back to x
  for (int i = 0; i < 6; i++) {
    _x[i] = x[i];
  }

  return i < LM_MAX_ITR;
}
