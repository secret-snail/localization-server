#pragma once
//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include <cleartext-ref/svd.hpp>

using namespace std;

// The matrix version
template <typename T>
int myinvert(cv::Mat M, cv::Mat result) {
  int m = M.rows;
  int n = M.cols;
  assert(m >= n);
  assert(n == result.rows);
  assert(m == result.cols);

  // need to copy M because svd overwrites with u matrix
  T** in = new T*[m];
  for (int i = 0; i < m; i++) {
    in[i] = new T[n];
    for (int j = 0; j < n; j++) {
      in[i][j] = M.at<T>(i, j);
    }
  }

  T** out = new T*[n];
  for (int i = 0; i < n; i++) {
    out[i] = new T[m];
  }

  T* w = new T[m];    // aka sigma, only diag
  T** v = new T*[n];  // nxn
  for (int i = 0; i < n; i++)
    v[i] = new T[n];

  // Mat u(m, nm, type, buf);
  // Mat w(nm, 1, type, u.ptr() + m*nm*esz);
  // Mat vt(nm, n, type, w.ptr() + nm*esz);
  // SVD::compute(src, w, u, vt);
  // SVD::backSubst(w, u, vt, Mat(), _dst);

  // compute_svd(in, m, n, w, v);// in overwritten with u
  // svd_backsubstitute(in, w, v, m, n, NULL, out[0]); // why is out is 1D? -
  // because of b :(

  if (svdcmp(in, m, n, w, v)) {
    return -1;
  }
  T** u = in;

  // cout << "w\n";
  // printVector(w, m);

  // cout << "v\n";
  // printMatrix(v, n, n);
  // cout << endl;

  // cout << "u\n";
  // printMatrix(u, m, n);
  // cout << endl;

  // result = inv(in) = v*inv(w)*uT
  // inv(w)*uT
  for (int j = 0; j < n; j++) {
    if (w[j]) {
      for (int i = 0; i < m; i++) {
        u[i][j] /= w[j];  // note u indices do transpose
      }
    } else {
      for (int i = 0; i < m; i++)
        u[i][j] = 0;
    }
  }

  // cout << "u\n";
  // printMatrix(u, m, n);
  // cout << endl;

  // v*(w*uT) (don't use matmult so we can do transpose ourselves)
  for (int j = 0; j < n; j++) {
    for (int jj = 0; jj < m; jj++) {
      out[j][jj] = 0.0;
      for (int k = 0; k < n; k++) {
        out[j][jj] += v[j][k] * u[jj][k];  // note u indices do transpose
      }
    }
  }

  // cout << "result\n";
  // printMatrix(out, n, m);
  // cout << endl;

  // copy into result
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      result.at<T>(i, j) = out[i][j];

  for (int i = 0; i < m; i++)
    delete[] in[i];
  delete[] in;
  for (int i = 0; i < n; i++)
    delete[] out[i];
  delete[] out;
  delete[] w;
  for (int i = 0; i < n; i++)
    delete[] v[i];
  delete[] v;
  return 0;
}

template <typename T>
int myinvert(T** M, int m, int n, T** res) {
  assert(m >= n);

  // need to copy M because svd overwrites with u matrix
  T** in = new T*[m];
  for (int i = 0; i < m; i++) {
    in[i] = new T[n];
    for (int j = 0; j < n; j++) {
      // in[i][j]=M.at<T>(i,j);
      in[i][j] = M[i][j];
    }
  }

  T* w = new T[m];    // aka sigma, only diag
  T** v = new T*[n];  // nxn
  for (int i = 0; i < n; i++)
    v[i] = new T[n];

  if (svdcmp<T>(in, m, n, w, v)) {
    for (int i = 0; i < m; i++)
      delete[] in[i];
    delete[] in;
    delete[] w;
    for (int i = 0; i < n; i++)
      delete[] v[i];
    delete[] v;
    return -1;
  }
  T** u = in;

  // printVector("w", w, m);
  // printMatrix("v", v, n, n);
  // printMatrix("u", u, m, n);

  // res = inv(in) = v*inv(w)*uT
  // inv(w)*uT
  for (int j = 0; j < n; j++) {
    if (w[j]) {
      for (int i = 0; i < m; i++) {
        u[i][j] /= w[j];  // note u indices do transpose
      }
    } else {
      for (int i = 0; i < m; i++)
        u[i][j] = 0;
    }
  }

  // cout << "u\n";
  // printMatrix(u, m, n);
  // cout << endl;

  // v*(w*uT) (don't use matmult so we can do transpose ourselves)
  for (int j = 0; j < n; j++) {
    for (int jj = 0; jj < m; jj++) {
      res[j][jj] = 0.0;
      for (int k = 0; k < n; k++) {
        res[j][jj] += v[j][k] * u[jj][k];  // note u indices do transpose
      }
    }
  }

  // cout << "result\n";
  // printMatrix(res, n, m);
  // cout << endl;

  for (int i = 0; i < m; i++)
    delete[] in[i];
  delete[] in;
  delete[] w;
  for (int i = 0; i < n; i++)
    delete[] v[i];
  delete[] v;
  return 0;
}
