#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

// #include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

#include <invert.h>
#include <svd.h>
#include <util.h>

using namespace emp;
using namespace std;

const uint32_t bitlen = 32;

void BuildInvertCircuit(Float* in[], int m, int n, Float* res[]) {
  assert(m >= n);
  Float zero_gate = Float(0.0, PUBLIC);

  Float* w = static_cast<Float*>(operator new[](n * sizeof(Float)));
  for (int i = 0; i < n; i++) {
    w[i] = Float(0.0, PUBLIC);
  }
  Float** v = new Float*[n];
  for (int i = 0; i < n; i++) {
    v[i] = static_cast<Float*>(operator new[](n * sizeof(Float)));
    for (int j = 0; j < n; j++) {
      v[i][j] = Float(0.0, PUBLIC);
    }
  }

  BuildSvdCircuit(in, m, n, w, v);

  Float** u = in;

  // cout << "w\n";
  // printFloatVector(w, m);

  // cout << "v\n";
  // printFloatMatrix(v, n, n);
  // cout << endl;

  // cout << "u\n";
  // printFloatMatrix(u, m, n);
  // cout << endl;

  // res = inv(in) = v*inv(w)*uT
  // inv(w)*uT
  for (int j = 0; j < n; j++) {
    // if (w[j]) {
    Bit ifw = w[j].equal(zero_gate);

    for (int i = 0; i < m; i++) {
      Float temp = u[i][j] / w[j];
      u[i][j] = temp.If(ifw, Float(0.0, PUBLIC));
    }
    //}
  }

  // v*(w*uT) (don't use matmult so we can do transpose ourselves)
  for (int j = 0; j < n; j++) {
    for (int jj = 0; jj < m; jj++) {
      res[j][jj] = Float(0.0, PUBLIC);
      for (int k = 0; k < n; k++) {
        res[j][jj] =
            res[j][jj] + (v[j][k] * u[jj][k]);  // note u indices do transpose
      }
    }
  }

  delete[] w;
  for (int i = 0; i < n; i++) {
    delete[] v[i];
  }
  delete[] v;
}
