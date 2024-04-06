#include <jlog.h>
#include <privacyconf.h>
#include <stdio.h>
#include <iostream>

#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

#include <invert.h>
#include <matmult.h>
#include <projectpoints.h>
#include <svd.h>
#include <twonormsq.h>
#include <util.h>

using namespace emp;
using namespace std;

// Note: All 2D matrix inputs passed as single dimension array
// Note: the constant 1's can be ignored (but space must be allocated)
//
// 3D World Points (4xnumPts) which look like:
//  _             _
// |  x1 x2 x3     |
// |  y1 y2 y3 ... |
// |  z1 z2 z3     |
// |_ 1  1  1     _|
// must be pass as single dimension:
// [x1, x2, ...;  y1, y2, ...;  z1, z2, ...; 1, 1, ... ]
//
// 2D Image Points (3xnumPts) which look like:
//  _             _
// |  u1 u2 u3     |
// |  v1 v2 v3 ... |
// |_ 1  1  1     _|
// must be passed as single dimension:
// [ u1, u2, ...;  v1, v2, ...;  1, 1, ... ]
//
// x (initial guess) : 6x1 [ rotation angles, translations ]
std::pair<bool, int> BuildGaussNewton(Float threeDPts[], Float y0[], int numPts,
                                      Float f, Float cx, Float cy, Float x[]) {
  // cout << "y0:\n";
  // printFloatVector(y0, 2*numPts);

  // add constant ones to threeDPts
  Float one_gate = Float(1.0, PUBLIC);
  for (int i = 0; i < numPts; i++) {
    threeDPts[3 * numPts + i] = one_gate;
  }

  // Camera params (3x3)
  Float K[9] = {f,
                Float(0.0, PUBLIC),
                cx,
                Float(0.0, PUBLIC),
                f,
                cy,
                Float(0.0, PUBLIC),
                Float(0.0, PUBLIC),
                Float(1.0, PUBLIC)};

  float jacob_epsilon =
      JACOB_EPSILON;  // this is needed, can't &JACOB_EPSILON directly

  Float min_er = Float(MIN_ER, PUBLIC);
  int i;
  for (i = 0; i < GN_MAX_ITR; i++) {

    // project points using x guess
    // share **yHomog = new share*[3*numPts];
    Float* yHomog =
        static_cast<Float*>(operator new[](3 * numPts * sizeof(Float)));
    BuildProjectPointsCircuit(threeDPts, x, K, yHomog, numPts, true);

    // throw away last "row" of result (constant ones),
    // interleave, and transpose y into 2nx1 vector
    // e.g. [x1; y1; x2; y2 ...]
    // share **y = new share*[2*numPts];
    Float* y = static_cast<Float*>(operator new[](2 * numPts * sizeof(Float)));
    for (int p = 0; p < numPts; ++p) {
      y[p * 2] = yHomog[p];
      y[(p * 2) + 1] = yHomog[p + numPts];
    }
    delete[] yHomog;
    // cout << "y:\n";
    // printFloatVector(y, 2*numPts);

    // calculate jacobian (real 2D matrix not linearized)
    //  _                                            _
    // | du1/dr1 du1/dr2 du1/dr3 du1/t1 du1/t2 du1/t3 |
    // | dv1/dr1 dv1/dr2 dv1/dr3 dv1/t1 dv1/t2 dv1/t3 |
    // | du2/dr1 du2/dr2 du2/dr3 du2/t1 du2/t2 du2/t3 |
    // |                      .                       |
    // |                      .                       |
    // |_                     . (2n)                 _|
    Float** jacob = new Float*[numPts * 2];
    for (int p = 0; p < numPts * 2; p++) {
      jacob[p] = static_cast<Float*>(operator new[](6 * sizeof(Float)));
    }
    for (int j = 0; j < 6; j++) {  // for each DOF
      // perturb x
      Float temp = x[j];
      x[j] = x[j] + Float(jacob_epsilon, PUBLIC);

      // project with epsilon
      Float* ytempHomog =
          static_cast<Float*>(operator new[](3 * numPts * sizeof(Float)));
      BuildProjectPointsCircuit(threeDPts, x, K, ytempHomog, numPts, true);

      // throw away last "row" of result,
      // reshape, and transpose y into 2nx1 vector
      // e.g. [x1; y1; x2; y2 ...]
      Float* ytemp =
          static_cast<Float*>(operator new[](2 * numPts * sizeof(Float)));
      for (int p = 0; p < numPts; ++p) {
        ytemp[p * 2] = ytempHomog[p];
        ytemp[(p * 2) + 1] = ytempHomog[p + numPts];
      }
      delete[] ytempHomog;

      // put into jacob
      for (int p = 0; p < 2 * numPts; ++p) {
        jacob[p][j] = (ytemp[p] - y[p]) / Float(jacob_epsilon, PUBLIC);
      }

      delete[] ytemp;

      // un-perturb x
      x[j] = temp;
    }
    // printFloatMatrix(jacob, 2*numPts, 6);

    // calculate error
    // dy = y0 - y;
    Float* dy = static_cast<Float*>(operator new[](2 * numPts * sizeof(Float)));
    for (int p = 0; p < 2 * numPts; p++) {
      dy[p] = y0[p] - y[p];
    }
    // cout << "dy:\n";
    // printFloatVector(dy, 2*numPts);
    delete[] y;

    // compute pseudo inverse to solve:
    // dy = J dx   ->   dx = Jdagger dy
    // Note: invert requires real
    // 2D arrays not linearized versions
    Float** jacobI = new Float*[6];  // 6x2n
    for (int p = 0; p < 6; p++) {
      jacobI[p] =
          static_cast<Float*>(operator new[](2 * numPts * sizeof(Float)));
    }

    BuildInvertCircuit(jacob, 2 * numPts, 6, jacobI);
    // printFloatMatrix(jacobI, 6, 2*numPts);

    // linearize for matmult
    Float* jacobILinear =
        static_cast<Float*>(operator new[](6 * 2 * numPts * sizeof(Float)));
    for (int p = 0; p < 6; p++) {
      for (int pp = 0; pp < 2 * numPts; pp++) {
        jacobILinear[(p * 2 * numPts) + pp] = jacobI[p][pp];
      }
    }
    Float* dx = static_cast<Float*>(operator new[](6 * sizeof(Float)));
    BuildMatmultCircuit(jacobILinear, 6, 2 * numPts, dy, 2 * numPts,
                        1,  // column vector
                        dx);
    // cout << "dx\n";
    // printFloatVector(dx, 6);
    delete[] dy;

    // cleanup jacobian
    for (int p = 0; p < numPts * 2; p++) {
      // individual jacob shares are deleted in invert->svd()
      delete[] jacob[p];
    }
    delete[] jacob;

    Float norm = BuildTwoNormSqCircuit(dx, 6);
    Float absnorm = norm.abs();

    for (int p = 0; p < 6; p++) {
      delete[] jacobI[p];
    }
    delete[] jacobI;
    delete[] jacobILinear;

    Bit er_flag = absnorm.less_than(min_er);

    // update pose - note: this is usually done after checking error
    for (int p = 0; p < 6; p++) {
#if PPL_FLOW == PPL_FLOW_DO  // set dx to zero if no error
      for (int bi = 0; bi < dx[p].size(); bi++)
        dx[p][i] = dx[p][i] & er_flag;
#endif
      x[p] = x[p] + dx[p];
    }
    delete[] dx;

#if PPL_FLOW != PPL_FLOW_DO
    // okay to reveal how many iterations it took
    // return if error under threshold
    //  abs(norm(dy, cv::NORM_L2SQR)) < MIN_ER
    if (er_flag.reveal<bool>()) {
      break;
    }
#endif
  }

  if (_verbosity & DBG_FLOW)
    cout << "Did " << i << " GN iterations\n";
  return {i < GN_MAX_ITR, i + 1};
}
