#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

#include <matmult.h>
#include <rodrigues.h>
#include <util.h>

using namespace emp;
using namespace std;

// Note: All 2D matrix inputs passed as single dimension array
// P=4xnumPoints [x1, x2, ...;  y1, y2, ...;  z1, z2, ...; 1, 1, ... ]
// x=6x1 [ rotation angles; translation] (column)
// K=3x3 camera intrinsic
// result=3xnumPoints preallocated, only first two rows are x,y.
//      last row needed for intermediate calculation (homogeneous)
//      [ u1, u2, ...;  v1, v2, ...;  1, 1, ... ]
void BuildProjectPointsCircuit(Float _P[], Float _x[], Float _K[], Float _res[],
                               int numPoints, bool skipOneGate = false) {
  Float one_gate = Float(1.0, PUBLIC);

  void* raw_memoryR = operator new[](12 * sizeof(Float));
  Float* R = static_cast<Float*>(raw_memoryR);

  BuildRodriguesCircuit(_x, R);

  // copy translation into last column of R
  R[3] = _x[3];
  R[7] = _x[4];
  R[11] = _x[5];

  // cout << "x:"<< endl;
  // printFloatVector(_x, 6);
  // cout <<"R:"<< endl;
  // printFloatVector(R, 12);
  // cout << endl;

  // res = K*R*P;
  void* raw_memoryscratch = operator new[](3 * 4 * sizeof(Float));
  Float* scratch = static_cast<Float*>(raw_memoryscratch);
  BuildMatmultCircuit(_K, 3, 3, R, 3, 4, scratch);
  BuildMatmultCircuit(scratch, 3, 4, _P, 4, numPoints, _res);

  // homogeneous coords to x,y
  for (int i = 0; i < numPoints; i++) {
    _res[i] = _res[i] / _res[2 * numPoints + i];
    _res[numPoints + i] = _res[numPoints + i] / _res[2 * numPoints + i];
    if (!skipOneGate) {
      _res[2 * numPoints + i] = one_gate;  // can be excluded for efficiency
    }
  }

  for (int i = 0; i < 12; i++) {
    R[i].~Float();
  }
  delete[] R;
  for (int i = 0; i < 3 * 4; i++) {
    scratch[i].~Float();
  }
  delete[] scratch;
}
