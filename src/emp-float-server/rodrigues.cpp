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

#include <rodrigues.h>
#include <trigfuncs.h>

using namespace emp;
using namespace std;

// r=3x1, R=3x4 (only 3x3 is used here)
// values are secret, sizes are public
void BuildRodriguesCircuit(Float r[], Float R[]) {
  Float one_gate = Float(1.0, PUBLIC);
  Float zero_gate = Float(0.0, PUBLIC);

  // calculate theta
  Float rzerosq = r[0] * r[0];
  Float ronesq = r[1] * r[1];
  Float rtwosq = r[2] * r[2];
  Float theta1 = rzerosq + ronesq;
  Float theta2 = theta1 + rtwosq;
  Float theta = theta2.sqrt();

  Float co = BuildCosCircuit(theta);
  Float si = BuildSinCircuit(theta);

  Float c1 = one_gate - co;

  // if theta == 0 then itheta = 0
  Bit thetaiszero = theta.equal(zero_gate);
  Float itheta = (one_gate / theta).If(thetaiszero, zero_gate);

  Float x = r[0] * itheta;
  Float y = r[1] * itheta;
  Float z = r[2] * itheta;

  Float xx = x * x;
  Float yy = y * y;
  Float zz = z * z;

  Float xy = x * y;
  Float xz = x * z;
  Float yz = y * z;

  Float sx = si * x;
  Float sy = si * y;
  Float sz = si * z;

  R[0] = (c1 * xx) + co;
  R[1] = (c1 * xy) - sz;
  R[2] = (c1 * xz) + sy;

  R[4] = (c1 * xy) + sz;
  R[5] = (c1 * yy) + co;
  R[6] = (c1 * yz) - sx;

  R[8] = (c1 * xz) - sy;
  R[9] = (c1 * yz) + sx;
  R[10] = (c1 * zz) + co;
}
