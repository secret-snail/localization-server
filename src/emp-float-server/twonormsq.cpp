#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "jlog.h"

//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

using namespace emp;
using namespace std;

// M mxn (M is secret, m n are public)
// N mmxnn
// res mxnn
Float BuildTwoNormSqCircuit(Float vect[], int sz) {
  Float sum = Float(0.0, PUBLIC);
  for (int i = 0; i < sz; i++) {
    sum = sum + (vect[i] * vect[i]);
  }
  return sum;
}
