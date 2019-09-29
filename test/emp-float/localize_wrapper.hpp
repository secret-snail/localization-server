#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>

#include "emp-sh2pc/emp-sh2pc.h"

using std::vector;

typedef std::pair<bool, int> (*emp_localizer)(Float threeDPts[], Float y0[],
                                              int numPts, Float f, Float cx,
                                              Float cy, Float x[]);

std::pair<bool, int> emp_localize_wrapper(int party, cv::Mat rvec, cv::Mat tvec,
                                          cv::Mat cameraMatrix,
                                          cv::Mat distCoeffs,
                                          vector<cv::Point3f> objectPoints,
                                          vector<cv::Point2f> imagePoints,
                                          vector<float>& res,
                                          emp_localizer secure_localize_func) {

  float f = cameraMatrix.at<float>(0, 0);
  float cx = cameraMatrix.at<float>(0, 2);
  float cy = cameraMatrix.at<float>(1, 2);

  int numPts = objectPoints.size();
  Float* sobjectPoints =
      static_cast<Float*>(operator new[](4 * numPts * sizeof(Float)));
  Float* simagePoints =
      static_cast<Float*>(operator new[](3 * numPts * sizeof(Float)));
  for (int i = 0; i < numPts; i++) {
    // [x1, x2... ; y1, y2... ; z1, z2...]
    sobjectPoints[i] = Float(objectPoints[i].x, ALICE);
    sobjectPoints[numPts + i] = Float(objectPoints[i].y, ALICE);
    sobjectPoints[2 * numPts + i] = Float(objectPoints[i].z, ALICE);
    // sobjectPoints[3*numPts + i] = Float(1.0, PUBLIC);

    // [x1, y1; x2, y2; ...]
    simagePoints[2 * i] = Float(imagePoints[i].x, ALICE);
    simagePoints[2 * i + 1] = Float(imagePoints[i].y, ALICE);
    // simagePoints[2*numPts + i] = Float(1.0, PUBLIC);
  }

  Float sf = Float(f, ALICE);
  Float scx = Float(cx, ALICE);
  Float scy = Float(cy, ALICE);

  Float* sx = static_cast<Float*>(operator new[](6 * sizeof(Float)));
  for (int i = 0; i < 3; i++) {
    sx[i] = Float(rvec.at<float>(i), ALICE);
    sx[i + 3] = Float(tvec.at<float>(i), ALICE);
  }

  auto sres = secure_localize_func(sobjectPoints, simagePoints, numPts, sf, scx,
                                   scy, sx);

  for (int i = 0; i < 6; ++i) {
    res[i] = sx[i].reveal<double>();
  }

  delete[] sx;
  delete[] sobjectPoints;
  delete[] simagePoints;
  return sres;
}