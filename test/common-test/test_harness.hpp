#include <opencv2/opencv.hpp>

#include <libgen.h>        // dirname
#include <linux/limits.h>  // PATH_MAX
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  // readlink

#include <iostream>
#include <random>
#include <string>

#include <jlog.h>
#include <printutil.h>
#include <test_params.h>
#include <util.h>
#include <eth3d_features.hpp>

using std::cout;
using std::vector;

constexpr const int seed = 0x666;

bool withinRel(float v, float t) {
  bool res = true;
  res &= fabs(v - t) <= std::max(localization_tol_rel,
                                 fabs(localization_tol_rel * std::max(v, t)));
  res &= v != std::numeric_limits<float>::max();
  res &= v != std::numeric_limits<float>::min();
  res &= v != std::numeric_limits<float>::infinity();
  res &= v != std::numeric_limits<float>::quiet_NaN();
  return res;
}

typedef std::function<bool(
    cv::Mat rvec, cv::Mat tvec, cv::Mat cameraMatrix, cv::Mat distCoeffs,
    vector<cv::Point3f> objectPointsSubset,
    vector<cv::Point2f> imagePointsSubset, vector<float>& res,
    const vector<float>& groundTruthRes)>
    localize_test_cb;

int eth3d_test_harness(int num_frames, int num_trials, int max_num_pts,
                       localize_test_cb cleartext_localize_func,
                       localize_test_cb secure_localize_func, bool silent) {
  uint32_t cv_successes = 0;
  uint32_t cleartext_successes = 0;
  uint32_t secure_successes = 0;
  uint32_t total_runs = 0;

  float f = 3408.57;
  float cx = 3114.7;
  float cy = 2070.92;
  float _cM[] = {f, 0, cx, 0, f, cy, 0, 0, 1};
  cv::Mat cameraMatrix = cv::Mat(3, 3, cv::DataType<float>::type, _cM);
  cv::Mat distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
  cv::Mat rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  cv::Mat tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
  vector<cv::Point2f> imagePoints;
  vector<cv::Point3f> objectPoints;

  char result[PATH_MAX];
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  const char* path;
  if (count != -1) {
    path = dirname(result);
  } else {
    std::cerr << "cant find path\n";
    return 1;
  }
  std::string base_path = string(path) + "/../../data-eth3d/";

  // Multiple locations affects the convergence properties of the various gd
  // algorithms Single location is better for comparison.
  // for (uint32_t l = 0; l < eth3d_locations.size(); l++) {
  for (int l = 0; l < 1; l++) {
    auto feats = ETH3DFeatures<float>(base_path, eth3d_locations[l]);

    int test_num_frames = MIN(feats.numberOfFrames(), num_frames);

    for (int frame = 0; frame < test_num_frames; frame++) {
      imagePoints.clear();
      objectPoints.clear();
      feats.imageFeatures(imagePoints);
      feats.worldFeatures(objectPoints);
      auto gtpose = feats.getGroundTruthPose();
      vector<float> initialGuess = {0, 0, 0, 0, 0, 1};
      // vector<float> initialGuess = {res.first[0], res.first[1], res.first[2],
      //         res.second[0], res.second[1], res.second[2]};
      // uint32_t numPts = imagePoints.size();
      // MSG("Found %d points\n", numPts);
      // for (int i=0; i<imagePoints.size(); i++) {
      //     MSG("Point %d 2d %f, %f\n", i, imagePoints[i].x, imagePoints[i].y);
      //     MSG("Point %d 3d %f, %f, %f\n", i, objectPoints[i].x,
      //     objectPoints[i].y, objectPoints[i].z);
      // }
      if (!silent) {
        printVector("Ground Truth Pose [rotation; translation]: ", &gtpose[0],
                    6);
      }

      int max_image_points = imagePoints.size();
      int test_num_pts = MIN(max_image_points, max_num_pts);

      // std::random_device dev;
      // std::mt19937 rng(dev());
      std::mt19937 rng(seed);
      std::uniform_int_distribution<std::mt19937::result_type> dist(
          0, max_image_points);

      for (int num_pts = 6; num_pts <= test_num_pts;
           num_pts += MAX((test_num_pts - 6) / 5, 1)) {  // at most 5 intervals

        for (int t = 0; t < num_trials; t++) {
          if (!silent) {
            cout << "Trial " << t << " with " << num_pts
                 << " randomly selected points on frame " << frame
                 << " from location " << l << "\n";
          }
          ++total_runs;

          // randomly sample i feature pairs from total
          vector<cv::Point2f> imagePointsSubset;
          vector<cv::Point3f> objectPointsSubset;
          for (int p = 0; p < num_pts; p++) {
            int r = dist(rng);
            imagePointsSubset.push_back(imagePoints[r]);
            objectPointsSubset.push_back(objectPoints[r]);
          }

          vector<float> opencvRes(6);
          {
            if (!silent) {
              cout << "\ntesting opencv\n";
            }
            distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
            rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0
            CLOCK(opencv);
            TIC(opencv);
            // OpenCV PnP method
            cv::solvePnP(objectPointsSubset, imagePointsSubset, cameraMatrix,
                         distCoeffs, rvec, tvec, false, cv::SOLVEPNP_ITERATIVE);
            TOC(opencv);
            if (!silent) {
              cout << "opencv result:" << endl;
              cout << rvec << endl;
              cout << tvec << endl << endl;
            }
            bool good = true;
            for (int i = 0; i < 3; ++i) {
              good &= withinRel(rvec.at<float>(i), gtpose[i]);
              good &= withinRel(tvec.at<float>(i), gtpose[i + 3]);
              opencvRes[i] = rvec.at<float>(i);
              opencvRes[i + 3] = tvec.at<float>(i);
            }
            if (good) {
              cout << "opencv converged!\n";
              ++cv_successes;
            } else {
              cout << "opencv did not converge.\n";
            }

            // If opencv does not converge to ground truth but does converge to
            // some reasonable value, we can proceed. Inf is not a reasonable
            // value and we should retry.
            if (std::any_of(opencvRes.begin(), opencvRes.end(), [](float x) {
                  return x == std::numeric_limits<float>::max();
                })) {
              continue;
            }
          }

          {
            if (!silent) {
              cout << "testing cleartext\n";
            }
            vector<float> res(6, 0);
            rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0
            CLOCK(cleartext);
            TIC(cleartext);
            cleartext_localize_func(rvec, tvec, cameraMatrix, distCoeffs,
                                    objectPointsSubset, imagePointsSubset, res,
                                    opencvRes);
            TOC(cleartext);
            if (!silent) {
              printVector("cleartext result:\n", &res[0], 6);
            }
            bool good = true;
            for (int i = 0; i < 6; ++i) {
              good &= withinRel(res[i], opencvRes[i]);
            }
            if (good) {
              cout << "cleartext converged!\n";
              ++cleartext_successes;
            } else {
              MSG("cleartext did not converge to same value as opencv.\n");
              continue;
            }
          }

          {
            if (!silent) {
              cout << "testing secure\n";
            }
            vector<float> res(6, 0);
            rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0
            if (secure_localize_func(rvec, tvec, cameraMatrix, distCoeffs,
                                     objectPointsSubset, imagePointsSubset, res,
                                     opencvRes)) {
              ++secure_successes;
            }
          }

          if (!silent) {
            cout << "\n\n";
          }
        }
      }
      feats.nextFrame();
    }
  }
  if (!silent) {
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "cv_successes", 0,
        cv_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "cleartext_successes", 0,
        cleartext_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "secure_successes", 0,
        secure_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "total_runs", 0, total_runs);
    MSG("OpenCV converged to the ground truth pose %d / %d (%f\%)\n",
        cv_successes, total_runs,
        static_cast<float>(cv_successes) / total_runs);
    MSG("cleartext converged to the opencv pose %d / %d (%f\%)\n",
        cleartext_successes, total_runs,
        static_cast<float>(cleartext_successes) / total_runs);
    MSG("secure localization converged to the opencv pose %d / %d (%f\%)\n",
        secure_successes, total_runs,
        static_cast<float>(secure_successes) / total_runs);
  }
  return 0;
}