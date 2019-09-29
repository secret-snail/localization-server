//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <jlog.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <random>
#include <string>

#include <libgen.h>        // dirname
#include <linux/limits.h>  // PATH_MAX
#include <unistd.h>        // readlink

//#include <catch2/catch_test_macros.hpp>
//#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <gaussnewtonlocalization.h>
#include <lmlocalization.h>
#include <printutil.h>
#include <util.h>
#include <cleartext-ref/gaussnewtonlocalization.hpp>
#include <cleartext-ref/lmlocalization.hpp>
#include <eth3d_features.hpp>

#include <fixed_point.h>

using std::cout;
using std::vector;

constexpr const int port = 8080;
constexpr const int seed = 0x666;
const static constexpr float localization_tol_abs = 0.05f;
const static constexpr float localization_tol_rel = 0.05f;

bool withinRel(float v, float t) {
  return fabs(v - t) <= std::max(localization_tol_rel,
                                 fabs(localization_tol_rel * std::max(v, t)));
}

template <typename T>
bool convert_and_localize(vector<cv::Point3_<float>> _objectPoints,
                          vector<cv::Point_<float>> _imagePoints, float _f,
                          float _cx, float _cy, T* _x,
                          auto cleartext_localize_func, vector<float>& gtpose,
                          string logprint = "") {

  T f = _f;
  T cx = _cx;
  T cy = _cy;
  vector<cv::Point3_<T>> objectPoints;
  vector<cv::Point_<T>> imagePoints;
  for (auto float_point : _objectPoints) {
    objectPoints.push_back({float_point.x, float_point.y, float_point.z});
  }
  for (auto float_point : _imagePoints) {
    imagePoints.push_back({float_point.x, float_point.y});
  }

  cleartext_localize_func(objectPoints, imagePoints, f, cx, cy, _x, logprint);

  bool good = true;
  for (int i = 0; i < 6; ++i) {
    float fx = _x[i];  // convert to float
    good &= withinRel(fx, gtpose[i]);
  }
  return good;
}

int eth3d_localize(NetIO* io, int party, int num_frames, int num_trials,
                   int max_num_pts, bool silent, std::string log_str = "") {
  setup_semi_honest(io, party);
  uint32_t cv_successes = 0;
  uint32_t gn_double_successes = 0;
  uint32_t lm_double_successes = 0;
  uint32_t gn_float_successes = 0;
  uint32_t lm_float_successes = 0;
  uint32_t gn_fixed32_successes = 0;
  uint32_t lm_fixed32_successes = 0;
  uint32_t gn_fixed64_successes = 0;
  uint32_t lm_fixed64_successes = 0;
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
    error("cant find path\n");
    return 1;
  }
  std::string base_path = string(path) + "/../../data-eth3d/";

  for (uint32_t l = 0; l < eth3d_locations.size(); l++) {
    // for (int l=0; l<1; l++) {
    auto feats = ETH3DFeatures<float>(base_path, eth3d_locations[l]);

    int test_num_frames = MIN(feats.numberOfFrames(), num_frames);

    for (int frame = 0; frame < test_num_frames; frame++) {
      imagePoints.clear();
      objectPoints.clear();
      feats.imageFeatures(imagePoints);
      feats.worldFeatures(objectPoints);
      auto gtpose = feats.getGroundTruthPose();
      vector<float> initialGuess = {0, 0, 0, 0, 0, 0};
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

          // randomly sample i feature pairs from total
          vector<cv::Point2f> imagePointsSubset;
          vector<cv::Point3f> objectPointsSubset;
          for (int p = 0; p < num_pts; p++) {
            int r = dist(rng);
            imagePointsSubset.push_back(imagePoints[r]);
            objectPointsSubset.push_back(objectPoints[r]);
          }

          {
            if (!silent) {
              cout << "\ntesting opencv\n";
            }
            distCoeffs = cv::Mat::zeros(4, 1, cv::DataType<float>::type);
            rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            // OpenCV PnP method
            cv::solvePnP(objectPointsSubset, imagePointsSubset, cameraMatrix,
                         distCoeffs, rvec, tvec, false, cv::SOLVEPNP_ITERATIVE);
            if (!silent) {
              cout << "opencv result:" << endl;
              cout << rvec << endl;
              cout << tvec << endl << endl;
            }
            bool fail = false;
            for (int i = 0; i < 3; ++i) {
              fail |= withinRel(rvec.at<float>(i), gtpose[i]);
              fail |= withinRel(tvec.at<float>(i), gtpose[i + 3]);
            }
            if (!fail) {
              ++cv_successes;
            } else {
              continue;  // if opencv does not converge, ignore this run
            }
          }

          {
            if (!silent) {
              cout << "testing cleartext - gn - double\n";
            }
            vector<double> initialGuessCopy;
            for (auto v : initialGuess) {
              initialGuessCopy.push_back(v);
            }
            bool good = convert_and_localize(objectPoints, imagePoints, f, cx,
                                             cy, &initialGuessCopy[0],
                                             gaussNewton<double>, gtpose);
            if (!silent) {
              printVector("result:\n", &initialGuessCopy[0], 6);
            }
            if (good) {
              ++gn_double_successes;
            }
          }

          {
            if (!silent) {
              cout << "testing cleartext - lm - double\n";
            }
            vector<double> initialGuessCopy;
            for (auto v : initialGuess) {
              initialGuessCopy.push_back(v);
            }
            bool good =
                convert_and_localize(objectPoints, imagePoints, f, cx, cy,
                                     &initialGuessCopy[0], lm<double>, gtpose);
            if (!silent) {
              printVector("cleartext result:\n", &initialGuessCopy[0], 6);
            }
            bool fail = false;
            for (int i = 0; i < 6; ++i) {
              fail |= withinRel(initialGuessCopy[i], gtpose[i]);
            }
            if (!fail) {
              ++lm_double_successes;
            }
          }

          {
            if (!silent) {
              cout << "testing cleartext - gn - float\n";
            }
            vector<float> initialGuessCopy;
            for (auto v : initialGuess) {
              initialGuessCopy.push_back(v);
            }
            bool good = convert_and_localize(objectPoints, imagePoints, f, cx,
                                             cy, &initialGuessCopy[0],
                                             gaussNewton<float>, gtpose);
            if (!silent) {
              printVector("result:\n", &initialGuessCopy[0], 6);
            }
            if (good) {
              ++gn_float_successes;
            }
          }

          {
            if (!silent) {
              cout << "testing cleartext - lm - float\n";
            }
            vector<float> initialGuessCopy;
            for (auto v : initialGuess) {
              initialGuessCopy.push_back(v);
            }
            bool good =
                convert_and_localize(objectPoints, imagePoints, f, cx, cy,
                                     &initialGuessCopy[0], lm<float>, gtpose);
            if (!silent) {
              printVector("cleartext result:\n", &initialGuessCopy[0], 6);
            }
            bool fail = false;
            for (int i = 0; i < 6; ++i) {
              fail |= withinRel(initialGuessCopy[i], gtpose[i]);
            }
            if (!fail) {
              ++lm_float_successes;
            }
          }

          {
            if (!silent) {
              cout << "testing cleartext - gn - fixed64\n";
            }
            vector<fixed_point<int64_t, 32>> initialGuessCopy;
            for (auto v : initialGuess) {
              initialGuessCopy.push_back(v);
            }
            bool good = convert_and_localize(
                objectPoints, imagePoints, f, cx, cy, &initialGuessCopy[0],
                gaussNewton<fixed_point<int64_t, 32>>, gtpose);
            if (!silent) {
              printVector("result:\n", &initialGuessCopy[0], 6);
            }
            if (good) {
              ++gn_fixed64_successes;
            }
          }
          {
            if (!silent) {
              cout << "testing cleartext - lm - fixed64\n";
            }
            vector<fixed_point<int64_t, 32>> initialGuessCopy;
            for (auto v : initialGuess) {
              initialGuessCopy.push_back(v);
            }
            bool good = convert_and_localize(
                objectPoints, imagePoints, f, cx, cy, &initialGuessCopy[0],
                lm<fixed_point<int64_t, 32>>, gtpose);
            if (!silent) {
              printVector("result:\n", &initialGuessCopy[0], 6);
            }
            if (good) {
              ++lm_fixed64_successes;
            }
          }

          {
            if (!silent) {
              cout << "testing cleartext - gn - fixed32\n";
            }
            vector<fixed_point<int32_t, 16>> initialGuessCopy;
            for (auto v : initialGuess) {
              initialGuessCopy.push_back(v);
            }
            bool good = convert_and_localize(
                objectPoints, imagePoints, f, cx, cy, &initialGuessCopy[0],
                gaussNewton<fixed_point<int32_t, 16>>, gtpose);
            if (!silent) {
              printVector("result:\n", &initialGuessCopy[0], 6);
            }
            if (good) {
              ++gn_fixed32_successes;
            }
          }
          {
            if (!silent) {
              cout << "testing cleartext - lm - fixed64\n";
            }
            vector<fixed_point<int32_t, 16>> initialGuessCopy;
            for (auto v : initialGuess) {
              initialGuessCopy.push_back(v);
            }
            bool good = convert_and_localize(
                objectPoints, imagePoints, f, cx, cy, &initialGuessCopy[0],
                lm<fixed_point<int32_t, 16>>, gtpose);
            if (!silent) {
              printVector("result:\n", &initialGuessCopy[0], 6);
            }
            if (good) {
              ++lm_fixed32_successes;
            }
          }

          ++total_runs;
        }
      }
      feats.nextFrame();
    }
  }
  if (!silent) {
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "cv_successes", 0,
        cv_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "gn_double_successes", 0,
        gn_double_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "lm_double_successes", 0,
        lm_double_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "gn_float_successes", 0,
        gn_float_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "lm_float_successes", 0,
        lm_float_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "gn_fixed32_successes", 0,
        gn_fixed32_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "lm_fixed32_successes", 0,
        lm_fixed32_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "gn_fixed64_successes", 0,
        gn_fixed64_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "lm_fixed64_successes", 0,
        lm_fixed64_successes);
    MSG("SeNtInAl,xy,%s,%s,%d,%u\n", __FUNCTION__, "total_runs", 0, total_runs);
  }
  finalize_semi_honest();
  return 0;
}

int main(int argc, char** argv) {
  if (argc != 4) {
    MSG("Usage: %s <num frames> <num trials> <max num pts>\n", argv[0]);
    return 1;
  }
  int num_frames = atoi(argv[1]);
  int num_trials = atoi(argv[2]);
  int max_num_pts = atoi(argv[3]);

  std::thread bob([&]() {
    NetIO* io = new NetIO("127.0.0.1", port);
    eth3d_localize(io, BOB, num_frames, num_trials, max_num_pts, true);
    delete io;
  });

  NetIO* io = new NetIO(nullptr, port);
  eth3d_localize(io, ALICE, num_frames, num_trials, max_num_pts, false);
  bob.join();
  delete io;
}

// using namespace std;
//
// typedef fixed_point<int64_t,32> baset;
//
// int main(int argc, char** argv) {
//     vector<cv::Point_<float>> fimagePoints;
//     vector<cv::Point3_<float>> fobjectPoints;
//     vector<cv::Point_<baset>> bimagePoints;
//     vector<cv::Point3_<baset>> bobjectPoints;
//
//     // ETH3D setup
//     baset ff=3408.57;
//     baset fcx=3114.7;
//     baset fcy=2070.92;
//     baset bf=3408.57;
//     baset bcx=3114.7;
//     baset bcy=2070.92;
//
//     char datadir[PATH_MAX];
//     strncpy(datadir, argv[0], sizeof(datadir));
//     dirname(datadir);
//     std::string bindir(datadir);
//     std::string base_path = bindir + "/../../data-eth3d/";
//
//     auto ffeats = ETH3DFeatures<float>(base_path, eth3d_locations[0]);
//     ffeats.imageFeatures(fimagePoints);
//     ffeats.worldFeatures(fobjectPoints);
//     auto fgtpose = ffeats.getGroundTruthPose();
//     vector<float> finitialGuess = {0, 0, 0, 0, 0, 0};
//
//     auto bfeats = ETH3DFeatures<baset>(base_path, eth3d_locations[0]);
//     bfeats.imageFeatures(bimagePoints);
//     bfeats.worldFeatures(bobjectPoints);
//     auto bgtpose = bfeats.getGroundTruthPose();
//     vector<baset> binitialGuess = {0, 0, 0, 0, 0, 0};
//     //vector<float> initialGuess = {res.first[0], res.first[1], res.first[2],
//     //        res.second[0], res.second[1], res.second[2]};
//     //uint32_t numPts = imagePoints.size();
//     //MSG("Found %d points\n", numPts);
//     //for (int i=0; i<5; i++) {
//         //cout << "Point "<<i<<" 2d "<<imagePoints[i].x<<"
//         "<<imagePoints[i].y<< endl;
//         //cout << "Point "<<i<<" 3d "<<objectPoints[i].x<<"
//         "<<objectPoints[i].y<<" "<<objectPoints[i].z <<endl;
//     //}
//     printVector("Ground Truth Pose [rotation; translation]: ", &fgtpose[0],
//     6);
//
//
//     int maxPoints=fimagePoints.size();
//
//     //for (int i=6; i<maxPoints; i++) {
//     for (int i=6; i<276; i++) { // gn breaks after 276
//         vector<cv::Point_<float>> fimagePointsSubset(fimagePoints.begin(),
//         fimagePoints.begin() + i); vector<cv::Point3_<float>>
//         fobjectPointsSubset(fobjectPoints.begin(), fobjectPoints.begin() +
//         i); vector<cv::Point_<baset>>
//         bimagePointsSubset(bimagePoints.begin(), bimagePoints.begin() + i);
//         vector<cv::Point3_<baset>> bobjectPointsSubset(bobjectPoints.begin(),
//         bobjectPoints.begin() + i); cout << "Testing with " <<
//         bimagePointsSubset.size() << " points\n";
//
//         {   // requires float, does not work with fixed point
//             cout << "\ntesting opencv\n";
//             float _cM[] = {ff, 0, fcx,
//                            0, ff, fcy,
//                            0, 0, 1};
//             cv::Mat cameraMatrix = cv::Mat(3, 3, cv::DataType<float>::type,
//             _cM); cv::Mat distCoeffs =
//             cv::Mat::zeros(4,1,cv::DataType<float>::type); cv::Mat rvec =
//             cv::Mat::zeros(3,1,cv::DataType<float>::type); cv::Mat tvec =
//             cv::Mat::zeros(3,1,cv::DataType<float>::type); CLOCK(opencv);
//             TIC(opencv);
//             // OpenCV PnP method
//             cv::solvePnP(fobjectPointsSubset, fimagePointsSubset,
//             cameraMatrix, distCoeffs, rvec, tvec, false,
//             cv::SOLVEPNP_ITERATIVE); TOC(opencv); cout << "opencv result:" <<
//             endl; cout << rvec << endl; cout << tvec << endl << endl;
//
//             float error[6];
//             error[0] = fgtpose[0] - rvec.at<float>(0);
//             error[1] = fgtpose[1] - rvec.at<float>(1);
//             error[2] = fgtpose[2] - rvec.at<float>(2);
//             error[3] = fgtpose[3] - tvec.at<float>(0);
//             error[4] = fgtpose[4] - tvec.at<float>(1);
//             error[5] = fgtpose[5] - tvec.at<float>(2);
//             printf("SeNtInAl,xy,%s,opencv_error_vs_numpts,%d,%f\n",
//             __FUNCTION__, i, twonormsq(error, 6));
//         }
//
//#ifdef PPL_GN
//         {
//             cout << "testing gn - cleartext - float\n";
//             vector<float> finitialGuessCopy = finitialGuess;
//             CLOCK(gn_float_cleartext);
//             TIC(gn_float_cleartext);
//             gaussNewton<float>(fobjectPointsSubset, fimagePointsSubset, ff,
//             fcx, fcy, &finitialGuessCopy[0], "float");
//             TOC(gn_float_cleartext); printVector("[rotation; translation]",
//             &finitialGuessCopy[0], 6);
//
//             float error[6];
//             error[0] = fgtpose[0] -  finitialGuessCopy[0];
//             error[1] = fgtpose[1] -  finitialGuessCopy[1];
//             error[2] = fgtpose[2] -  finitialGuessCopy[2];
//             error[3] = fgtpose[3] - finitialGuessCopy[3];
//             error[4] = fgtpose[4] - finitialGuessCopy[4];
//             error[5] = fgtpose[5] - finitialGuessCopy[5];
//             printf("SeNtInAl,xy,%s,gn_float_error_vs_numpts,%d,%f\n",
//             __FUNCTION__, i, twonormsq(error, 6));
//         }
//         {
//             cout << "testing gn - cleartext - baset\n";
//             vector<baset> binitialGuessCopy = binitialGuess;
//             CLOCK(gn_baset_cleartext);
//             TIC(gn_baset_cleartext);
//             gaussNewton<baset>(bobjectPointsSubset, bimagePointsSubset, bf,
//             bcx, bcy, &binitialGuessCopy[0], "baset");
//             TOC(gn_baset_cleartext); printVector("[rotation; translation]",
//             &binitialGuessCopy[0], 6);
//
//             baset error[6];
//             error[0] = bgtpose[0] -  binitialGuessCopy[0];
//             error[1] = bgtpose[1] -  binitialGuessCopy[1];
//             error[2] = bgtpose[2] -  binitialGuessCopy[2];
//             error[3] = bgtpose[3] - binitialGuessCopy[3];
//             error[4] = bgtpose[4] - binitialGuessCopy[4];
//             error[5] = bgtpose[5] - binitialGuessCopy[5];
//             printf("SeNtInAl,xy,%s,gn_baset_error_vs_numpts,%d,%f\n",
//             __FUNCTION__, i, (double)twonormsq(error, 6));
//         }
//#endif
//#ifdef PPL_LM
//         {
//             cout << "testing lm - cleartext - float\n";
//             vector<float> finitialGuessCopy = finitialGuess;
//             CLOCK(lm_float_cleartext);
//             TIC(lm_float_cleartext);
//             lm<float>(fobjectPointsSubset, fimagePointsSubset, ff, fcx, fcy,
//             &finitialGuessCopy[0]); TOC(lm_float_cleartext);
//             printVector("[rotation; translation]", &finitialGuessCopy[0], 6);
//
//             float error[6];
//             error[0] = fgtpose[0] -  finitialGuessCopy[0];
//             error[1] = fgtpose[1] -  finitialGuessCopy[1];
//             error[2] = fgtpose[2] -  finitialGuessCopy[2];
//             error[3] = fgtpose[3] - finitialGuessCopy[3];
//             error[4] = fgtpose[4] - finitialGuessCopy[4];
//             error[5] = fgtpose[5] - finitialGuessCopy[5];
//             printf("SeNtInAl,xy,%s,lm_float_error_vs_numpts,%d,%f\n",
//             __FUNCTION__, i, twonormsq(error, 6));
//         }
//         {
//             cout << "testing lm - cleartext - baset\n";
//             vector<baset> binitialGuessCopy = binitialGuess;
//             CLOCK(lm_baset_cleartext);
//             TIC(lm_baset_cleartext);
//             lm<baset>(bobjectPointsSubset, bimagePointsSubset, bf, bcx, bcy,
//             &binitialGuessCopy[0]); TOC(lm_baset_cleartext);
//             printVector("[rotation; translation]", &binitialGuessCopy[0], 6);
//
//             baset error[6];
//             error[0] = bgtpose[0] -  binitialGuessCopy[0];
//             error[1] = bgtpose[1] -  binitialGuessCopy[1];
//             error[2] = bgtpose[2] -  binitialGuessCopy[2];
//             error[3] = bgtpose[3] - binitialGuessCopy[3];
//             error[4] = bgtpose[4] - binitialGuessCopy[4];
//             error[5] = bgtpose[5] - binitialGuessCopy[5];
//             printf("SeNtInAl,xy,%s,lm_baset_error_vs_numpts,%d,%f\n",
//             __FUNCTION__, i, (double)twonormsq(error, 6));
//         }
//#endif
//         cout << "\n\n";
//     }
//
//     return 0;
// }
