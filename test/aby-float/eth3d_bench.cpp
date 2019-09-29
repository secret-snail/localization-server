#include <opencv2/opencv.hpp>

#include <libgen.h>        // dirname
#include <linux/limits.h>  // PATH_MAX
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  // readlink

#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <type_traits>

#include <gaussnewtonlocalization.h>
#include <jlog.h>
#include <lmlocalization.h>
//#include <util.h>
#include <test_params.h>
#include <cleartext-ref/gaussnewtonlocalization.hpp>
#include <cleartext-ref/lmlocalization.hpp>
#include <localize_wrapper.hpp>
#include <test_harness.hpp>

using std::vector;

static_assert(std::is_same_v<decltype(&BuildAndRunGaussNewton), aby_localizer>,
              "ABY GN localizer function sig does not match what "
              "localize_wrapper expects. Tests will fail.");

static_assert(std::is_same_v<decltype(&BuildAndRunLM), aby_localizer>,
              "ABY LM localizer function sig does not match what "
              "localize_wrapper expects. Tests will fail.");

// These tests require cleartext functions to have this signature
typedef bool (*cleartext_localizer)(
    vector<cv::Point3f> objectPointsSubset,
    vector<cv::Point2f> imagePointsSubset, float f, float cx, float xy,
    float* x, /* initial guess for { r1, r2, r3, t1, t2, t3 } */
    string logprint);

static_assert(
    std::is_same_v<decltype(&gaussNewton<float>), cleartext_localizer>,
    "cleartext GN localizer function sig does not match what "
    "this test expects.");

static_assert(std::is_same_v<decltype(&lm<float>), cleartext_localizer>,
              "cleartext LM localizer function sig does not match what "
              "this test expects.");

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cout << "Usage: " << argv[0]
              << " <lm or gn> <num frames> <num trials> <max num pts>\n";
    return 1;
  }
  int num_frames = atoi(argv[2]);
  int num_trials = atoi(argv[3]);
  int max_num_pts = atoi(argv[4]);

  std::string lm_str("lm");
  cleartext_localizer clear_func =
      lm_str == argv[1] ? lm<float> : gaussNewton<float>;
  aby_localizer secure_func =
      lm_str == argv[1] ? BuildAndRunLM : BuildAndRunGaussNewton;
  std::string log_str = lm_str == argv[1] ? "lm" : "gn";
  bool silent = false;

  eth3d_test_harness(
      num_frames, num_trials, max_num_pts,
      [&](cv::Mat rvec, cv::Mat tvec, cv::Mat cameraMatrix, cv::Mat distCoeffs,
          vector<cv::Point3f> objectPointsSubset,
          vector<cv::Point2f> imagePointsSubset, vector<float>& res,
          const vector<float>& groundTruthRes) {
        vector<float> initialGuess(6);
        for (int i = 0; i < 3; ++i) {
          initialGuess[i] = rvec.at<float>(i);
          initialGuess[i + 3] = tvec.at<float>(i);
        }
        float f = cameraMatrix.at<float>(0, 0);
        float cx = cameraMatrix.at<float>(0, 2);
        float cy = cameraMatrix.at<float>(1, 2);

        bool out = clear_func(objectPointsSubset, imagePointsSubset, f, cx, cy,
                              &initialGuess[0], "clearlm");
        for (int i = 0; i < 6; ++i) {
          res[i] = initialGuess[i];
        }
        return out;
      },
      [&](cv::Mat rvec, cv::Mat tvec, cv::Mat cameraMatrix, cv::Mat distCoeffs,
          vector<cv::Point3f> objectPointsSubset,
          vector<cv::Point2f> imagePointsSubset, vector<float>& res,
          const vector<float>& groundTruthRes) {
        bool result_matches = true;
        int num_pts = objectPointsSubset.size();

        for (auto& share_type : ctypes) {
          if (!silent) {
            std::cout << "testing secure " << cnames[share_type] << '\n';
          }

          std::thread bob([&]() {
            vector<float> res(6, 0);
            rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
            tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0
            aby_localize_wrapper(CLIENT, share_type, rvec, tvec, cameraMatrix,
                                 distCoeffs, objectPointsSubset,
                                 imagePointsSubset, res, secure_func);
          });

          vector<float> res(6, 0);
          rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
          tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
          tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0

          time_point<high_resolution_clock> tic = high_resolution_clock::now();

          aby_localize_wrapper(SERVER, share_type, rvec, tvec, cameraMatrix,
                               distCoeffs, objectPointsSubset,
                               imagePointsSubset, res, secure_func);
          bob.join();

          time_point<high_resolution_clock> toc = high_resolution_clock::now();

          if (!silent) {
            cout << "secure result: ";
            for (auto const& f : res)
              cout << f << ' ';
          }

          for (int i = 0; i < 6; ++i) {
            result_matches &= withinRel(res[i], groundTruthRes[i]);
          }

          if (result_matches) {
            if (!silent) {
              cout << "aby converged!\n";
            }
            std::string timer_name = "aby_" + cnames[share_type] + "_float_" +
                                     log_str + "_time_vs_points";
#if PPL_FLOW == PPL_FLOW_DO
            timer_name += "_dataobl";
#endif
            MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, timer_name.c_str(),
                num_pts,
                std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                        .count() /
                    1000000.0);
            // MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__,
            //     (timer_name + "_per_loc_itr").c_str(), num_pts,
            //     std::chrono::duration_cast<std::chrono::microseconds>(toc -
            //                                                           tic)
            //             .count() /
            //         (1000000.0 * num_loc_iterations));

          } else {
            if (!silent) {
              MSG("Not printing timing - not converged.\n");
            }
          }
        }
        return result_matches;
      },
      silent);
}