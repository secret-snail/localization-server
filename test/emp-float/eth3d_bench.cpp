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

#include <gaussnewtonlocalization.h>
#include <lmlocalization.h>
#include <printutil.h>
#include <test_params.h>
#include <util.h>
#include <cleartext-ref/gaussnewtonlocalization.hpp>
#include <cleartext-ref/lmlocalization.hpp>
#include <localize_wrapper.hpp>
#include <test_harness.hpp>

using std::cout;
using std::vector;

constexpr const int port = 8080;

static_assert(std::is_same_v<decltype(&BuildGaussNewton), emp_localizer>,
              "EMP GN localizer function sig does not match what "
              "localize_wrapper expects. Tests will fail.");

static_assert(std::is_same_v<decltype(&BuildLM), emp_localizer>,
              "EMP LM localizer function sig does not match what "
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
    MSG("Usage: %s <lm or gn> <num frames> <num trials> <max num pts>\n",
        argv[0]);
    return 1;
  }
  int num_frames = atoi(argv[2]);
  int num_trials = atoi(argv[3]);
  int max_num_pts = atoi(argv[4]);

  std::string lm_str("lm");
  auto clear_func = lm_str == argv[1] ? lm<float> : gaussNewton<float>;
  auto secure_func = lm_str == argv[1] ? BuildLM : BuildGaussNewton;
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
        int num_pts = objectPointsSubset.size();

        if (!silent) {
          std::cout << "testing secure\n";
        }

        std::thread bob([&]() {
          NetIO* io = new NetIO("127.0.0.1", port);
          setup_semi_honest(io, ALICE);

          rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
          tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
          tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0
          auto [converged, num_loc_iterations] = emp_localize_wrapper(
              ALICE, rvec, tvec, cameraMatrix, distCoeffs, objectPointsSubset,
              imagePointsSubset, res, secure_func);
          if (!silent) {
            MSG("EMP reported %s after %d iterations.\n",
                converged ? "convergence" : "failure to converge",
                num_loc_iterations);
          }
          delete io;
        });

        NetIO* io = new NetIO(NULL, port);
        setup_semi_honest(io, BOB);
        rvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
        tvec = cv::Mat::zeros(3, 1, cv::DataType<float>::type);
        tvec.at<float>(2) = 1;  // cleartext has bug where it fails if z=0

        time_point<high_resolution_clock> tic = high_resolution_clock::now();

        auto [converged, num_loc_iterations] = emp_localize_wrapper(
            BOB, rvec, tvec, cameraMatrix, distCoeffs, objectPointsSubset,
            imagePointsSubset, res, secure_func);

        bob.join();

        time_point<high_resolution_clock> toc = high_resolution_clock::now();

        if (!silent) {
          cout << "secure result ";
          for (auto const& f : res)
            cout << f << ' ';
        }

        bool result_matches = true;
        for (int i = 0; i < 6; ++i) {
          result_matches &= withinRel(res[i], groundTruthRes[i]);
        }

        if (result_matches) {
          if (!silent) {
            cout << "emp converged!\n";
          }

          std::string timer_name = "emp_float_" + log_str + "_time_vs_points";
#if PPL_FLOW == PPL_FLOW_DO
          timer_name += "_dataobl";
#endif
          MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, timer_name.c_str(),
              num_pts,
              std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                      .count() /
                  1000000.0);
          MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__,
              (timer_name + "_per_loc_itr").c_str(), num_pts,
              std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                      .count() /
                  (1000000.0 * num_loc_iterations));

          MSG("SeNtInAl,grouped_bar,%s,%s%s,%d,%d\n", __FUNCTION__,
              log_str.c_str(), "_additions", num_pts, num_additions);
          MSG("SeNtInAl,grouped_bar,%s,%s%s,%d,%d\n", __FUNCTION__,
              log_str.c_str(), "_subtractions", num_pts, num_subtractions);
          MSG("SeNtInAl,grouped_bar,%s,%s%s,%d,%d\n", __FUNCTION__,
              log_str.c_str(), "_multiplications", num_pts,
              num_multiplications);
          MSG("SeNtInAl,grouped_bar,%s,%s%s,%d,%d\n", __FUNCTION__,
              log_str.c_str(), "_divisions", num_pts, num_divisions);

          MSG("SeNtInAl,xy,%s,%s%s,%d,%lu\n", __FUNCTION__, log_str.c_str(),
              "_bytes_tx", num_pts, io->size_tx);
          MSG("SeNtInAl,xy,%s,%s%s,%d,%lu\n", __FUNCTION__, log_str.c_str(),
              "_bytes_tx_per_itr", num_pts, io->size_tx / num_loc_iterations);
          MSG("SeNtInAl,xy,%s,%s%s,%d,%lu\n", __FUNCTION__, log_str.c_str(),
              "_bytes_tx_per_itr_per_feat", num_pts,
              io->size_tx / num_loc_iterations / num_pts);

          MSG("SeNtInAl,xy,%s,%s%s,%d,%lu\n", __FUNCTION__, log_str.c_str(),
              "_bytes_rx", num_pts, io->size_rx);
          MSG("SeNtInAl,xy,%s,%s%s,%d,%lu\n", __FUNCTION__, log_str.c_str(),
              "_bytes_rx_per_itr", num_pts, io->size_rx / num_loc_iterations);
          MSG("SeNtInAl,xy,%s,%s%s,%d,%lu\n", __FUNCTION__, log_str.c_str(),
              "_bytes_rx_per_itr_per_feat", num_pts,
              io->size_rx / num_loc_iterations / num_pts);

        } else {
          if (!silent) {
            MSG("Not printing timing - not converged.\n");
          }
        }

        // reset stats
        num_additions = 0;
        num_subtractions = 0;
        num_multiplications = 0;
        num_divisions = 0;
        delete io;

        return result_matches;
      },
      silent);
}
