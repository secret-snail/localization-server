#pragma once
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <vector>

constexpr const float ms = 0.65;

template <typename T>
void AprilSnail3DPoints(std::vector<cv::Point3_<T>>& points) {
  points.push_back(cv::Point3_<T>(-ms, -ms, 0));
  points.push_back(cv::Point3_<T>(-ms, ms, 0));
  points.push_back(cv::Point3_<T>(ms, ms, 0));
  points.push_back(cv::Point3_<T>(ms, -ms, 0));
}

template <typename T>
void AprilSnail2DPoints(std::vector<cv::Point_<T>>& points) {
  points.push_back(cv::Point_<T>(207.71281433, 159.10252380));
  points.push_back(cv::Point_<T>(230.89843750, 334.47018433));
  points.push_back(cv::Point_<T>(403.77438354, 314.49819946));
  points.push_back(cv::Point_<T>(384.17871094, 138.54327393));
}
