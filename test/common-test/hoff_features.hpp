#pragma once
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <vector>

template <typename T>
void Hoffs3DPoints(std::vector<cv::Point3_<T>>& points) {
  points.push_back(cv::Point3_<T>(0, 10, 6));
  points.push_back(cv::Point3_<T>(0, 2, 6));
  points.push_back(cv::Point3_<T>(2, 0, 6));
  points.push_back(cv::Point3_<T>(0, 10, 2));
  points.push_back(cv::Point3_<T>(0, 2, 2));
  points.push_back(cv::Point3_<T>(2, 0, 2));
}

template <typename T>
void Hoffs2DPoints(std::vector<cv::Point_<T>>& points) {
  points.push_back(cv::Point_<T>(183, 147));
  points.push_back(cv::Point_<T>(350, 133));
  points.push_back(cv::Point_<T>(454, 144));
  points.push_back(cv::Point_<T>(176, 258));
  points.push_back(cv::Point_<T>(339, 275));
  points.push_back(cv::Point_<T>(444, 286));
}
