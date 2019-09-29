//#include <opencv2/opencv.hpp>
#include <jlog.h>
#include <opencv2/core/dualquaternion.hpp>

#include <libgen.h>
#include <linux/limits.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <string>
#include <vector>

using cv::DualQuatd;
using cv::Mat;
using cv::Point2f;
using cv::Point3_;
using cv::Point3f;
using cv::Point_;
using cv::Quatd;
using cv::Vec3f;

// vector<string> eth3d_indoor_locations = {
////"courtyard/",
//"delivery_area/",
////"electro/",
////"facade/",
//"kicker/",
////"meadow/",
//"office/",
//"pipes/",
////"playground/",
//"relief/",
//"relief_2/",
////"terrace/",
//"terrains/",
//};
//
// vector<string> eth3d_outdoor_locations = {
//"courtyard/",
////"delivery_area/",
//"electro/",
//"facade/",
////"kicker/",
//"meadow/",
////"office/",
////"pipes/",
//"playground/",
////"relief/",
////"relief_2/",
//"terrace/",
////"terrains/",
//};

std::vector<std::string> eth3d_locations = {
    "courtyard", "delivery_area", "electro",  "facade",     "kicker",
    "meadow",    "office",        "pipes",    "playground", "relief",
    "relief_2",  "terrace",       "terrains",
};

// typedef Eigen::Transform<float, 3, Eigen::Affine> SE3f;

struct ColmapFeatureObservation {
  // EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Sub-pixel coordinates of the observation in its image, given in pixels.
  // Eigen::Vector2f xy;
  Point2f xy;

  // Id of the corresponding 3D point or -1 if no 3D point is associated to this
  // observation.
  int point3d_id;
};

// Holds data of a COLMAP image.
struct ColmapImage {
  // Unique image id.
  int image_id;

  // Id of the camera model for this image.
  int camera_id;

  // Path to the image file, may be a relative path.
  std::string file_path;

  // Global-to-image transformation.
  // SE3f image_T_global;
  DualQuatd image_T_global;

  // Image-to-global transformation.
  // SE3f global_T_image;
  DualQuatd global_T_image;

  // 2D feature observations in this image.
  // std::vector<ColmapFeatureObservation,
  // Eigen::aligned_allocator<ColmapFeatureObservation>> observations;
  std::vector<ColmapFeatureObservation> observations;
};

typedef std::shared_ptr<ColmapImage> ColmapImagePtr;
typedef std::shared_ptr<const ColmapImage> ColmapImageConstPtr;

typedef std::vector<ColmapImagePtr> ColmapImagePtrVector;
typedef std::unordered_map<int, ColmapImagePtr> ColmapImagePtrMap;

template <class T>
class ETH3DFeatures {
 private:
  std::unordered_map<int, Point3_<T>> all_world_features;
  ColmapImagePtrMap images;
  bool error = false;
  std::vector<int> frame_to_id;
  int frame_number = 0;

 public:
  ETH3DFeatures(const std::string& base_path, const std::string& location) {
    const std::string calib_path = "/dslr_calibration_undistorted/";
    std::string images_txt_path =
        base_path + location + calib_path + "images.txt";
    std::string points3d_txt_path =
        base_path + location + calib_path + "points3D.txt";

    // open images.txt file
    std::ifstream images_file_stream(images_txt_path, std::ios::in);
    if (!images_file_stream) {
      ERROR_MSG("Could not open images.txt file\n");
      error = true;
      return;
    }

    // read all 2d points
    while (!images_file_stream.eof() && !images_file_stream.bad()) {
      std::string line;
      std::getline(images_file_stream, line);
      if (line.size() == 0 || line[0] == '#')
        continue;

      // Read image info line.
      ColmapImage* new_image = new ColmapImage();
      // Eigen::Quaternionf image_R_global;
      float qw, qx, qy, qz, tx, ty, tz;
      std::istringstream image_stream(line);
      image_stream >> new_image->image_id >> qw >> qx >> qy >> qz >> tx >> ty >>
          tz >> new_image->camera_id >> new_image->file_path;

      Quatd iTg_rot(qw, qx, qy, qz);
      // new_image->image_T_global.linear() = image_R_global.toRotationMatrix();
      new_image->image_T_global = DualQuatd::createFromAngleAxisTrans(
          iTg_rot
              .norm(),  // angle is the norm of rotation vector or quaternion?
          iTg_rot.toRotVec(),  // rvec describes rotation axis
          {tx, ty, tz}         // transation vector is easy
      );
      // new_image->global_T_image = new_image->image_T_global.inverse();
      new_image->global_T_image = new_image->image_T_global.inv();

      // Read feature observations line.
      std::getline(images_file_stream, line);
      std::istringstream observations_stream(line);
      while (!observations_stream.eof() && !observations_stream.bad()) {
        new_image->observations.emplace_back();
        ColmapFeatureObservation* new_observation =
            &new_image->observations.back();
        // observations_stream >> new_observation->xy.x()
        //                     >> new_observation->xy.y()
        //                     >> new_observation->point3d_id;
        observations_stream >> new_observation->xy.x >> new_observation->xy.y >>
            new_observation->point3d_id;
      }

      images.insert(
          std::make_pair(new_image->image_id, ColmapImagePtr(new_image)));
      frame_to_id.push_back(new_image->image_id);
    }

    // open points3d.txt file
    std::ifstream points3d_file_stream(points3d_txt_path, std::ios::in);
    if (!points3d_file_stream) {
      ERROR_MSG("Could not open points3d.txt file\n");
      error = true;
      return;
    }

    while (!points3d_file_stream.eof() && !points3d_file_stream.bad()) {
      std::string line;
      std::getline(points3d_file_stream, line);
      if (line.size() == 0 || line[0] == '#')
        continue;

      // Read 3d point line.
      int id;
      Point3_<T> p;
      std::istringstream image_stream(line);
      image_stream >> id >> p.x >> p.y >> p.z;
      //>> red
      //>> green
      //>> blue
      //>> reprojection error
      //>> in these 2d images
      all_world_features[id] = p;
    }
  }

  void imageFeatures(std::vector<Point_<T>>& points) {
    if (error) {
      ERROR_MSG("ETH3DFeatures not properly configured\n");
      return;
    }

    points.resize(0);
    int id = frame_to_id[frame_number];

    for (auto& obs : images[id]->observations) {
      if (obs.point3d_id == -1)
        continue;
      points.push_back({obs.xy.x, obs.xy.y});
    }
  }

  void worldFeatures(std::vector<Point3_<T>>& points) {
    if (error) {
      ERROR_MSG("ETH3DFeatures not properly configured\n");
      return;
    }

    int id = frame_to_id[frame_number];

    // for each 2d point find 3d correspondance
    for (auto& obs : images[id]->observations) {
      if (obs.point3d_id == -1)
        continue;  // no correspondance found by dataset owners

      const auto& f = all_world_features.find(obs.point3d_id);
      if (f == all_world_features.end()) {
        ERROR_MSG("2d-3d point correspondence could not be made\n");
        MSG("offending 2d point ID: %d\n", obs.point3d_id);

      } else {
        points.push_back({f->second.x, f->second.y, f->second.z});
      }
    }
  }

  std::vector<T> getGroundTruthPose() {
    if (error) {
      ERROR_MSG("ETH3DFeatures not properly configured\n");
      return {};
    }

    int id = frame_to_id[frame_number];
    auto m = images[id]->image_T_global;  // not global_T_image, right?
    auto r = m.getRotation().toRotVec();
    auto t = m.getTranslation();
    std::vector<T> pose = {static_cast<T>(r[0]), static_cast<T>(r[1]),
                           static_cast<T>(r[2]), static_cast<T>(t[0]),
                           static_cast<T>(t[1]), static_cast<T>(t[2])};
    return pose;
  }

  std::string getRelImageFilePath() {
    int id = frame_to_id[frame_number];
    return images[id]->file_path;
  }

  void nextFrame() {
    if (error) {
      ERROR_MSG("ETH3DFeatures not properly configured\n");
      return;
    }
    frame_number++;
  }
  int numberOfFrames() {
    if (error) {
      ERROR_MSG("ETH3DFeatures not properly configured\n");
      return 0;
    }
    return frame_to_id.size();
  }
};
