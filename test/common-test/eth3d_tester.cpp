#include <libgen.h>
#include <linux/limits.h>
#include <stdio.h>
#include <string.h>

#include <string>
#include <vector>

#include <opencv2/opencv.hpp>

#include <eth3d_features.hpp>

using cv::Point2f;

int main(int argc, char** argv) {
  char datadir[PATH_MAX];
  strncpy(datadir, argv[0], sizeof(datadir));
  dirname(datadir);
  std::string bindir(datadir);

  std::string base_path = bindir + "/../../data-eth3d/";
  auto feats = ETH3DFeatures<float>(base_path, eth3d_locations[0]);

  std::string rel_image_path = feats.getRelImageFilePath();
  std::string image_path = base_path + "images/" + rel_image_path;

  std::cout << "First image directory:" << std::endl;
  std::cout << image_path << std::endl;

  std::cout << "Getting points..." << std::flush;
  std::vector<Point2f> image_points;
  std::vector<Point3f> world_points;
  feats.imageFeatures(image_points);
  feats.worldFeatures(world_points);
  auto res = feats.getGroundTruthPose();
  std::cout << " done\n";

  std::cout << "first point\n";
  std::cout << "2d " << image_points[0].x << " " << image_points[0].y
            << std::endl;
  std::cout << "3d " << world_points[0].x << " " << world_points[0].y << " "
            << world_points[0].z << std::endl;

  std::cout << "pose: ";
  for (int i = 0; i < 6; i++)
    std::cout << res[i] << " ";
  std::cout << std::endl;

  return 0;
}
