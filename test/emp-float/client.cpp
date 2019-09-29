//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <jlog.h>
#include <printutil.h>
#include <stdio.h>
#include <hoff_features.hpp>
#include <iostream>
#include <string>

#include <stdint.h>
#include <string.h>

#include <eth3d_features.hpp>
#include "emp-sh2pc/emp-sh2pc.h"

// to get path name
#include <libgen.h>
#include <linux/limits.h>
#include <stdio.h>
#include <string.h>

using namespace emp;
using namespace std;

PRG prg;
block delta;

void p128_hex_u32(__m128i in) {
  alignas(16) uint32_t v[4];
  _mm_store_si128((__m128i*)v, in);
  printf("v4_u32: %x %x %x %x\n", v[0], v[1], v[2], v[3]);
}

void sendFloatBlock(NetIO* io, float f) {
  int* in = (int*)&f;
  bool* b = new bool[32];
  int_to_bool<int>(b, *in, 32);

  block label[32];
  prg.random_block(label, 32);
  uint32_t mask = 1;
  for (int i = 0; i < 32; ++i) {
    if (b[i]) {
      label[i] = label[i] ^ delta;
    }
  }
  io->send_block(label, 32);
  // p128_hex_u32(label[0]);
  io->flush();
}

void recvFloatBlock(NetIO* aliceio, NetIO* bobio, float* f) {
  bool alicelsb = 0;
  bool boblsb = 0;
  uint32_t res = 0;
  for (int i = 0; i < 32; ++i) {
    aliceio->recv_data(&alicelsb, 1);
    bobio->recv_data(&boblsb, 1);
    uint32_t tmp = alicelsb ^ boblsb;
    res |= (tmp << i);
  }
  float* fp = (float*)(&res);
  *f = *fp;
}

int main(int argc, char** argv) {
  if (argc != 1) {
    MSG("Usage: %s\n", argv[0]);
    return 1;
  }
  MSG("Running\n");

  int baseport = 8080;
  NetIO* aliceio = new NetIO("127.0.0.1", baseport + ALICE * 17);
  NetIO* bobio = new NetIO("127.0.0.1", baseport + BOB * 17);
  MSG("Connected to alice port: %d, bob port: %d\n", baseport + ALICE * 17,
      baseport + BOB * 17);

  vector<cv::Point2f> imagePoints;
  vector<cv::Point3f> objectPoints;

  // float f=715;
  // float cx=354;
  // float cy=245;
  // Hoffs2DPoints(imagePoints);
  // Hoffs3DPoints(objectPoints);
  // vector<float> initialGuess = {1.5, -1.0, 0.0, 0, 0, 30};
  // uint32_t numPts = imagePoints.size();

  float f = 3408.57;
  float cx = 3114.7;
  float cy = 2070.92;
  char datadir[PATH_MAX];
  strncpy(datadir, argv[0], sizeof(datadir));
  dirname(datadir);
  std::string bindir(datadir);
  std::string base_path = bindir + "/../../data-eth3d/";
  auto feats = ETH3DFeatures<float>(base_path, eth3d_locations[0]);
  feats.imageFeatures(imagePoints);
  feats.worldFeatures(objectPoints);
  imagePoints.resize(6);
  objectPoints.resize(6);
  auto gtpose = feats.getGroundTruthPose();
  vector<float> initialGuess = {0, 0, 0, 0, 0, 0};
  // vector<float> initialGuess = {1, 1, -1, 0, 0, 4};
  // vector<float> initialGuess = {res.first[0], res.first[1], res.first[2],
  //         res.second[0], res.second[1], res.second[2]};
  uint32_t numPts = imagePoints.size();
  MSG("Found %d points\n", numPts);
  for (int i = 0; i < imagePoints.size(); i++) {
    MSG("Point %d 2d %f, %f\n", i, imagePoints[i].x, imagePoints[i].y);
    MSG("Point %d 3d %f, %f, %f\n", i, objectPoints[i].x, objectPoints[i].y,
        objectPoints[i].z);
  }
  printVector("Ground Truth Pose [rotation; translation]: ", &gtpose[0], 6);

  // Setup
  aliceio->send_data(&numPts, sizeof(uint32_t));
  bobio->send_data(&numPts, sizeof(uint32_t));
  aliceio->flush();
  bobio->flush();
  MSG("sent alice and bob numpoints\n");

  block seed;
  prg.random_block(&seed, 1);  // used for alice's inputs only
  prg.reseed(&seed);

  aliceio->recv_block(&delta, 1);
  MSG("got delta from alice\n");
  // p128_hex_u32(delta);

  // Distribute secret input data via 3-way OT
  // send prg seed directly to bob so he can generate his own labels
  bobio->send_block(&seed, 1);
  bobio->flush();
  MSG("sent bob seed\n");
  // p128_hex_u32(seed);

  // must send all labels ^ gc->delta to alice, she cannot learn seed
  for (int i = 0; i < numPts; i++) {
    sendFloatBlock(aliceio, objectPoints[i].x);
    sendFloatBlock(aliceio, objectPoints[i].y);
    sendFloatBlock(aliceio, objectPoints[i].z);

    sendFloatBlock(aliceio, imagePoints[i].x);
    sendFloatBlock(aliceio, imagePoints[i].y);
  }
  sendFloatBlock(aliceio, f);
  sendFloatBlock(aliceio, cx);
  sendFloatBlock(aliceio, cy);
  for (auto g : initialGuess) {
    sendFloatBlock(aliceio, g);
  }
  aliceio->flush();
  MSG("sent alice labels\n");

  // Offload servers run computation...
  MSG("Waiting for computation to finish...\n");

  // Collect result
  for (int i = 0; i < 6; i++) {
    recvFloatBlock(aliceio, bobio, &initialGuess[i]);
    cout << initialGuess[i] << " ";
  }
  cout << endl;

  return 0;
}
