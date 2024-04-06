#pragma once
#include <stdint.h>

// Algorithm control flow
// DO = data oblivious running a maximal, fixed upper bound of iterations
// of gradient descent and svd.
#define PPL_FLOW_DO 2
// Loop Leak is where the number of iterations is revealed to offload servers.
#define PPL_FLOW_LOOP_LEAK 3
// Single shot localization (SiSL) runs a single optimization iteration.
// It also uses a data oblivious SVD.
#define PPL_FLOW_SiSL 4
// Fakes localization at the same rate as plaintext for power testing
// on the raspberry pi snail.
#define PPL_FLOW_POWER_TESTING 5

#define PPL_FLOW PPL_FLOW_SiSL

const float GT_MIN_ER =
    1;  // (pose - ground truth) L2 norm less than this considered correct

const float JACOB_EPSILON = 0x0.0000c3p0;  // 25 mpc iterations
const float MIN_ER = 1e-2;

// Another set of values that works for the ETH3D dataset.
// const float JACOB_EPSILON = 0x0.000050p0;
// const float MIN_ER = 5e-2;

// maximum number of iterations
// from opencv calibration.cpp - cvFindExtrinsicCameraParams2() = 20
const int GN_MAX_ITR = 30;
const int LM_MAX_ITR = 30;
const float LM_LAMBDA_INIT = 1e-3;
const float LM_LAMBDA_MAX = 1e5;  // opencv default, works with hoff and ETH3D
// const float LM_LAMBDA_MAX = 1e2;  // snail with aruco markers may require smaller lambda max
const float LM_LAMBDA_MIN = 1e-5;

// Debugging tools
#define DBG_FLOW 0x01
#define DBG_PROJECT 0x02
#define DBG_JACOB 0x04
#define DBG_POSE_UPDATE 0x08
#define DBG_ER 0x10
#define DBG_ARGS 0x20
#define DBG_LAMBDA 0x40
const uint32_t _verbosity = 0x00;
const bool printints =
    false;  // used for ABY, their debug prints all floats as ints
