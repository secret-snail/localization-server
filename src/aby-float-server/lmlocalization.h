#include <iostream>
#include <vector>

#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <abycore/aby/abyparty.h>
#include <abycore/circuit/arithmeticcircuits.h>
#include <abycore/circuit/booleancircuits.h>
#include <abycore/circuit/circuit.h>
#include <abycore/sharing/sharing.h>

using namespace std;

void BuildAndRunLM(share* s_threeDPts[], share* s_twoDPts[], int numPts,
                   share* s_f, share* s_cx, share* s_cy, share* s_x[],
                   BooleanCircuit* c, ABYParty* party, e_role role);

uint32_t test_lm_circuit(
    e_role role, const std::string& address, uint16_t port, seclvl seclvl,
    uint32_t nthreads, e_mt_gen_alg mt_alg, e_sharing sharing,
    vector<cv::Point3f> threeDPts, vector<cv::Point2f> twoDPts, float f,
    float cx, float cy,
    float* x /* initial guess for { r1, r2, r3, t1, t2, t3 } */);
