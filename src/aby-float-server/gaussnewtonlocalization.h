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

void BuildGaussNewtonCircuit(share* threeDPts[], int numThreeD,
                             share* twoDPts[], int numTwoD, share* f, share* cx,
                             share* cy, share* x[], BooleanCircuit* c);

void BuildAndRunGaussNewton(share* s_threeDPts[], share* s_twoDPts[],
                            int numPts, share* s_f, share* s_cx, share* s_cy,
                            share* s_x[], BooleanCircuit* c, ABYParty* party,
                            e_role role);
