#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

//#include <opencv2/opencv.hpp>
//#include <opencv2/core/core.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/calib3d/calib3d.hpp>
//#include <opencv2/highgui/highgui.hpp>

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

// r=3x1, R=3x4 (only 3x3 is used here)
// values are secret, sizes are public
void BuildRodriguesCircuit(share* r[], share* R[], BooleanCircuit* c);
