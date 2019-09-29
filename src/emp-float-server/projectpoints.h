#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

using namespace emp;
using namespace std;

// Note: All 2D matrix inputs passed as single dimension array
// P=4xnumPoints [x1, x2, ...;  y1, y2, ...;  z1, z2, ...; 1, 1, ... ]
// x=6x1 [ rotation angles; translation] (column)
// K=3x3 camera intrinsic
// result=3xnumPoints preallocated, only first two rows are x,y.
//      last row needed for intermediate calculation (homogeneous)
//      [ u1, u2, ...;  v1, v2, ...;  1, 1, ... ]
void BuildProjectPointsCircuit(Float _P[], Float _x[], Float _K[], Float _res[],
                               int numPoints, bool skipOneGate);
