#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

using namespace emp;
using namespace std;

// r=3x1, R=3x4 (only 3x3 is used here)
// values are secret, sizes are public
Float BuildSinCircuit(Float a);
Float BuildCosCircuit(Float a);
