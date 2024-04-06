#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

#include <trigfuncs.h>

using namespace emp;
using namespace std;

Float BuildSinCircuit(Float a) {
  Float pi = Float(M_PI, PUBLIC);
  return (a / pi).sin();  // circuit mutliplies by pi
}

Float BuildCosCircuit(Float a) {
  Float pi = Float(M_PI, PUBLIC);
  return (a / pi).cos();  // circuit mutliplies by pi
}
