#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

using namespace emp;

// M mxn (M is secret, m n are public)
// N mmxnn
// res mxnn
void BuildMatmultCircuit(Float M[], int m, int n, Float N[], int mm, int nn,
                         Float res[]);

void BuildMatmult2DwTransposeCircuit(Float** M, int m, int n,
                                     const bool transposeM, Float** N, int mm,
                                     int nn, const bool transposeN,
                                     Float** res);
