#pragma once
#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

using namespace emp;

void BuildSignCircuit(Float* a, Float* b);

Float BuildPythagCircuit(Float* a, Float* b);

int BuildSvdCircuit(Float** a, int nRows, int nCols, Float* w, Float** v);
