#pragma once

#include <stdio.h>
#include <vector>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

using namespace emp;

void print_float32_bits(Float a, bool mute = false);

void printFloatMatrix(Float** a, int rows, int cols, bool mute = false);

void printFloatVector(Float* v, int size, bool mute = false);
