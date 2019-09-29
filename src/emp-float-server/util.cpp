#include <stdio.h>
#include <iostream>
#include <vector>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

#include <jlog.h>
#include "util.h"

using namespace emp;
using namespace std;

#define BUILD_TIMING 0
#define EXEC_TIMING 0

void print_float32_bits(Float a, bool mute) {
  for (int i = 31; i >= 0; i--) {
    if (mute) {
      a[i].reveal<bool>();
    } else {
      printf("%d", a[i].reveal<bool>());
    }
  }
  if (!mute) {
    printf("\n\n");
  }
}

// prints an arbitrary size vector to the standard output
void printFloatVector(Float* v, int size, bool mute) {
  int i;

  for (i = 0; i < size; i++) {
    if (mute) {
      v[i].reveal<double>();
    } else {
      printf("%.8lf ", v[i].reveal<double>(PUBLIC));
    }
  }
  if (!mute) {
    printf("\n\n");
  }
}

// prints an arbitrary size matrix to the standard output
void printFloatMatrix(Float** a, int rows, int cols, bool mute) {
  int i, j;

  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      if (mute) {
        a[i][j].reveal<double>();
      } else {
        printf("%.8lf ", a[i][j].reveal<double>(PUBLIC));
      }
    }
    if (!mute) {
      printf("\n");
    }
  }
  printf("\n");
}
