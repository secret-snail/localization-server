#pragma once
#include <math.h>
#include <privacyconf.h>
#include <typeinfo>

template <typename T>
T mysin(T x) {
  // return x - (pow(x,3)/(3*2)) + (pow(x,5)/(5*4*3*2));// -
  // (pow(x,7)/(7*6*5*4*3*2));
  return sin(x);
}

template <typename T>
T mycos(T x) {
  // return 1 - (pow(x,2)/2) + (pow(x,4)/(4*3*2));// - (pow(x,6)/(6*5*4*3*2));
  return cos(x);
}
