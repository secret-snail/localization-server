#pragma once
// M mxn
// N mmxnn
// res mxnn
template <typename T>
T twonormsq(T* vect, int sz) {
  T sum = 0.0;
  for (int i = 0; i < sz; i++) {
    sum += vect[i] * vect[i];
  }
  return sum;
}
