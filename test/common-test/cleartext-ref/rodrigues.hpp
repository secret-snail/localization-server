#pragma once
#include <math.h>
#include <cleartext-ref/trigfuncs.hpp>

// r=3x1, R=3x4 (only 3x3 is used here)
template <typename T>
void rodrigues(T* r, T* R) {
  T theta = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

  T c = mycos<T>(theta);
  T s = mysin<T>(theta);
  T one = 1;
  T zero = 0;
  T c1 = one - c;
  T itheta = theta ? one / theta : zero;

  T x = r[0] * itheta;
  T y = r[1] * itheta;
  T z = r[2] * itheta;

  // R = cos(theta)*I + (1 - cos(theta))*r*rT + sin(theta)*[r_x]
  R[0] = c + c1 * x * x;
  R[1] = c1 * x * y - s * z;
  R[2] = c1 * x * z + s * y;

  R[4] = c1 * x * y + s * z;
  R[5] = c + c1 * y * y;
  R[6] = c1 * y * z - s * x;

  R[8] = c1 * x * z - s * y;
  R[9] = c1 * y * z + s * x;
  R[10] = c + c1 * z * z;
}
