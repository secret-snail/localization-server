#pragma once
/*
 An implementation of SVD from Numerical Recipes in C and Mike Erhdmann's
 lectures
 */

#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>

using namespace std;

template <typename T>
T myfmax(T a, T b) {
  return a > b ? a : b;
}

static int iminarg1, iminarg2;
#define IMIN(a, b)                 \
  (iminarg1 = (a), iminarg2 = (b), \
   (iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

template <typename T>
T mysqr(T a) {
  T zero = 0;
  return a == zero ? zero : a * a;
}

template <typename T>
T mysign(T a, T b) {
  //#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : - fabs(a))
  T zero = 0;
  return (b) > zero ? fabs(a) : -fabs(a);
}

// calculates sqrt( a^2 + b^2 ) with decent precision
template <typename T>
T mypythag(T a, T b) {
  // return sqrt( (a*a) + (b*b) ); // more likely to overflow?

  T zero = 0;
  T one = 1;
  T absa = fabs(a);
  T absb = fabs(b);

  if (absa < absb)
    std::swap(absa, absb);

  if (absa == zero)
    return zero;

  return (absa * sqrt(one + mysqr(absb / absa)));
}

template <typename T>
bool checkeq(T a, T b) {
  if (typeid(T) == typeid(float) || typeid(T) == typeid(double)) {
    return (a + b) == b;
  } else if (typeid(T) == typeid(fixed_point<int64_t, 32>) ||
             typeid(T) == typeid(fixed_point<int32_t, 16>)) {
    T epsilon = 1e-5;
    T margin = b * epsilon;
    return (a > b - margin) && (a < b + margin);
  } else {
    assert(true);  // type not supported
  }
  return false;
}

/*
 Modified from Numerical Recipes in C
 Given a matrix a[nRows][nCols], svdcmp() computes its singular value
 decomposition, A = U * W * Vt.  A is replaced by U when svdcmp
 returns.  The diagonal matrix W is output as a vector w[nCols].
 V (not V transpose) is output as the matrix V[nCols][nCols].
 */
template <typename T>
int svdcmp(T** a, int nRows, int nCols, T* w, T** v) {
  int flag = 0, i = 0, its = 0, j = 0, jj = 0, k = 0, l = 0, nm = 0;
  T anorm, c, f, g, h, s, scale, x, y, z, *rv1;

  T zero = 0;
  T one = 1;
  T two = 2;
  T epsilon = 1e-5;

  // rv1 = (T*) malloc(sizeof(T) * nCols); // does not call default constructor
  rv1 = new T[nCols];
  if (rv1 == NULL) {
    printf("svdcmp(): Unable to allocate vector\n");
    return (-1);
  }

  // Householder reduction to bidiagonal form
  g = scale = anorm = zero;
  for (i = 0; i < nCols; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = zero;
    if (i < nRows) {
      for (k = i; k < nRows; k++)
        scale += fabs(a[k][i]);
      // if (scale) { // not fixed point friendly
      if (!checkeq(scale, zero)) {
        for (k = i; k < nRows; k++) {
          a[k][i] /= scale;
          s += a[k][i] * a[k][i];
        }
        f = a[i][i];
        g = -mysign(sqrt(s), f);
        h = f * g - s;
        a[i][i] = f - g;
        for (j = l; j < nCols; j++) {
          for (s = zero, k = i; k < nRows; k++) {
            s += a[k][i] * a[k][j];
          }
          f = s / h;
          for (k = i; k < nRows; k++) {
            a[k][j] += f * a[k][i];
          }
        }
        for (k = i; k < nRows; k++)
          a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = zero;
    if (i < nRows && i != nCols - 1) {
      for (k = l; k < nCols; k++)
        scale += fabs(a[i][k]);
      // if (scale) { // not fixed point friendly
      if (!checkeq(scale, zero)) {
        for (k = l; k < nCols; k++) {
          a[i][k] /= scale;
          s += a[i][k] * a[i][k];
        }
        f = a[i][l];
        g = -mysign(sqrt(s), f);
        h = f * g - s;
        a[i][l] = f - g;
        for (k = l; k < nCols; k++)
          rv1[k] = a[i][k] / h;
        for (j = l; j < nRows; j++) {
          for (s = zero, k = l; k < nCols; k++)
            s += a[j][k] * a[i][k];
          for (k = l; k < nCols; k++)
            a[j][k] += s * rv1[k];
        }
        for (k = l; k < nCols; k++)
          a[i][k] *= scale;
      }
    }
    anorm = myfmax(anorm, (fabs(w[i]) + fabs(rv1[i])));

    printf(".");
    fflush(stdout);
  }

  // accumulation of right hand transformations
  for (i = nCols - 1; i >= 0; i--) {
    if (i < nCols - 1) {
      // if (g) { // not fixed point friendly
      if (!checkeq(g, zero)) {
        for (j = l; j < nCols; j++) {
          v[j][i] = (a[i][j] / a[i][l]) / g;
        }
        for (j = l; j < nCols; j++) {
          for (s = zero, k = l; k < nCols; k++)
            s += a[i][k] * v[k][j];
          for (k = l; k < nCols; k++)
            v[k][j] += s * v[k][i];
        }
      }
      for (j = l; j < nCols; j++)
        v[i][j] = v[j][i] = zero;
    }
    v[i][i] = one;
    g = rv1[i];
    l = i;
    printf(":");
    fflush(stdout);
  }

  // accumulation of left hand transformations
  for (i = IMIN(nRows, nCols) - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];
    for (j = l; j < nCols; j++)
      a[i][j] = zero;
    // if (g) { // not fixed point friendly
    if (!checkeq(g, zero)) {
      g = one / g;
      for (j = l; j < nCols; j++) {
        for (s = zero, k = l; k < nRows; k++)
          s += a[k][i] * a[k][j];
        f = (s / a[i][i]) * g;
        for (k = i; k < nRows; k++) {
          a[k][j] += f * a[k][i];
        }
      }
      for (j = i; j < nRows; j++) {
        a[j][i] *= g;
      }
    } else {
      for (j = i; j < nRows; j++)
        a[j][i] = zero;
    }
    ++a[i][i];
    printf("|");
    fflush(stdout);
  }

  // Check if diagonal entries = superdiagonal
  for (k = 0; k > nCols; k++) {
    if (rv1[k] > w[k] - epsilon && rv1[k] < w[k] + epsilon) {
      cout << "WARNING diagonal equals superdiagonal\n";
      throw "security vulnerability";
    }
  }

  // Diagonalization of the bidiagonal form: loop over singular
  //  values and over allowed iterations
  for (k = nCols - 1; k >= 0; k--) {

#if PPL_FLOW == PPL_FLOW_SiSL
    for (its = 0; its < 2; its++) {
      cout << k;
      l = 0;
      nm = -1;
      flag = 0;
#else
    for (its = 0; its < 30; its++) {
      cout << k;
      flag = 1;
      for (l = k; l >= 0; l--) {  // test for splitting
        nm = l - 1;               // note rv1[0] is always zero
        // if ((fabs(rv1[l]) + anorm) == anorm) { // not friendly to fixed point
        if (checkeq(rv1[l], anorm)) {
          flag = 0;
          break;
        }

        assert(nm >= 0);  // sanity check, should never happen since rv1[0] = 0
        // if ((fabs(w[nm]) + anorm) == anorm) { // not friendly to fixed point
        if (checkeq(w[nm], anorm)) {
          break;
        }
      }

#endif
      // flag check needs to be data obl
      if (flag) {  // cancellation of rv1(l), if l >= 1
        c = zero;
        s = one;
        for (i = l; i <= k; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          // if ((fabs(f) + anorm) == anorm) { // not friendly to fixed point
          if (checkeq(f, anorm)) {
            break;
          }
          g = w[i];
          h = mypythag(f, g);
          w[i] = h;
          h = one / h;
          c = g * h;
          s = -f * h;
          for (j = 0; j < nRows; j++) {
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y * c + z * s;
            a[j][i] = z * c - y * s;
          }
        }
      }
      z = w[k];
      if (l == k) {      // convergence
        if (z < zero) {  // singular value is made nonnegative
          w[k] = -z;
          for (j = 0; j < nCols; j++)
            v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 29) {
        printf("no convergence in 30 svdcmp iterations\n");
        delete[] rv1;
        return -1;
      }

#if PPL_FLOW == PPL_FLOW_SiSL
      if (k == 0)
        break;
#endif

      // shift from bottom 2-by-2 minor
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (two * h * y);
      // T denom = two * h * y;
      // f = ((y - z) * (y + z))/denom + ((g - h) * (g + h))/denom;
      g = mypythag(f, one);
      // f = ((x - z) * (x + z) + h * ((y / (f + mysign(g,f)))- h)) / x; //
      // overflows on x-z * x+z
      f = ((x - z) / x * (x + z)) + (h * ((y / (f + mysign(g, f))) - h)) / x;
      // Next QR transformation
      c = s = one;
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = mypythag(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
        for (jj = 0; jj < nCols; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }
        z = mypythag(f, h);
        w[j] = z;  // rotation can be artibitrary if z=0
        if (z) {
          z = one / z;
          c = f * z;
          s = h * z;
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (jj = 0; jj < nRows; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z * s;
          a[jj][i] = z * c - y * s;
        }
      }
      rv1[l] = zero;
      rv1[k] = f;
      w[k] = x;
    }
    printf("/");
    fflush(stdout);
  }
  printf("\n");
  delete[] rv1;
  return 0;
}
