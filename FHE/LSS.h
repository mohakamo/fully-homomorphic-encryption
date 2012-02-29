#ifndef _LSS_H_
#define _LSS_H_

#include <vector>
#include "FHE.h"

template<class T>
void Sum_Up(std::vector<T> &x) { // summing up in a way to get logarithmic depth
  int size = x.size();
  
  /*
  while (size > 1) {
    for (int i = 0; i < size / 2; i++) {
      x[i] = x[2 * i] + x[2 * i + 1];
    }
    if (size % 2 != 0) {
      x[size / 2] = x[size - 1];
    }
    size /= 2;
  }
  */
  int step = 1;

  while (step < size) {
    for (int i = 0; i < size; i += 2 * step) {
      if (i + step < size) {
	x[i] = x[i] + x[i + step]; // TODO: add operator +=
      }
    }
    step *= 2;
  }
}

template<class T>
Pair<T, T> Compute_LSS(std::vector<T> xs, std::vector<T> ys, T * den = NULL) {
  int size = xs.size();
  assert(xs.size() == ys.size());
  std::vector<T> x_x(size), x_y(size), sum_x(size), sum_y(size);
  // if xs(i) and ys(i) are bounded by say M, then...
  for (int i = 0; i < size; i++) {
    x_x[i] = xs[i] * xs[i]; // bounded by M^2
    x_y[i] = xs[i] * ys[i]; // M^2
    sum_x[i] = xs[i];
    sum_y[i] = ys[i];
  }

  Sum_Up<T>(x_x); // bounded by size * M^2
  Sum_Up<T>(x_y); // M^2
  Sum_Up<T>(sum_x); // bounded by size * M
  Sum_Up<T>(sum_y); // size * M

  T a, b;
  int t1, t2, t3;

  if (den != NULL) {
    T tmp_den1 = x_x[0] * ZZ(INIT_VAL, xs.size());
    T tmp_den2 = sum_x[0] * sum_x[0];
    (*den) = tmp_den1 - tmp_den2;
  }


  T tmp1 = sum_x[0] * sum_y[0]; // bounded by (size * M)^2
  T tmp2 = x_y[0] * ZZ(INIT_VAL, xs.size()); // bounded by (size * M)^2
  a = tmp2 - tmp1; // bounded by 2 * (size * M)^2

  T tmp3 = x_x[0] * sum_y[0]; // bounded by size^2 * M^3
  T tmp4 = sum_x[0] * x_y[0]; // bounded by size^2 * M^3
  b = tmp3 - tmp4; // bounded by 2 * size^2 * M^3

  return Pair<T, T>(a, b);
}

#endif /* _LSS_H_ */
