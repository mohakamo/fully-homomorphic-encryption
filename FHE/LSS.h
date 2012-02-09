#ifndef _LSS_H_
#define _LSS_H_

#include <vector>
#include "FHE.h"

template<class T>
void Sum_Up(std::vector<T> &x) {
  int size = x.size();
  
  while (size > 1) {
    for (int i = 0; i < size / 2; i++) {
      x[i] = x[2 * i] + x[2 * i + 1];
    }
    if (size % 2 != 0) {
      x[size / 2] = x[size - 1];
    }
    size /= 2;
  }
}

template<class T>
Pair<T, T> Compute_LSS(std::vector<T> xs, std::vector<T> ys) {
  // TODO: need to add in a way to get logarithmic depth!!!!!
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

  if (typeid(T) == typeid(FHE_Cipher_Text)) {
    t2 = (reinterpret_cast<FHE_Cipher_Text *>(&sum_x[0]))->Get_Cipher().second;
    t3 = (reinterpret_cast<FHE_Cipher_Text *>(&sum_y[0]))->Get_Cipher().second;
  }

  T tmp1 = sum_x[0] * sum_y[0]; // bounded by (size * M)^2
  if (typeid(T) == typeid(FHE_Cipher_Text)) {
    t1 = (reinterpret_cast<FHE_Cipher_Text *>(&tmp1))->Get_Cipher().second + 1;
    if (t1 != t2) {
      std::cout << "#1 tmp1 = " << t1 << ", sum_x = " << t2 << ", sum_y = " << t3 << std::endl;
      exit(1);
    }
  }
  T tmp2 = x_y[0] * ZZ(INIT_VAL, xs.size()); // bounded by (size * M)^2
  if (typeid(T) == typeid(FHE_Cipher_Text)) {
    t1 = (reinterpret_cast<FHE_Cipher_Text *>(&tmp2))->Get_Cipher().second + 1;
    t2 = (reinterpret_cast<FHE_Cipher_Text *>(&x_y[0]))->Get_Cipher().second;
    if (t1 != t2) {
      std::cout << "#2 t1 = " << t1 << ", t2 = " << t2 << std::endl;
      exit(1);
    }
  }
  a = tmp2 - tmp1; // bounded by 2 * (size * M)^2
  T tmp3 = x_x[0] * sum_y[0]; // bounded by size^2 * M^3
  T tmp4 = sum_x[0] * x_y[0]; // bounded by size^2 * M^3
  b = tmp3 - tmp4; // bounded by 2 * size^2 * M^3

  return Pair<T, T>(a, b);
}

#endif /* _LSS_H_ */