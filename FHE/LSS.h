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
  for (int i = 0; i < size; i++) {
    x_x[i] = xs[i] * xs[i];
    x_y[i] = xs[i] * ys[i];
    sum_x[i] = xs[i];
    sum_y[i] = ys[i];
  }

  Sum_Up<T>(x_x);
  Sum_Up<T>(x_y);
  Sum_Up<T>(sum_x);
  Sum_Up<T>(sum_y);

  T a, b;
  int t1, t2, t3;

  if (typeid(T) == typeid(FHE_Cipher_Text)) {
    t2 = (reinterpret_cast<FHE_Cipher_Text *>(&sum_x[0]))->Get_Cipher().second;
    t3 = (reinterpret_cast<FHE_Cipher_Text *>(&sum_y[0]))->Get_Cipher().second;
  }

  T tmp1 = sum_x[0] * sum_y[0];
  if (typeid(T) == typeid(FHE_Cipher_Text)) {
    t1 = (reinterpret_cast<FHE_Cipher_Text *>(&tmp1))->Get_Cipher().second + 1;
    if (t1 != t2) {
      std::cout << "#1 tmp1 = " << t1 << ", sum_x = " << t2 << ", sum_y = " << t3 << std::endl;
      exit(1);
    }
  }
  T tmp2 = x_y[0] * xs.size();
  if (typeid(T) == typeid(FHE_Cipher_Text)) {
    t1 = (reinterpret_cast<FHE_Cipher_Text *>(&tmp2))->Get_Cipher().second + 1;
    t2 = (reinterpret_cast<FHE_Cipher_Text *>(&x_y[0]))->Get_Cipher().second;
    if (t1 != t2) {
      std::cout << "#2 t1 = " << t1 << ", t2 = " << t2 << std::endl;
      exit(1);
    }
  }
  a = tmp2 - tmp1;
  T tmp3 = x_x[0] * sum_y[0];
  T tmp4 = sum_x[0] * x_y[0];
  b = tmp3 - tmp4;

  return Pair<T, T>(a, b);
}

#endif /* _LSS_H_ */
