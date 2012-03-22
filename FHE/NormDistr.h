#ifndef _NORMDISTR_H_
#define _NORMDISTR_H_

#include "R_Ring_Vector.h"
#include <math.h>

class NormDistr {
 public:
  static long double sample_standard(long double center, int deviation) {
    static bool spare_ready = false;
    static long double spare;
    long double S;

    if (spare_ready) {
      spare_ready = false;
      return spare;
    }

    long double U, V;

    do {
      U = (RandomBnd(RAND_MAX) / (long double)RAND_MAX - 0.5) * 2;
      V = (RandomBnd(RAND_MAX) / (long double)RAND_MAX - 0.5) * 2;
      //      std::cout << "(" << U << ", " << V << ")" << std::endl;

      S = sqrt(U * U + V * V);
    } while (S >= 1 || S == 0);

    S = sqrt(-2 * log(S) / S) * deviation;
    spare_ready = true;
    spare = center + V * S;

    return center + U * S;
  }
  
  static R_Ring_Number sample(ZZ q, int d, int deviation) {
    R_Ring_Number res(q, d);
    for (int i = 0; i < d; i++) {
      res[i] = sample_standard(0, deviation);
      //      std::cout << "res[" << i << "] = " << res[i] << std::endl;
      if (res[i] >= q) {
	res[i] = q - 1;
	//	std::cout << "ALERT! :" << res[i] << " >= " << q << std::endl;
	//	assert(res[i] < q);
      }
    }
    return res;
  }

  static R_Ring_Vector sample_vector(ZZ q, int d, int dimension, int deviation) {
    R_Ring_Vector v(q, d, dimension);

    for (int i = 0; i < dimension; i++) {
      v[i] = sample(q, d, deviation);
    }
    return v;
  }
};

#endif /* _NORMDISTR_H_ */
