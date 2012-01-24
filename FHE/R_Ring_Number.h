/*
 *  R_Ring_Number.h
 *  Created by valerini on 1/5/12.
 */
#ifndef _R_RING_NUMBER_H_
#define _R_RING_NUMBER_H_
#include <cassert>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

class R_Ring_Number {
  int q, d;
  int* vec;
public:
  ~R_Ring_Number() {
    if (vec != NULL) {
      delete [] vec;
    }
  }
  
  R_Ring_Number() {
    q = d = 0;
    vec = NULL;
  }
  
  R_Ring_Number(int q_, int d_) {
    Initialize(q_, d_);
  }
  
  R_Ring_Number(int q_, int d_, int array[]) {
    Initialize(q_, d_);
    for (int i = 0; i < d; i++) {
      vec[i] = array[i];
    }
  }
  
  R_Ring_Number(const R_Ring_Number &v) {
    Initialize(v.q, v.d);
    for (int i = 0; i < d; i++) {
      vec[i] = v.vec[i];
    }
  }
  
  // Assume Initialize is called on blank instance, otherwise there will be memory leak for vec
  void Initialize(int q_, int d_) {
    assert(q_ != 0 && d_ != 0);
    q = q_;
    d = d_;
    vec = new int [2 * d];
    for (int i = 0; i < 2 * d; i++) {
      vec[i] = 0;
    }
  }

  int Fast_Reduce(int n) const {
    int q_half_b = q / 2;
    int q_half_a = -(q - 1) / 2;
    if (n > q_half_b) {
      n -= q;
    } else if (n < q_half_a) {
      n += q;
    }
    assert(n <= q_half_b && n >= q_half_a);
    return n;
  }    

  static int Reduce(int n, int modul) {
    // for odd module reduce to range [-(q - 1) / 2, (q - 1) / 2]
    // for even module reduce to range [(q - 2) / 2, q / 2]
    int q_half_b = modul / 2;
    int q_half_a = -(modul - 1) / 2;
    if (n > q_half_b) {
      n -= ((n - q_half_b - 1) / modul + 1) * modul;
    } else if (n < q_half_a) {
      n += (-(n - q_half_a + 1) / modul + 1) * modul;
    }
    assert(n <= q_half_b && n >= q_half_a);
    return n;
  }    

  int Reduce(int n) const {
    return Reduce(n, q);
  }

  const R_Ring_Number operator -() const {
    R_Ring_Number res_v(q, d);
    for (int i = 0; i < d; i++) {
      // res_v.vec[i] = (q - vec[i]) % q;
      res_v.vec[i] = Fast_Reduce(-vec[i]);
    }
    return res_v;
  }
  
  const R_Ring_Number operator +(const R_Ring_Number &v) const {
    assert(q == v.q && d == v.d);
    R_Ring_Number res_v(q, d);
    for (int i = 0; i < d; i++) {
      // res_v.vec[i] = (vec[i] + v.vec[i]) % q;
      res_v.vec[i] = Fast_Reduce(vec[i] + v.vec[i]);
    }
    return res_v;
  }

  const R_Ring_Number operator -(const R_Ring_Number &v) const {
    assert(q == v.q && d == v.d);
    R_Ring_Number res_v(q, d);
    for (int i = 0; i < d; i++) {
      // res_v.vec[i] = (vec[i] + v.vec[i]) % q;
      res_v.vec[i] = Fast_Reduce(vec[i] - v.vec[i]);
    }
    return res_v;
  }
  
  R_Ring_Number& operator +=(const R_Ring_Number &v) {
    assert(d == v.d && q == v.q);
    for (int i = 0; i < d; i++) {
      // vec[i] = (vec[i] + v.vec[i]) % q;
      vec[i] = Fast_Reduce(vec[i] + v.vec[i]);
    }
    return *this;
  }
  
  R_Ring_Number& operator =(const R_Ring_Number &v) {
    if (d != v.d) {
      if (vec != NULL) {
	delete [] vec;
      }
      vec = new int [2 * v.d];
      d = v.d;
    }
    for (int i = 0; i < d; i++) {
      vec[i] = v.vec[i];
      if (q != 0 && q != v.q) {
	vec[i] = Reduce(vec[i]);
      }
    }
    for (int i = d; i < 2 * d; i++) {
      vec[i] = 0;
    }
    if (q == 0) q = v.q;
    return *this;
  }
  
  R_Ring_Number& operator =(int n) {
    assert(d != 0 && q > 0 && vec != NULL);
    for (int i = 0; i < d; i++) {
      vec[i] = 0;
    }
    vec[0] = Reduce(n);
    return *this;
  }
  
  R_Ring_Number operator *(int constant) const {
    R_Ring_Number res(q, d);
    for (int i = 0; i < d; i++) {
      res.vec[i] = Reduce(vec[i] * constant);
    }
    return res;
  }
  
  R_Ring_Number operator *(const R_Ring_Number &v) const {
    assert(q == v.q && d == v.d);
    R_Ring_Number res_v(q, d);
    
    for (int i = 0; i < d; i++) {
      if (v.vec[i] != 0) {
	for (int j = 0; j < d; j++) {
	  res_v.vec[i + j] += vec[j] * v.vec[i];
	}
      }
    }
    
    for (int i = 0; i < 2 * d; i++) {
      // res_v.vec[i] %= q;
      res_v.vec[i] = Reduce(res_v.vec[i]);
    }
    
    // reduction by polynomial x^d + 1
    for (int i = 2 * d - 2; i >= d; i--) {
      if (res_v[i] != 0) {
	res_v.vec[i - d] = Fast_Reduce(res_v.vec[i - d] - res_v.vec[i]);
	res_v.vec[i] = 0;
      }
    }
    return res_v;
  }
  
  R_Ring_Number operator -() {
    R_Ring_Number result(Get_q(), Get_d());
    for (int i = 0; i < d; i++) {
      // result.vec[i] = (-vec[i] + q) % q;
      result.vec[i] = Fast_Reduce(-vec[i]); // in fact reduction here is unnessesary
    }
    return result;
  }

  bool operator ==(const R_Ring_Number &r) {
    if (d != r.d) {
      return false;
    }
    for (int i = 0; i < d; i++) {
      if (vec[i] != r.vec[i]) {
	return false;
      }
    }
    return true;
  }

  bool operator !=(const R_Ring_Number &r) {
    return !((*this) == r);
  }
  
  int& operator [] (int index) {
    return vec[index];
  }
  
  R_Ring_Number Clamp(int modul) {
    for (int i = 0; i < d; i++) {
      vec[i] = Reduce(vec[i], modul);
    }
    return *this;
  }

  R_Ring_Number Get_Clamped(int modul) {
    R_Ring_Number result(*this);
    result.Clamp(modul);
    return result;
  }

  void Increase_Modul(int new_q) {
    assert(q <= new_q);
    q = new_q;
  }

  R_Ring_Number Scale(int q_, int p, int r) {
    assert(q_ == q);
    assert(p < q);
    assert(p % 2 == 1 && q % 2 == 1);
    assert(r == 2); // for simplicity, in future should be assert(r < p)
    float fraq = (p - 1) / (float)(q - 1); // (q - 1) / 2 should become (p - 1) / 2
    R_Ring_Number res_v(p, d);
    for (int i = 0; i < d; i++) {
      int desired_module = Reduce(vec[i], r);
      float tmp = vec[i] * fraq;
      int value[] = {tmp, tmp + 1, tmp - 1, tmp + 2, tmp - 2}; // TODO: think about better approach
      int max_dist = 2 * q;
      for (int j = 0; j < 5; j++) {
	value[j] = Reduce(value[j], p);
	int dist = fabs(tmp - value[j]);
	if (Reduce(value[j], r) == desired_module && dist < max_dist) {
	  max_dist = dist;
	  res_v[i] = value[j];
	}
      }
    }
    return res_v;
  }

  void print(void) {
    std::cout << "(";
    int i = 2 * d - 1;
    while (i >= 0 && vec[i] == 0) {
      i--;
    }
    if (i == -1) {
      std::cout << "0)";
    }
    
    for (; i >= 0; i--) {
      std::cout << vec[i];
      if (i == 0) {
	std::cout << ")";
      } else {
	std::cout << ", ";
      }
    }
  }
    
  /*  friend std::ostream & operator<< (const std::ostream &cout, const R_Ring_Number &v) const {
    cout << "(";
    int i = 2 * Get_d();
    while (i >= 0 && vec[i] == 0) {
      i--;
    }
    if (i == -1) {
      cout << "0)";
      return cout;
    }
    
    for (; i >= 0; i--) {
      cout << vec[i];
      if (i == 0) {
	cout << ")";
      } else {
	cout << ", ";
      }
    }
    return cout;
    }*/

  static R_Ring_Number Uniform_Rand(int q, int d, int bound = -1) {
    // ALLERT!!! use Salso20 PRG - it is believed to be PRG, it is based on some hard problem
    //srand(time(NULL));
    if (bound != -1) {
      bound = bound > q ? q : bound;
    } else {
      bound = q;
    }
    R_Ring_Number r(q, d);
    for (int i = 0; i < d; i++) {
      r.vec[i] = Reduce(rand(), bound);
    }
    return r;
  }
  
  int Get_d(void) const {
    return d;
  }
  
  int Get_q(void) const {
    return q;
  }

  int* Get_vec(void) const {
    return vec;
  }
};


#endif /* _R_RING_NUMBER_H_ */
