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
  long long q;
  int d;
  long long* vec;
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
  
  R_Ring_Number(long long q_, int d_) {
    Initialize(q_, d_);
  }
  
  R_Ring_Number(long long q_, int d_, long long array[]) {
    Initialize(q_, d_);
    for (int i = 0; i < d; i++) {
      vec[i] = array[i];
    }
  }
  
  R_Ring_Number(const R_Ring_Number &v) {
    if (v.q == 0 || v.d == 0) {
      q = d = 0;
      vec = NULL;
      return;
    }
    Initialize(v.q, v.d);
    for (int i = 0; i < d; i++) {
      vec[i] = v.vec[i];
    }
  }
  
  // Assume Initialize is called on blank instance, otherwise there will be memory leak for vec
  void Initialize(long long q_, int d_) {
    assert(q_ != 0 && d_ != 0);
    q = q_;
    d = d_;
    vec = new long long [2 * d];
    for (int i = 0; i < 2 * d; i++) {
      vec[i] = 0;
    }
  }

  long long Fast_Reduce(long long n) const {
    return Reduce(n);
    /*
    long long q_half_b = q / 2;
    long long q_half_a = -(q - 1) / 2;
    if (n > q_half_b) {
      n -= q;
    } else if (n < q_half_a) {
      n += q;
    }
    assert(n <= q_half_b && n >= q_half_a);
    return n; */
  }

  static long long Reduce(long long n, long long modul) {
    // for odd module reduce to range [-(q - 1) / 2, (q - 1) / 2]
    // for even module reduce to range [(q - 2) / 2, q / 2]
    long long q_half_b = modul / 2;
    long long q_half_a = -(modul - 1) / 2;
    if (n > q_half_b) {
      n -= ((n - q_half_b - 1) / modul + 1) * modul;
    } else if (n < q_half_a) {
      n += (-(n - q_half_a + 1) / modul + 1) * modul;
    }
    assert(n <= q_half_b && n >= q_half_a);
    return n;
  }    

  long long Reduce(long long n) const {
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
      vec = new long long [2 * v.d];
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
  
  R_Ring_Number& operator =(long long n) {
    assert(d != 0 && q > 0 && vec != NULL);
    for (int i = 0; i < d; i++) {
      vec[i] = 0;
    }
    vec[0] = Reduce(n);
    return *this;
  }
  
  R_Ring_Number operator *(long long constant) const {
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
	  res_v.vec[i + j] = Reduce(res_v.vec[i + j] + Reduce(vec[j] * v.vec[i]));
	}
      }
    }
    
    //    for (int i = 0; i < 2 * d; i++) {
      // res_v.vec[i] %= q;
    //      res_v.vec[i] = Reduce(res_v.vec[i]);
    //  }
    
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
  
  long long& operator [] (int index) {
    return vec[index];
  }

  static long long Clamp(long long n, long long modul) {
    // for odd module reduce to range [-(q - 1) / 2, (q - 1) / 2]
    // for even module reduce to range [(q - 2) / 2, q / 2]
    long long q_half_b = modul / 2;
    long long q_half_a = -(modul - 1) / 2;
    if (n > q_half_b) {
      n -= ((n - q_half_b - 1) / modul + 1) * modul;
    } else if (n < q_half_a) {
      n += (-(n - q_half_a + 1) / modul + 1) * modul;
    }
    assert(n <= q_half_b && n >= q_half_a);
    return n;
  }
  
  R_Ring_Number Clamp(long long modul) {
    for (int i = 0; i < d; i++) {
      vec[i] = Clamp(vec[i], modul);
    }
    return *this;
  }

  R_Ring_Number Get_Clamped(long long modul) {
    R_Ring_Number result(*this);
    result.Clamp(modul);
    return result;
  }

  void Change_Modul(long long new_q) {
    assert(new_q > 0);
    q = new_q;
    for (int i = 0; i < d; i++) {
      vec[i] = Reduce(vec[i]);
    }
  }

  void Increase_Modul(long long new_q) {
    assert(q <= new_q);
    q = new_q;
  }

  void Decrease_Modul(long long new_q) {
    assert(q >= new_q);
    for (int i = 0; i < d; i++) {
      vec[i] = Reduce(vec[i]);
    }
    q = new_q;
  }

  R_Ring_Number Scale(long long q_, long long p, int r) {
    if (r == 2) {
    assert(q_ == q);
    assert(p < q);
    assert(p % 2 == 1 && q % 2 == 1);
    assert(r == 2); // for simplicity, in future should be assert(r < p)
    long double fraq = (p - 1) / (long double)(q - 1); // (q - 1) / 2 should become (p - 1) / 2
    R_Ring_Number res_v(p, d);
    for (int i = 0; i < d; i++) {
      int desired_module = Clamp(vec[i], r);
      long double tmp = vec[i] * fraq;
      long long value[] = {tmp, tmp + ((tmp > 0) ? -1 : 1), tmp + ((tmp > 0) ? 1 : -1), tmp + 2, tmp - 2}; // TODO: think about better approach
      double max_dist = LLONG_MAX;
      for (int j = 0; j < 5; j++) {
	value[j] = Reduce(value[j], p);
	double dist = fabs(tmp - value[j]);
	if (Clamp(value[j], r) == desired_module && dist < max_dist) {
	  max_dist = dist;
	  res_v[i] = value[j];
	}
      }
      assert(fabs(res_v[i] - tmp) <= 0.5 * r);
      assert(max_dist != LLONG_MAX);
    }
    return res_v;
    } else {
    assert(q_ == q);
    assert(p < q);
    assert(p % 2 == 1 && q % 2 == 1);
    double fraq = (p - 1) / (double)(q - 1); // (q - 1) / 2 should become (p - 1) / 2
    R_Ring_Number res_v(p, d);
    for (int i = 0; i < d; i++) {
      int desired_module = Clamp(vec[i], r);
      long double tmp_d = (vec[i] / (double)(q - 1)) * (p - 1);
      long long tmp = tmp_d;
      tmp -= Clamp(tmp, r);
      tmp += desired_module;
      long long value[] = {tmp, tmp + (tmp > 0 ? -r : r), tmp + (tmp > 0 ? r : -r), tmp - 2 * r, tmp + 2 * r}; // TODO: think about better approach
      double max_dist = LLONG_MAX;
      for (int j = 0; j < 5; j++) {
	value[j] = Reduce(value[j], p);
	double dist = fabs(tmp_d - value[j]);
	if (dist < max_dist && Clamp(value[j], r) == desired_module) {
	  max_dist = dist;
	  res_v[i] = value[j];
	}
      }
      if (fabs(res_v[i] - tmp_d) > r + 1e-5) {
	std::cout << "q = " << q << ", p = " << p << ", tmp_d = " << tmp_d << ", vec[i] = " << vec[i] << ", tmp = " << tmp << ", res_v[i] = " << res_v[i] << ", r = " << r << ", fraq = " << fraq << std::endl;
	std::cout << "res_v[i] - tmp_d = " << res_v[i] - tmp_d << std::endl;
	std::cout << "value[0] - tmp_d = " << value[0] - tmp_d << std::endl;
	std::cout << "value[1] - tmp_d = " << value[1] - tmp_d << std::endl;
	std::cout << "value[2] - tmp_d = " << value[2] - tmp_d << std::endl;
	std::cout << "value[3] - tmp_d = " << value[3] - tmp_d << std::endl;
	std::cout << "value[4] - tmp_d = " << value[4] - tmp_d << std::endl;

	std::cout << "0.5 * r = " << 0.5 * r << std::endl;
	
	std::cout << "value[0] = " << value[0] << std::endl;
	std::cout << "value[1] = " << value[1] << std::endl;
	std::cout << "value[2] = " << value[2] << std::endl;
	std::cout << "value[3] = " << value[3] << std::endl;
	std::cout << "value[4] = " << value[4] << std::endl;

	std::cout << "sizeof(double) = " << sizeof(double) << std::endl;
	exit(1);
      }
      // assert(fabs(res_v[i] - tmp) <= 0.5 * r);
      assert(max_dist != LLONG_MAX);
    }
    return res_v;
    }      
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

  static R_Ring_Number Uniform_Rand(long long q, int d, long long bound = -1) {
    // ALLERT!!! use Salso20 PRG - it is believed to be PRG, it is based on some hard problem
    //srand(time(NULL));
    if (bound != -1) {
      bound = bound > q ? q : bound;
    } else {
      bound = q;
    }
    R_Ring_Number r(q, d);
    for (int i = 0; i < d; i++) {
      r.vec[i] = /*rand() % bound;*/Clamp(rand(), bound);
    }
    return r;
  }
  
  int Get_d(void) const {
    return d;
  }
  
  long long Get_q(void) const {
    return q;
  }

  long long* Get_vec(void) const {
    return vec;
  }

  double Get_Field_Expansion(void) const {
    return sqrt((double)Get_d());
  }
  
  long long Get_Norm(void) const {
    long long result = 0;
    for (int i = 0; i < d; i++) {
      result += abs(vec[i]);
    }
    return result;
  }
};


#endif /* _R_RING_NUMBER_H_ */

