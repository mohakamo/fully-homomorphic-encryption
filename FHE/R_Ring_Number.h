/**
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
#include "NTL/ZZ.h"
#include "NTL/ZZX.h"

NTL_CLIENT

class R_Ring_Number {
  ZZ q;
  int d;

public:
  ZZX vec;

  R_Ring_Number() {
    q = d = 0;
  }
  
  R_Ring_Number(ZZ q_, int d_) {
    d = d_;
    q = q_;
    vec.SetMaxLength(2 * d);
    vec.rep.SetLength(d);
  }
  
  void Initialize(ZZ q_, int d_) {
    d = d_;
    q = q_;
    vec.SetMaxLength(2 * d);
    vec.rep.SetLength(d);
  }

  R_Ring_Number(ZZ q_, int d_, ZZ array[]) {
    q = q_;
    d = d_;
    vec_ZZ v;
    v.SetLength(d);
    for (int i = 0; i < d; i++) {
      v[i] = Reduce(array[i]);
    }
    vec = to_ZZX(v);
    vec.normalize();
    vec.SetMaxLength(2 * d);
  }
  
  R_Ring_Number(ZZ q_, int d_, ZZX v_vec) {
    q = q_;
    d = d_;
    vec = v_vec;
    Reduce();
    vec.normalize();
  }

  R_Ring_Number(const R_Ring_Number &v) {
    q = v.q;
    d = v.d;
    vec = v.vec;
    vec.normalize();
  }
  
  ZZ Fast_Reduce(ZZ n) const {
    return Reduce(n);
    /*
    ZZ q_half_b = q / 2;
    ZZ q_half_a = -(q - 1) / 2;
    if (n > q_half_b) {
      n -= q;
    } else if (n < q_half_a) {
      n += q;
    }
    // assert(n <= q_half_b && n >= q_half_a);
    return n; */
  }

  ZZX& Fast_Reduce(ZZX &v) const {
    for (int i = 0; i < d && i < v.rep.length(); i++) {
      v.rep[i] = Fast_Reduce(v.rep[i]);
    }
    v.normalize();
    return v;
  }

  void Fast_Reduce() {
    Fast_Reduce(vec);
  }

  static ZZ Reduce(ZZ n, ZZ modul) {
    // for odd module reduce to range [-(q - 1) / 2, (q - 1) / 2]
    // for even module reduce to range [(q - 2) / 2, q / 2]
    ZZ q_half_b, q_half_a;
    if (modul % 2 == 1) {
      q_half_b = (modul - 1) / 2;
      q_half_a = -(modul - 1) / 2;
    } else {
      q_half_b = modul / 2;
      q_half_a = -(modul - 2) / 2;
    }
    if (n > q_half_b) {
      n -= ((n - q_half_b - 1) / modul + 1) * modul;
    } else if (n < q_half_a) {
      n += (-(n - q_half_a + 1) / modul + 1) * modul;
    }
    // assert(n <= q_half_b && n >= q_half_a);
    return n;
  }    

  ZZ Reduce(ZZ n) const {
    return Reduce(n, q);
  }

  ZZX& Reduce(ZZX &v) const {
    for (int i = 0; i < d && i < v.rep.length(); i++) {
      v.rep[i] = Reduce(v.rep[i], q);
    }
    v.normalize();
    return v;
  }

  void Reduce() {
    Reduce(vec);
  }

  const R_Ring_Number operator -() const {
    return R_Ring_Number(q, d, -vec);
  }
  
  const R_Ring_Number operator +(const R_Ring_Number &v) const {
    assert(q == v.q && d == v.d);
    ZZX res = vec + v.vec;
    return R_Ring_Number(q, d, Fast_Reduce(res));
  }

  const R_Ring_Number operator -(const R_Ring_Number &v) const {
    assert(q == v.q && d == v.d);
    ZZX res = vec - v.vec;
    return R_Ring_Number(q, d, Fast_Reduce(res));
  }
  
  R_Ring_Number& operator +=(const R_Ring_Number &v) {
    if (d != v.d || q != v.q) {
      std::cout << d << " " << v.d << " " << q << " " << v.q << std::endl;
    }
    assert(d == v.d && q == v.q);
    vec += v.vec;
    Reduce();

    return *this;
  }
  
  R_Ring_Number& operator =(const R_Ring_Number &v) {
    if (q == 0) {
      q = v.q;
    }
    d = v.d;
    vec = v.vec;

    return *this;
  }
  
  R_Ring_Number& operator =(ZZ n) {
    assert(d != 0 && q > 0);
    clear(vec);
    vec.rep.SetLength(d);
    vec.rep[0] = Reduce(n, q);
    vec.normalize();

    return *this;
  }
  
  R_Ring_Number& operator =(int n) {
    return *this = to_ZZ(n);
  }

  R_Ring_Number operator *(ZZ constant) const {
    ZZX res = vec * constant;
    return R_Ring_Number(q, d, Reduce(res));
  }
  
  R_Ring_Number operator *(const R_Ring_Number &v) const {
    ZZX result = vec * v.vec;
    
    // reduction modulo x^d + 1
    for (int i = d - 1; i >= 0; i--) {
      if (i < result.rep.length() && i + d < result.rep.length()) {
	result.rep[i] -= result.rep[i + d];
	result.rep[i + d] = 0;
      }
    }

    return R_Ring_Number(q, d, Reduce(result));
  }
  

  bool operator ==(const R_Ring_Number &r) const {
    return (d == r.d && vec == r.vec);
  }

  bool operator !=(const R_Ring_Number &r) const {
    return !((*this) == r);
  }
  
  ZZ& operator [] (int index) {
    if (vec.rep.length() <= index) {
      vec.rep.SetLength(index + 1);
    }
    return vec.rep[index];
  }

  static ZZ Clamp(ZZ n, ZZ modul) {
    return Reduce(n, modul);
  }
  
  R_Ring_Number Clamp(ZZ modul) {
    for (int i = 0; i < d && i < vec.rep.length(); i++) {
      vec.rep[i] = Clamp(vec.rep[i], modul);
    }
    vec.normalize();
    return *this;
  }

  R_Ring_Number Get_Clamped(ZZ modul) {
    R_Ring_Number result(*this);
    result.Clamp(modul);
    return result;
  }

  void Change_Modul(ZZ new_q) {
    // assert(new_q > 0);
    q = new_q;
    Reduce();
  }

  void Increase_Modul(ZZ new_q) {
    // assert(q <= new_q);
    q = new_q;
  }

  void Decrease_Modul(ZZ new_q) {
    // assert(q >= new_q);
    q = new_q;
    Reduce();
  }

  R_Ring_Number Scale(ZZ q_, ZZ p, ZZ r) {
    if (r == 2) {
    // assert(q_ == q);
    // assert(p < q);
    if (p % 2 != 1 || q % 2 != 1) {
      std::cout << "p % 2 = " << p % 2 << std::endl;
      std::cout << "q % 2 = " << q % 2 << std::endl;
      
      // assert(p % 2 == 1 && q % 2 == 1);
    }
    // assert(r == 2); // for simplicity, in future should be // assert(r < p)
    //    long double fraq = (p - 1) / (long double)(q - 1); // (q - 1) / 2 should become (p - 1) / 2
    R_Ring_Number res_v(p, d);
    for (int i = 0; i < d && i < vec.rep.length(); i++) {
      ZZ desired_module = Clamp(vec.rep[i], r);
      ZZ tmp = (vec.rep[i] * (p - 1)) / (q - 1);
      ZZ value[] = {tmp, tmp + ((tmp > 0) ? -1 : 1), tmp + ((tmp > 0) ? 1 : -1), tmp + 2, tmp - 2}; // TODO: think about better approach
      ZZ max_dist;
      max_dist = -1;
      for (int j = 0; j < 5; j++) {
	value[j] = Reduce(value[j], p);
	ZZ dist = abs(tmp - value[j]);
	if (Clamp(value[j], r) == desired_module && (dist < max_dist || max_dist == -1)) {
	  max_dist = dist;
	  res_v[i] = value[j];
	}
      }
      // assert(2 * abs(res_v[i] - tmp) <= r);
      // assert(max_dist != -1);
    }
    res_v.vec.normalize();
    return res_v;
    } else {
    // assert(q_ == q);
    // assert(p < q);
    if (p % r != 1 || q % r != 1) {
      std::cout << "p = " << p << ", q = " << q << ", r = " << r << std::endl;
      std::cout << "p % r = " << p % r << std::endl;
      std::cout << "q % r = " << q % r << std::endl;
      
      assert(p % r == 1 && q % r == 1);
    }
    //    double fraq = (p - 1) / (double)(q - 1); // (q - 1) / 2 should become (p - 1) / 2
    R_Ring_Number res_v(p, d);
    for (int i = 0; i < d && i < vec.rep.length(); i++) {
      ZZ desired_module = Clamp(vec.rep[i], r);
      ZZ tmp_d = (vec.rep[i] * (p - 1)) / (q - 1);
      ZZ tmp = tmp_d;
      tmp -= Clamp(tmp, r);
      tmp += desired_module;
      ZZ value[] = {tmp, tmp + (tmp > 0 ? -r : r), tmp + (tmp > 0 ? r : -r), tmp - 2 * r, tmp + 2 * r}; // TODO: think about better approach
      ZZ max_dist;
      max_dist = -1;
      for (int j = 0; j < 5; j++) {
	value[j] = Reduce(value[j], p);
	ZZ dist = abs(tmp_d - value[j]);
	if ((dist < max_dist || max_dist == -1) && Clamp(value[j], r) == desired_module) {
	  max_dist = dist;
	  res_v[i] = value[j];
	}
      }
      if (abs(res_v[i] - tmp_d) > r + 1) {
	std::cout << "q = " << q << ", p = " << p << ", tmp_d = " << tmp_d << ", vec[i] = " << vec.rep[i] << ", tmp = " << tmp << ", res_v[i] = " << res_v[i] << ", r = " << r << std::endl;
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
      // assert(max_dist != -1);
    }
    res_v.vec.normalize();
    return res_v;
    }      
  }

  void print(void) {
    std::cout << vec.rep;
  }
    
  static R_Ring_Number Uniform_Rand(ZZ q, int d, ZZ bound = ZZ(INIT_VAL, -1)) {
    if (bound != -1) {
      bound = bound > q ? q : bound;
    } else {
      bound = q;
    }
    R_Ring_Number r(q, d);
    ZZ sub_entr = (bound % 2 == 0) ? (bound - 2) / 2 : (bound - 1) / 2;
    for (int i = 0; i < d; i++) {
      r.vec.rep[i] = RandomBnd(bound) - sub_entr;
    }
    r.vec.normalize();
    return r;
  }
  
  int Get_d(void) const {
    return d;
  }
  
  ZZ Get_q(void) const {
    return q;
  }


  double Get_Field_Expansion(void) const {
    return sqrt((double)Get_d());
  }
  
  ZZ Get_Norm(void) const {
    ZZ result;
    result = 0;
    for (int i = 0; i < d; i++) {
      result += abs(vec.rep[i]);
    }
    return result;
  }

  operator ZZ() {
    return vec.rep[0];
  }
};

/*ostream& operator<<(ostream& s, const R_Ring_Number& a) {
  s << a.vec;
  return s;
  }*/

#endif /* _R_RING_NUMBER_H_ */
