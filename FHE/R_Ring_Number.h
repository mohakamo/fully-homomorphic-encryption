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
#include "defs.h"

NTL_CLIENT

class R_Ring_Number {
  ZZ q;

public:
  ZZX vec;

  R_Ring_Number() : q(), vec() {}
  
  R_Ring_Number(ZZ q_, int d_) {
    assert(Ring_Number_d != 0);
    q = q_;
    vec.SetMaxLength(2 * Ring_Number_d);
    vec.rep.SetLength(Ring_Number_d);
  }
  
  void Initialize(ZZ q_, int d_) {
    q = q_;
    vec.SetMaxLength(2 * Ring_Number_d);
    vec.rep.SetLength(Ring_Number_d);
  }

  R_Ring_Number(ZZ q_, int d_, ZZ array[]) {
    assert(Ring_Number_d != 0);
    q = q_;
    vec_ZZ v;
    v.SetLength(Ring_Number_d);
    for (int i = 0; i < Ring_Number_d; i++) {
      v[i] = Reduce(array[i]);
    }
    vec = to_ZZX(v);
    vec.normalize();
    vec.SetMaxLength(2 * Ring_Number_d);
  }
  
  R_Ring_Number(ZZ q_, int d_, ZZX v_vec) {
    assert(Ring_Number_d != 0);
    q = q_;
    vec = v_vec;
    // NB!
    //    Reduce();
    vec.normalize();
  }

  R_Ring_Number(const R_Ring_Number &v) {
    assert(Ring_Number_d != 0);
    q = v.q;
    vec = v.vec;
    vec.normalize();
  }

  static ZZ Fast_Reduce(ZZ n, ZZ modul) {
    return Reduce(n, modul);
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

  static ZZX& Fast_Reduce(ZZX &v, ZZ modul) {
    for (int i = 0; i < Ring_Number_d && i < v.rep.length(); i++) {
      v.rep[i] = Fast_Reduce(v.rep[i], modul);
    }
    v.normalize();
    return v;
  }

  ZZX& Fast_Reduce(ZZX &v) const {
    return Fast_Reduce(v, q);
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

  static ZZX& Reduce(ZZX &v, ZZ modul) {
    for (int i = 0; i < Ring_Number_d && i < v.rep.length(); i++) {
      v.rep[i] = Reduce(v.rep[i], modul);
    }
    v.normalize();
    return v;
  }

  ZZX& Reduce(ZZX &v) const {
    return Reduce(v, q);
  }

  void Reduce() {
    Reduce(vec);
  }

  const R_Ring_Number operator -() const {
    return R_Ring_Number(q, Ring_Number_d, -vec);
  }
  
  const R_Ring_Number operator +(const R_Ring_Number &v) const {
    if (q != v.q) {
      std::cout << "q = " << q << std::endl;
      std::cout << "v.q = " << v.q << std::endl;
      assert(q == v.q);
    }
    ZZX res = vec + v.vec;
    return R_Ring_Number(q, Ring_Number_d, Fast_Reduce(res));
  }

  static R_Ring_Number sum(const R_Ring_Number &v1, const R_Ring_Number &v2, bool doReduce = true) {
    assert(v1.q == v2.q);
    ZZX res = v1.vec + v2.vec;
    return R_Ring_Number(v1.q, Ring_Number_d, doReduce ? Fast_Reduce(res, v1.q) : res);
  }

  const R_Ring_Number operator +(const ZZ &v) const {
    R_Ring_Number res(q, Ring_Number_d);
    res = v;
    return (*this) + res;
  }

  const R_Ring_Number operator -(const R_Ring_Number &v) const {
    assert(q == v.q);
    ZZX res = vec - v.vec;
    return R_Ring_Number(q, Ring_Number_d, Fast_Reduce(res));
  }
  
  R_Ring_Number& operator +=(const R_Ring_Number &v) {
    if (q != v.q) {
      std::cout << q << " " << v.q << std::endl;
    }
    assert(q == v.q);
    vec += v.vec;
    Reduce();

    return *this;
  }
  
  R_Ring_Number& operator =(const R_Ring_Number &v) {
    if (q == 0) {
      q = v.q;
    }
    vec = v.vec;

    return *this;
  }
  
  R_Ring_Number& operator =(ZZ n) {
    assert(Ring_Number_d != 0 && q > 0);
    clear(vec);
    vec.rep.SetLength(Ring_Number_d);
    vec.rep[0] = Reduce(n, q);
    vec.normalize();

    return *this;
  }
  
  R_Ring_Number& operator =(int n) {
    return *this = to_ZZ(n);
  }

  R_Ring_Number operator /(ZZ constant) const {
    if (constant == 1) {
      return *this;
    }
    ZZX res = vec;
    for (int i = 0; i < res.rep.length(); i++) {
      if (res.rep[i] < 0) {
	res.rep[i] = (res.rep[i] + q) / constant;
      } else {
	res.rep[i] = res.rep[i] / constant;
      }
    }
    return R_Ring_Number(q, Ring_Number_d, Reduce(res));
  }

  R_Ring_Number operator *(ZZ constant) const {
    ZZX res = vec * constant;
    return R_Ring_Number(q, Ring_Number_d, Reduce(res));
  }
  
  R_Ring_Number operator *(const R_Ring_Number &v) const {
    ZZX result = vec * v.vec;
    
    // reduction modulo x^d + 1
    for (int i = Ring_Number_d - 1; i >= 0; i--) {
      if (i < result.rep.length() && i + Ring_Number_d < result.rep.length()) {
	result.rep[i] -= result.rep[i + Ring_Number_d];
	result.rep[i + Ring_Number_d] = 0;
      }
    }

    return R_Ring_Number(q, Ring_Number_d, Reduce(result));
  }
  
  static R_Ring_Number mul(const R_Ring_Number &v1, const R_Ring_Number &v2, bool doReduction = true) {
    ZZX result = v1.vec * v2.vec;
    
    // reduction modulo x^d + 1
    for (int i = Ring_Number_d - 1; i >= 0; i--) {
      if (i < result.rep.length() && i + Ring_Number_d < result.rep.length()) {
	result.rep[i] -= result.rep[i + Ring_Number_d];
	result.rep[i + Ring_Number_d] = 0;
      }
    }

    return R_Ring_Number(v1.Get_q(), Ring_Number_d, doReduction ? Reduce(result, v1.Get_q()) : result);
  }

  bool operator ==(const R_Ring_Number &r) const {
    return (vec == r.vec);
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
    for (int i = 0; i < Ring_Number_d && i < vec.rep.length(); i++) {
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
    R_Ring_Number res_v(p, Ring_Number_d);
    for (int i = 0; i < Ring_Number_d && i < vec.rep.length(); i++) {
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
    R_Ring_Number res_v(p, Ring_Number_d);
    for (int i = 0; i < Ring_Number_d && i < vec.rep.length(); i++) {
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
    R_Ring_Number r(q, Ring_Number_d);
    ZZ sub_entr = (bound % 2 == 0) ? (bound - 2) / 2 : (bound - 1) / 2;
    for (int i = 0; i < Ring_Number_d; i++) {
      r.vec.rep[i] = RandomBnd(bound) - sub_entr;
    }
    r.vec.normalize();
    return r;
  }
  
  int Get_d(void) const {
    return Ring_Number_d;
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
    for (int i = 0; i < Ring_Number_d; i++) {
      result += abs(vec.rep[i]);
    }
    return result;
  }

  operator ZZ() {
    return vec.rep[0];
  }
};

ostream& operator<<(ostream& s, const R_Ring_Number& a);

#endif /* _R_RING_NUMBER_H_ */
