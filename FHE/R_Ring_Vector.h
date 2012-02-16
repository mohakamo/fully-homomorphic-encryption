/*
 *  FHE_R_Ring_Vector.h
 *  Created by valerini on 1/5/12.
 */
#ifndef _R_RING_VECTOR_H_
#define _R_RING_VECTOR_H_

#include "R_Ring_Number.h"
#include <cassert>

class R_Ring_Matrix;

class R_Ring_Vector {
  R_Ring_Number* vec;
  int dimension;
public:
  
  ~R_Ring_Vector() {
    if (vec != NULL) {
      delete [] vec;
    }
  }
  
  R_Ring_Vector() {
    dimension = 0;
    vec = NULL;
  }
  
  R_Ring_Vector(ZZ q, int d, int dimension)  {
    vec = NULL;
    Initialize(q, d, dimension);
  }

  R_Ring_Vector(int dimension_, R_Ring_Number* vec_) {
    dimension = dimension_;
    vec = new R_Ring_Number [dimension];
    for (int i = 0; i < dimension; i++) {
      vec[i] = vec_[i];
    }
  }
      
  
  R_Ring_Vector(const R_Ring_Vector &r) {
    dimension = 0;
    vec = NULL;
    *this = r;
  }  
  
  R_Ring_Vector& operator =(const R_Ring_Vector &r) {
    if (dimension != r.dimension && r.dimension != 0) {
      if (vec != NULL) {
	delete [] vec;
      }
      vec = new R_Ring_Number[r.dimension];
      dimension = r.dimension;
    }

    for (int i = 0; i < dimension; i++) {
      vec[i] = r.vec[i];
    }
    if (dimension != 0 && Get_q() != r.Get_q()) {
      for (int i = 0; i < dimension; i++) {
	vec[i].Change_Modul(r.Get_q());
      }
    }
    return *this;
  }

  void Initialize(ZZ q, int d, int dimension_) { 
   if (vec != NULL) {
      delete [] vec;
    }
    dimension = dimension_;
    vec = new R_Ring_Number [dimension];
    for (int i = 0; i < dimension; i++) {
      vec[i].Initialize(q, d);
    }
  }
  
  R_Ring_Number& operator[] (unsigned int index) const {
    /* // ommit corectness check for speed purposes
       if (index >= 0 && index < vec.size()) {
       return Ring_Number.Invalid_Number;
       }
    */
    return vec[index];
  }
  
  int Get_Dimension() const {
    return dimension;
  }
  
  R_Ring_Vector operator + (const R_Ring_Vector &r) const {
    // assert(dimension == r.dimension && r.dimension != 0);
    R_Ring_Vector res(Get_q(), Get_d(), dimension);
    for (int i = 0; i < dimension; i++) {
      res[i] = vec[i] + r.vec[i];
    }
    
    return res;
  }

  R_Ring_Vector operator - (const R_Ring_Vector &r) const {
    // assert(dimension == r.dimension && r.dimension != 0);
    R_Ring_Vector res(Get_q(), Get_d(), dimension);
    for (int i = 0; i < dimension; i++) {
      res[i] = vec[i] - r.vec[i];
    }
    
    return res;
  }
  
  /*R_Ring_Matrix operator R_Ring_Matrix() const {
    return R_Ring_Matrix(Get_q(), Get_d(), vec, dimension, 1);
    }*/
  
  R_Ring_Number Dot_Product(const R_Ring_Vector &r) const {
    // assert(dimension == r.dimension && r.dimension != 0);
    R_Ring_Number res(Get_q(), Get_d()); // the number will be zero by default
    for (int i = 0; i < dimension; i++) {
      res += vec[i] * r.vec[i];
    }
    return res;
  }

  R_Ring_Vector operator *(const R_Ring_Vector &v) const {
    // assert(Get_Dimension() == v.Get_Dimension());
    R_Ring_Vector res(Get_q(), Get_d(), Get_Dimension());
    for (int i = 0; i < Get_Dimension(); i++) {
      res.vec[i] = vec[i] * v.vec[i];
    }
    return res;
  }

  R_Ring_Vector operator *(const R_Ring_Matrix &m) const;

  
  R_Ring_Vector operator *(ZZ number) const {
    R_Ring_Vector res(Get_q(), Get_d(), Get_Dimension());
    for (int i = 0; i < Get_Dimension(); i++) {
      res[i] = vec[i] * number;
    }
    return res;
  }

  R_Ring_Vector operator *(int number) const {
    return *this * ZZ(INIT_VAL, number);
  }

  bool operator ==(const R_Ring_Vector &c) const {
    if (c.Get_Dimension() != Get_Dimension()) {
      return false;
    }
    for (int i = 0; i < Get_Dimension(); i++) {
      if (vec[i] != c.vec[i]) {
	return false;
      }
    }
    return true;
  }

  bool operator !=(const R_Ring_Vector &c) const {
    return !(*this == c);
  }
  
  ZZ Get_q(void) const {
    // assert(dimension != 0);
    return vec[0].Get_q();
  }
  
  int Get_d(void) const {
    // assert(dimension != 0);
    return vec[0].Get_d();
  }
  
  double Get_Field_Expansion(void) const {
    // assert(dimension != 0);
    return vec[0].Get_Field_Expansion();
  }

  // Returns the vector of elements start_index, ..., end_index
  R_Ring_Vector Get_Sub_Vector(int start_index, int end_index) const {
    // assert(end_index >= start_index && start_index >= 0 && end_index < dimension);
    R_Ring_Vector res_v(Get_q(), Get_d(), end_index - start_index + 1);
    for (int i = start_index; i <= end_index; i++) {
      res_v[i - start_index] = vec[i];
    }
    return res_v;
  }
  
  R_Ring_Matrix Transpose() const;
  
  R_Ring_Vector Tensor_Product(const R_Ring_Vector &r) const {
    // assert(Get_Dimension() == r.Get_Dimension());
    int new_dimension = (Get_Dimension() * (Get_Dimension() + 1)) / 2;
    R_Ring_Vector res(Get_q(), Get_d(), new_dimension);
    int index = 0;
    for (int i = 0; i < Get_Dimension(); i++) {
      for (int j = i; j < Get_Dimension(); j++) {
	res[index++] = vec[i] * r.vec[j];
      }
    }
    // assert(index == new_dimension);
    return res;
  }

  void Increase_Modul(ZZ new_q) {
    for (int i = 0; i < dimension; i++) {
      vec[i].Increase_Modul(new_q);
    }
  }

  void Decrease_Modul(ZZ new_q) {
    for (int i = 0; i < dimension; i++) {
      vec[i].Decrease_Modul(new_q);
    }
  }

  static R_Ring_Vector Uniform_Rand(ZZ  q, int d, int dimension, ZZ bound = ZZ(INIT_VAL, -1)) {
    R_Ring_Vector res(q, d, dimension);
    for (int i = 0; i < dimension; i++) {
      res[i] = R_Ring_Number::Uniform_Rand(q, d, bound);
    }
    return res;
  }

  void Clamp(ZZ modul) {
    for (int i = 0; i < dimension; i++) {
      vec[i].Clamp(modul);
    }
  }

  R_Ring_Vector Get_Clamped(ZZ modul) {
    R_Ring_Vector res = *this;
    res.Clamp(modul);
    return res;
  }

  void print(void) const {
    std::cout << "{";
    for (int i = 0; i < dimension; i++) {
      vec[i].print();
      if (i + 1 < dimension) {
	std::cout << ", ";
      }
    }
    std::cout << "}";
  }
};

#endif /* _R_RING_VECTOR_H_ */
