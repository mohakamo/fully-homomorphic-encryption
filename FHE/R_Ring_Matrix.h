/*  matrix.h
 *  Created by valerini on 1/5/12
 */
#ifndef _FHE_R_RING_MATRIX_H_
#define _FHE_R_RING_MATRIX_H_

#include <stdio.h>
#include <vector>
#include "R_Ring_Vector.h"

class R_Ring_Number;
class R_Ring_Vector;

class R_Ring_Matrix {
  R_Ring_Number *matrix;
  long long noof_columns, noof_rows;
  
public:
  
  ~R_Ring_Matrix() {
    if (matrix != NULL) {
	delete [] matrix;
      }
  }
  
  // Added for array allocation
  R_Ring_Matrix() {
    noof_columns = noof_rows = 0;
    matrix = NULL;
  }
  
  R_Ring_Matrix(ZZ q, int d, const R_Ring_Number *matrix, long long noof_rows, long long noof_columns) {
    matrix = NULL;
    Initialize(q, d, noof_rows, noof_columns);
    memcpy(this->matrix, matrix, sizeof(R_Ring_Number) * noof_columns * noof_rows);
  }
  
  R_Ring_Matrix(ZZ q, int d, long long noof_rows, long long noof_columns) {
    matrix = NULL;
    Initialize(q, d, noof_rows, noof_columns);
  }
  
  R_Ring_Matrix(const R_Ring_Vector &v) {
    matrix = NULL;
    Initialize(v.Get_q(), v.Get_d(), v.Get_Dimension(), 1);
    for (long long i = 0; i < v.Get_Dimension(); i++) {
      matrix[i] = v[i];
    }
  }
  
  R_Ring_Matrix(const R_Ring_Matrix &m) {
    noof_columns = noof_rows = 0;
    matrix = NULL;
    *this = m;
  }

 R_Ring_Matrix& operator =(const R_Ring_Matrix &m) {
    long long dimension = m.noof_rows * m.noof_columns;
    
    if (noof_rows != m.noof_rows || noof_columns != m.noof_columns) {
      if (matrix != NULL) {
	delete [] matrix;
      }
      
      matrix = new R_Ring_Number [dimension];
      noof_columns = m.noof_columns;
      noof_rows = m.noof_rows;
    }
    for (long long i = 0; i < dimension; i++) {
      matrix[i] = m.matrix[i];
    }
    if (dimension != 0 && Get_q() != m.Get_q()) {
      for (long long i = 0; i < dimension; i++) {
	matrix[i].Change_Modul(m.Get_q());
      }
    }
    return *this;
  }
  
  // Added for initialization during array allocation
  void Initialize(ZZ q, int d, long long _noof_rows, long long _noof_columns) {
    if (matrix != NULL) {
      delete [] matrix;
    }
    noof_rows = _noof_rows;
    noof_columns = _noof_columns;
    long long dimension = noof_rows * noof_columns;
    matrix = new R_Ring_Number [dimension];
    for (long long i = 0; i < dimension; i++) {
      matrix[i].Initialize(q, d);
    }
  }
  
  // this operator allows access to the matrix as to two dimensional array
  /*  // !!!! this operator is not sufficient  since it returns not an array of references
  // thus I can't write matrix[i][j] = something;
  // for write access use () operator
  R_Ring_Number * operator[] (int row_index) const {
    return &matrix[row_index * noof_columns];
    }*/
  
  R_Ring_Number& operator() (long long row_number, long long column_number) const {
    return matrix[row_number * noof_columns + column_number];
  }

  R_Ring_Vector Get_Vector() const {
    assert(noof_rows == 1);
    return R_Ring_Vector(noof_columns, matrix);
  }
  
  R_Ring_Matrix Get_Transpose() const {
    R_Ring_Matrix res_m (Get_q(), Get_d(), noof_columns, noof_rows);
    
    for (long long r = 0; r < noof_rows; r++) {
      for (long long c = 0; c < noof_columns; c++) {
	res_m(c, r) = (*this)(r, c);
      }
    }
    
    return res_m;
  }
     
  R_Ring_Matrix operator +(const R_Ring_Matrix &m) const {
    assert(noof_rows == m.noof_rows && noof_columns == m.noof_columns);
    
    R_Ring_Matrix res_m(Get_q(), Get_d(), noof_rows, noof_columns);
    
    for (long long i = 0; i < noof_columns * noof_rows; i++) {
      res_m.matrix[i] = matrix[i] + m.matrix[i];
    }
    return res_m;
  }
  
  R_Ring_Matrix operator -() const {
    R_Ring_Matrix res_m(Get_q(), Get_d(), noof_rows, noof_columns);
    
    for (long long i = 0; i < noof_rows * noof_columns; i++) {
      res_m.matrix[i] = -matrix[i];
    }
    return res_m;
  }
  
  R_Ring_Vector operator *(const R_Ring_Vector &v) const {
    assert(noof_columns == v.Get_Dimension());
    R_Ring_Vector res_v(Get_q(), Get_d(), noof_rows);
    
    for (long long i = 0; i < noof_rows; i++) {
      for (long long j = 0; j < noof_columns; j++) {
	res_v[i] += (*this)(i, j) * v[j];
      }
    }
    return res_v;
  }
  
  R_Ring_Matrix operator *(const R_Ring_Matrix &m) const {
    assert(noof_columns == m.noof_rows);
    R_Ring_Matrix res_m(Get_q(), Get_d(), noof_rows, m.noof_columns);
    R_Ring_Number res(Get_q(), Get_d());
	
    for (long long i = 0; i < noof_rows; i++) {
      for (long long j = 0; j < m.noof_columns; j++) {
	res = 0;
	for (long long k = 0; k < noof_columns; k++) {
	  res += (*this)(i, k) * m(k, j);
	}
	res_m(i, j) = res;
      }
    }
    return res_m;
  }
  
  void Set_Block(long long r, long long c, const R_Ring_Matrix &m) {
    assert(noof_rows == r + m.noof_rows && noof_columns == c + m. noof_columns);
    for (long long i = r; i < noof_rows; i++) {
      for (long long j = c; j < noof_columns; j++) {
	matrix[i * noof_columns + j] = m.matrix[(i - r) * m.noof_columns + j - c];
      }
    }
  }
  
  void Set_Column(long long c, const R_Ring_Vector &r) {
    assert(c >= 0 && c < noof_columns && noof_rows == r.Get_Dimension());
    for (long long i = 0; i < noof_rows; i++) {
      matrix[i * noof_columns + c] = r[i];
    }
  }
  
  R_Ring_Matrix & Add_To_Column(long long column_number, const R_Ring_Vector &r) {
    assert(column_number >= 0 && column_number < noof_columns && noof_rows == r.Get_Dimension());
    for (long long i = 0; i < noof_rows; i++) {
      (*this)(i, column_number) += r[i];
    }
    return *this;
  }

  void Increase_Modul(ZZ new_q) {
    // assert(new_q >= Get_q());
    for (long long i = 0; i < noof_rows * noof_columns; i++) {
      matrix[i].Increase_Modul(new_q);
    }
  }
  
  static R_Ring_Matrix Uniform_Rand(ZZ q, int d, long long noof_rows,
				    long long noof_columns,
				    ZZ bound = ZZ(INIT_VAL, -1)) {
    R_Ring_Matrix res(q, d, noof_rows, noof_columns);
    for (long long i = 0; i < noof_rows * noof_columns; i++) {
      res.matrix[i] = R_Ring_Number::Uniform_Rand(q, d, bound);
    }
    return res;
  }

  ZZ Get_q(void) const {
    assert(noof_rows != 0 && noof_columns != 0);
    return matrix[0].Get_q();
  }
  
  int Get_d(void) const {
    assert(noof_rows != 0 && noof_columns != 0);
    return matrix[0].Get_d();
  }
  
  long long Get_Noof_Rows(void) const {
    return noof_rows;
  }
  
  long long Get_Noof_Columns(void) const {
    return noof_columns;
  }

  void print(void) const {
    for (long long i = 0; i < noof_rows; i++) {
      for (long long j = 0; j < noof_columns; j++) {
	matrix[i * noof_columns + j].print();
	if (j + 1 < noof_columns) {
	  std::cout << ", ";
	} else {
	  std::cout << std::endl;
	}
      }
    }
  }
};

ostream& operator<<(ostream& s, const R_Ring_Matrix& a);

#endif /* _FHE_R_RING_MATRIX_H_ */
