#include "R_Ring_Vector.h"
#include "R_Ring_Matrix.h"

int Ring_Number_d = 0;
ZZ Ring_Number_q[10];

R_Ring_Matrix R_Ring_Vector::Transpose() const {
  R_Ring_Matrix res_m(Get_q(), Get_d(), 1, Get_Dimension());

  for (int i = 0; i < Get_Dimension(); i++) {
    res_m(0, i) = vec[i];
  }
  return res_m;
}

R_Ring_Vector R_Ring_Vector::operator *(const R_Ring_Matrix &m) const {
  assert(m.Get_Noof_Rows() == dimension);
  R_Ring_Vector res_v(Get_q(), Get_d(), m.Get_Noof_Columns());
    
  for (int i = 0; i < m.Get_Noof_Columns(); i++) {
    for (int j = 0; j < dimension; j++) {
      res_v.vec[i] += vec[j] * m(j, i);
    }
  }
  return res_v;
}
