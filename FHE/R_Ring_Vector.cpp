#include "R_Ring_Vector.h"
#include "R_Ring_Matrix.h"

R_Ring_Matrix R_Ring_Vector::Transpose() const {
  R_Ring_Matrix res_m(Get_q(), Get_d(), 1, Get_Dimension());

  for (int i = 0; i < Get_Dimension(); i++) {
    res_m(0, i) = vec[i];
  }
  return res_m;
}
