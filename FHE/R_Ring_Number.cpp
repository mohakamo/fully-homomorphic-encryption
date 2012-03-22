#include "R_Ring_Number.h"
#include "R_Ring_Vector.h"
#include "R_Ring_Matrix.h"

ostream& operator<<(ostream& s, const R_Ring_Number& a) {
  s << a.vec;
  return s;
}

ostream& operator<<(ostream& s, const R_Ring_Vector& a) {
  s << "{";
  for (int i = 0; i < a.Get_Dimension(); i++) {
    s << a[i];
    if (i + 1 < a.Get_Dimension()) {
      s << ", ";
    }
  }
  s << "}";

  return s;
}

ostream& operator<<(ostream& s, const R_Ring_Matrix& a) {
  for (long long i = 0; i < a.Get_Noof_Rows(); i++) {
    for (long long j = 0; j < a.Get_Noof_Columns(); j++) {
      s << a(i, j);
      if (j + 1 < a.Get_Noof_Columns()) {
	s << ", ";
      } else {
	s << "\n";
      }
    }
  }

  return s;
}

/*R_Ring_Number R_Ring_Number::operator +(const R_Ring_Number &v) const {
  // assert(q == v.q && d == v.d);
  R_Ring_Number res_v;
  res_v.Initialize(v.q, v.d);
  for (int i = 0; i < v.d; i++) {
    res_v.vec[i] = vec[i] + v.vec[i];
  }
	return res_v;
	}*/
