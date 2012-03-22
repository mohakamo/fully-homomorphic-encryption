#include "R_Ring_Vector.h"
#include "defs.h"

void Bit_Decomposition(const R_Ring_Vector &x, ZZ q, R_Ring_Vector &result) {
  assert(q == x.Get_q());
  int noof_vectors = NumBits(q);
  if (result.Get_Dimension() != noof_vectors * x.Get_Dimension() ||
      result.Get_d() != x.Get_d() ||
      result.Get_q() != x.Get_q()) {
    result = R_Ring_Vector(q, x.Get_d(), noof_vectors * x.Get_Dimension());
  }
  int index;
  for (int i = 0; i < x.Get_Dimension(); i++) {
    index = i * noof_vectors;
    for (int p = 0; p < noof_vectors; p++) {
      for (int j = 0; j < x.Get_d(); j++) {
	// ZZ n = R_Ring_Number::Clamp(x[i][j], q);
	ZZ n = x[i][j]; // do not need clamp, since q == x.Get_q() due to first assertion
	if (n < 0) {
	  n += q;
	}
	//	  // assert(n >= 0 && n < q);
	result[p + index][j] = (n >> p) & 1; 
      }
    }
  }
}

void Powersof2(const R_Ring_Vector &x, ZZ q, R_Ring_Vector &result) {
  static ZZ TWO = ZZ(INIT_VAL, 2);
  int noof_vectors = NumBits(q);
  if (result.Get_Dimension() != noof_vectors * x.Get_Dimension() ||
      result.Get_d() != x.Get_d() ||
      result.Get_q() != q) {
    result = R_Ring_Vector(q, x.Get_d(), noof_vectors * x.Get_Dimension());
    assert(result.Get_q() == q);
  }
  int index;
  
  for (int i = 0; i < x.Get_Dimension(); i++) {
    result[i * noof_vectors] = x[i];
    for (int p = 1; p < noof_vectors; p++) {
      index = p + i * noof_vectors;
      result[index] = result[index - 1] * TWO; // not to have overflowing
    }
  }
}

