#include "FHE.h"

// The result is modul c1.Get_q(), to switch to other modul, use Scale
R_Ring_Vector FHE_Cipher_Text::Switch_Key(R_Ring_Matrix A, R_Ring_Vector c1) {
  assert(A.Get_q() == c1.Get_q());
  R_Ring_Vector v = FHE::Bit_Decomposition(c1, A.Get_q());
  R_Ring_Matrix m = v.Transpose();
  return (m * A).Get_Vector();
}

R_Ring_Vector FHE_Cipher_Text::Scale(R_Ring_Vector &x, int q, int p, int r) {
  assert(p < q);
  R_Ring_Vector rn(p, x.Get_d(), x.Get_Dimension());
  for (int i = 0; i < x.Get_Dimension(); i++) {
    rn[i] = x[i].Scale(q, p, r);
  }
  return rn;
}
  
// A_j_j_1 - matrix that transforms key from level j to level j - 1
// q_j_1 < q_j
// Moves c from level c.second to level c.second - 1
void FHE_Cipher_Text::Refresh(Pair<R_Ring_Vector, int> &c, FHE_Public_Key_Type *pk) {
  assert(c.second - 1 >= 0);
  int L = ((*my_pk).size() - 1) / 2;
  R_Ring_Vector c1 = FHE::Powersof2(c.first, (*pk)[c.second].Get_q());
  R_Ring_Vector c2 = Scale(c1, (*pk)[c.second].Get_q(), (*pk)[c.second - 1].Get_q(), 2);
  assert(c2.Get_q() == (*pk)[L + c.second].Get_q());
  R_Ring_Vector c3 = Switch_Key((*pk)[L + c.second], c2);

  c = Pair<R_Ring_Vector, int>(c3, c.second - 1);
}
  
/*
void FHE_Cipher_Text::Refresh(Pair<R_Ring_Vector, int> &c, R_Ring_Matrix A_j_j_1, int q_j, int q_j_1) {
  assert(c.first.Get_q() == q_j);
  R_Ring_Vector c1 = FHE::Powersof2(c.first, q_j);
  R_Ring_Vector c2 = Switch_Key(A_j_j_1, c2);

  assert(c2.Get_q() == c1.Get_q());
  R_Ring_Vector c3 = Scale(c1, q_j, q_j_1, 2);
  assert(c3.Get_q() == q_j_1);
  c = Pair<R_Ring_Vector, int>(c3, c.second - 1);
  } */

void FHE_Cipher_Text::Update_To_Same_Level(Pair<R_Ring_Vector, int> &c1, Pair<R_Ring_Vector, int> &c2) {
  int L = ((*my_pk).size() - 1) / 2;
  if (c1.second != c2.second) {
    if (c1.second < c2.second) {
      while (c1.second != c2.second) {
	assert(c2.second - 1 >= 0);
	// switch c2 to next level
	// c2 goes from level c2.second to level c2.second - 1
	Refresh(c2, my_pk);
      }
    } else {
      while (c1.second != c2.second) {
	// switch c1 to higher level
	Refresh(c1, my_pk);
      }
    }
  }
}

Pair<R_Ring_Vector, int> FHE_Cipher_Text::Add(Pair<R_Ring_Vector, int> &c1, Pair<R_Ring_Vector, int> &c2) {
  int L = ((*my_pk).size() - 1) / 2;
  Update_To_Same_Level(c1, c2);
  R_Ring_Vector c3_p = c1.first + c2.first;
  assert(c1.first.Get_Dimension() == c2.first.Get_Dimension());
  R_Ring_Vector c3(c3_p.Get_q(), c3_p.Get_d(), (c1.first.Get_Dimension() * (c1.first.Get_Dimension() + 1)) / 2);
  for (int i = 0; i < c1.first.Get_Dimension(); i++) {
    c3[i] = c3_p[i];
  }
  assert(c1.second >= 1);
  // going from c1.second level to c1.second - 1 level
  Pair<R_Ring_Vector, int> result(c3, c1.second);
  Refresh(result, my_pk);
  return result;
}

Pair<R_Ring_Vector, int> FHE_Cipher_Text::Mult(Pair<R_Ring_Vector, int> &c1, Pair<R_Ring_Vector, int> &c2) {
  int L = ((*my_pk).size() - 1) / 2;
  Update_To_Same_Level(c1, c2);
  assert(c1.second == c2.second);
  assert(c1.first.Get_Dimension() == c2.first.Get_Dimension());
  assert(c1.first.Get_q() == c2.first.Get_q());
  int dimension = c1.first.Get_Dimension();
  int new_dimension = (dimension * (dimension + 1)) / 2;
  R_Ring_Vector c3(c1.first.Get_q(), c1.first.Get_d(), new_dimension);
  int index = 0;
  for (int i = 0; i < dimension; i++) {
    c3[index++] = c1.first[i] * c2.first[i];
    for (int j = i + 1; j < dimension; j++) {
      c3[index++] = c1.first[i] * c2.first[j] + c1.first[j] * c2.first[i];
    }
  }
  assert(index == new_dimension);
  assert(c1.second >= 1);
  // going from c1.second level to c1.second - 1 level
  Pair<R_Ring_Vector, int> result(c3, c1.second);
  Refresh(result, my_pk);
  return result;
}
