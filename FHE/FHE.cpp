#include "FHE.h"

// The result is modul c1.Get_q(), to switch to other modul, use Scale
R_Ring_Vector FHE_Cipher_Text::Switch_Key(R_Ring_Matrix A, R_Ring_Vector c1) {
  assert(A.Get_q() == c1.Get_q());
  R_Ring_Vector v = FHE::Bit_Decomposition(c1, A.Get_q());
  R_Ring_Matrix m = v.Transpose();
  return (m * A).Get_Vector();
}

R_Ring_Vector FHE_Cipher_Text::Scale(R_Ring_Vector &x, ZZ q, ZZ p, int r) {
  assert(p < q);
  R_Ring_Vector rn(p, x.Get_d(), x.Get_Dimension());
  for (int i = 0; i < x.Get_Dimension(); i++) {
    rn[i] = x[i].Scale(q, p, ZZ(INIT_VAL, r));
  }
  return rn;
}
  
// A_j_j_1 - matrix that transforms key from level j to level j - 1
// q_j_1 < q_j
// Moves c from level c.second to level c.second - 1
void FHE_Cipher_Text::Refresh(Pair<R_Ring_Vector, int> &c, FHE_Public_Key_Type *pk, FHE_Secret_Key_Type *sk) {
  assert(c.second - 1 >= 0);
  int L = ((*my_pk).size() - 1) / 2;
  R_Ring_Vector c1 = FHE::Powersof2(c.first, (*pk)[c.second].Get_q());

  /**** debug code ****/
  // noise check
  ZZ B;
  R_Ring_Vector s_pp;
  if (sk != NULL) {
    s_pp = FHE::Bit_Decomposition((*sk)[c.second].Tensor_Product((*sk)[c.second]), (*sk)[c.second].Get_q());
    B = c1.Dot_Product(s_pp).Get_Norm();
    ZZ required_upper_bound = (*sk)[c.second].Get_q() / 2 - (*sk)[c.second - 1].Get_q() * c.first.Get_d() * (int)(c.first.Get_Field_Expansion() + 1) * (((*sk)[c.second].Get_Dimension() - 1) * (*sk)[c.second].Get_Dimension()) / 2 * NumBits(c.first.Get_q()) / (*sk)[c.second].Get_q();
    assert(B < required_upper_bound);
    ZZ expected_B = (*sk)[c.second - 1].Get_q() * B / (*sk)[c.second].Get_q() + c.first.Get_d() * (int)(c.first.Get_Field_Expansion() + 1) * (((*sk)[c.second].Get_Dimension() - 1) * (*sk)[c.second].Get_Dimension()) / 2 * NumBits(c.first.Get_q());
    assert(expected_B * 2 < c.first.Get_q());
  }
  /**** endof debug code ****/
  
  R_Ring_Vector c2 = Scale(c1, (*pk)[c.second].Get_q(), (*pk)[c.second - 1].Get_q(), my_p);
  
  /**** debug code ****/
  // noise check
  if (sk != NULL) {
    s_pp.Decrease_Modul(c2.Get_q());
    B = c2.Dot_Product(s_pp).Get_Norm();
    assert(B * 2 < c2.Get_q());
  }
  /**** endof debug code ****/
  
  assert(c2.Get_q() == (*pk)[L + c.second].Get_q());
  R_Ring_Vector c3 = Switch_Key((*pk)[L + c.second], c2);

  /**** debug code ****/
  // noise check
  if (sk != NULL) {
    ZZ new_B = c3.Dot_Product((*sk)[c.second - 1]).Get_Norm();
    // 2 multiplier before second summand stands for B_{\xi}
    ZZ expected_B = B + 2 * sqrt((double)c.first.Get_d()) * 2 * c.first.Get_Field_Expansion() * (((*sk)[c.second].Get_Dimension() - 1) * (*sk)[c.second].Get_Dimension()) / 2 * NumBits(c.first.Get_q()) * NumBits(c.first.Get_q());
    assert(expected_B * 2 < c3.Get_q());
    assert(new_B * 2 < c2.Get_q());
  }
  /**** endof debug code ****/

  /*  if (my_p == 3) {
    std::cout << "**************** Refresh operation *******************" << std::endl;
    std::cout << "c = ";
    c.first.print();
    std::cout << std::endl;
    std::cout << "c1 = ";
    c1.print();
    std::cout << std::endl;
    std::cout << "Scale(" << (*pk)[c.second].Get_q() << ", " << (*pk)[c.second - 1].Get_q() << ", " <<  my_p << ")" << std::endl;
    std::cout << "c2 = ";
    c2.print();
    std::cout << std::endl;
    std::cout << "c3 = ";
    c3.print();
    std::cout << std::endl;
    }*/

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
	R_Ring_Vector res_c(c2.first.Get_q(), c2.first.Get_d(), (c2.first.Get_Dimension() * (c2.first.Get_Dimension() + 1)) / 2);
	for (int i = 0; i < c2.first.Get_Dimension(); i++) {
	  res_c[i] = c2.first[i];
	}
	Pair<R_Ring_Vector, int> result(res_c, c2.second);
	assert(res_c.Get_q() == (*my_pk)[c2.second].Get_q());
	Refresh(result, my_pk);
	c2 = result;
	assert(result.first.Get_q() == (*my_pk)[c2.second].Get_q());
	assert(c2.first.Get_q() == (*my_pk)[c2.second].Get_q());

	//	Refresh(c2, my_pk);
      }
    } else {
      while (c1.second != c2.second) {
	// switch c1 to higher level
	R_Ring_Vector res_c(c1.first.Get_q(), c1.first.Get_d(), (c1.first.Get_Dimension() * (c1.first.Get_Dimension() + 1)) / 2);
	for (int i = 0; i < c1.first.Get_Dimension(); i++) {
	  res_c[i] = c1.first[i];
	}
	Pair<R_Ring_Vector, int> result(res_c, c1.second);
	assert(res_c.Get_q() == (*my_pk)[c1.second].Get_q());
	Refresh(result, my_pk);
	c1 = result;
	assert(result.first.Get_q() == (*my_pk)[c1.second].Get_q());
	assert(c1.first.Get_q() == (*my_pk)[c1.second].Get_q());
	//	Refresh(c1, my_pk);
      }
    }
  }
}

Pair<R_Ring_Vector, int> FHE_Cipher_Text::Mult(Pair<R_Ring_Vector, int> c, ZZ n, FHE_Secret_Key_Type *sk) {
  /**** debug code ****/
  // noise check
  ZZ B;
  if (sk != NULL) {
    B = c.first.Dot_Product((*sk)[c.second]).Get_Norm();
    assert(c.first.Get_Field_Expansion() * B * n * 2 < c.first.Get_q());
  }
  /**** endof debug code ****/

  int L = ((*my_pk).size() - 1) / 2;
  R_Ring_Vector res_c(c.first.Get_q(), c.first.Get_d(), (c.first.Get_Dimension() * (c.first.Get_Dimension() + 1)) / 2);
  for (int i = 0; i < c.first.Get_Dimension(); i++) {
    res_c[i] = c.first[i] * n;
  }
  
  /**** debug code ****/
  if (sk != NULL) {
    B = res_c.Dot_Product((*sk)[c.second].Tensor_Product((*sk)[c.second])).Get_Norm();
    assert(B * 2 < c.first.Get_q());
  }
  /**** endof debug code ****/

  Pair<R_Ring_Vector, int> result(res_c, c.second);
  Refresh(result, my_pk, sk);
  return result;
}

Pair<R_Ring_Vector, int> FHE_Cipher_Text::Add(Pair<R_Ring_Vector, int> c1, Pair<R_Ring_Vector, int> c2, bool sign, FHE_Secret_Key_Type *sk) { // sign = true for addition, sign = false for substraction
  int L = ((*my_pk).size() - 1) / 2;
  Update_To_Same_Level(c1, c2);

  /**** debug code ****/
  ZZ B1, B2;
  if (sk != NULL) {
    B1 = c1.first.Dot_Product((*sk)[c1.second]).Get_Norm();
    B2 = c2.first.Dot_Product((*sk)[c1.second]).Get_Norm();
    assert((B1 + B2) * 2 < c1.first.Get_q());
  }
  /**** endof debug code ****/

  R_Ring_Vector c3_p = sign ? (c1.first + c2.first) : (c1.first - c2.first);

  /* // TODO: to figure out when do we need to make refresh
  assert(c1.first.Get_Dimension() == c2.first.Get_Dimension());
  R_Ring_Vector c3(c3_p.Get_q(), c3_p.Get_d(), (c1.first.Get_Dimension() * (c1.first.Get_Dimension() + 1)) / 2);
  for (int i = 0; i < c1.first.Get_Dimension(); i++) {
    c3[i] = c3_p[i];
  }
  assert(c1.second >= 1);
  // going from c1.second level to c1.second - 1 level
  Pair<R_Ring_Vector, int> result(c3, c1.second);
  Refresh(result, my_pk);
  return result; */
  return Pair<R_Ring_Vector, int>(c3_p, c1.second);
}

Pair<R_Ring_Vector, int> FHE_Cipher_Text::Mult(Pair<R_Ring_Vector, int> c1, Pair<R_Ring_Vector, int> c2, FHE_Secret_Key_Type *sk) {
  int L = ((*my_pk).size() - 1) / 2;
  Update_To_Same_Level(c1, c2);

  /**** debug code ****/
  ZZ B1, B2;
  if (sk != NULL) {
    B1 = c1.first.Dot_Product((*sk)[c1.second]).Get_Norm();
    B2 = c2.first.Dot_Product((*sk)[c1.second]).Get_Norm();
    assert(c1.first.Get_Field_Expansion() * B1 * B2 * 2 < c1.first.Get_q());
  }
  /**** endof debug code ****/

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

  /**** debug code ****/
  if (sk != NULL) {
    B1 = c3.Dot_Product((*sk)[c1.second].Tensor_Product((*sk)[c1.second])).Get_Norm();
    assert(B1 * 2 < c3.Get_q());
  }
  /**** endof debug code ****/

  // going from c1.second level to c1.second - 1 level
  Pair<R_Ring_Vector, int> result(c3, c1.second);
  Refresh(result, my_pk, sk);
  return result;
}
