#include "FHE.h"

// COMMENT: along with each cipher element, I keep it's maximum noise bound and
// recomute it with every operation

// The result is modul c1.Get_q(), to switch to other modul, use Scale
R_Ring_Vector FHE_Cipher_Text::Switch_Key(const R_Ring_Matrix &A, const R_Ring_Vector &c1) const {
  // assert(A.Get_q() == c1.Get_q());
  return FHE::Bit_Decomposition(c1, A.Get_q()) * A;
}

R_Ring_Vector FHE_Cipher_Text::Scale(R_Ring_Vector &x, ZZ q, ZZ p, ZZ r) {
  // assert(p < q);
  R_Ring_Vector rn(p, x.Get_d(), x.Get_Dimension());
  for (int i = 0; i < x.Get_Dimension(); i++) {
    rn[i] = x[i].Scale(q, p, r);
  }
  return rn;
}
  
// A_j_j_1 - matrix that transforms key from level j to level j - 1
// q_j_1 < q_j
// Moves c from level c.second to level c.second - 1
void FHE_Cipher_Text::Refresh(FHE_Cipher_Text &c_cipher) const {
  Pair<R_Ring_Vector, int> &c = c_cipher.my_cipher;
  FHE_Public_Key_Type *pk = c_cipher.my_pk;
  FHE_Secret_Key_Type *sk = c_cipher.my_sk;

  // assert(c.second - 1 >= 0);
  int L = ((*my_pk).size() - 1) / 2;
  R_Ring_Vector c1 = FHE::Powersof2(c.first, (*pk)[c.second].Get_q());

  ZZ ThNoise = c_cipher.ThNoise;
  // noise check
  // check upper bound is not exceeded
#ifdef _CHECK_NOISE_
  ZZ required_upper_bound = (*pk)[c.second].Get_q() / 2 - (*pk)[c.second - 1].Get_q() * c.first.Get_d() * (int)(c.first.Get_Field_Expansion() + 1) * (((*pk)[c.second].Get_Noof_Columns() - 1) * (*pk)[c.second].Get_Noof_Columns()) / 2 * NumBits(c.first.Get_q()) / (*pk)[c.second].Get_q();
  // assert(ThNoise < required_upper_bound);
  // calculate expected amount of noise after Scaling
  ThNoise = (*pk)[c.second - 1].Get_q() * ThNoise / (*pk)[c.second].Get_q() + c.first.Get_d() * (int)(c.first.Get_Field_Expansion() + 1) * (((*pk)[c.second].Get_Noof_Columns() - 1) * (*pk)[c.second].Get_Noof_Columns()) / 2 * NumBits(c.first.Get_q());
  // assert(ThNoise * 2 < c.first.Get_q());
#endif /* _CHECK_NOISE_ */

#ifdef _FHE_DEBUG_MODE_
  /**** debug code ****/
  // noise check
  ZZ B;
  R_Ring_Vector s_pp;
  if (sk != NULL) {
    s_pp = FHE::Bit_Decomposition((*sk)[c.second].Tensor_Product((*sk)[c.second]), (*sk)[c.second].Get_q());
    B = c1.Dot_Product(s_pp).Get_Norm();
    ZZ required_upper_bound = (*sk)[c.second].Get_q() / 2 - (*sk)[c.second - 1].Get_q() * c.first.Get_d() * (int)(c.first.Get_Field_Expansion() + 1) * (((*sk)[c.second].Get_Dimension() - 1) * (*sk)[c.second].Get_Dimension()) / 2 * NumBits(c.first.Get_q()) / (*sk)[c.second].Get_q();
    // assert(B < required_upper_bound);
    ZZ expected_B = (*sk)[c.second - 1].Get_q() * B / (*sk)[c.second].Get_q() + c.first.Get_d() * (int)(c.first.Get_Field_Expansion() + 1) * (((*sk)[c.second].Get_Dimension() - 1) * (*sk)[c.second].Get_Dimension()) / 2 * NumBits(c.first.Get_q());
    // assert(expected_B * 2 < c.first.Get_q());
  }
  /**** endof debug code ****/
#endif /* _FHE_DEBUG_MODE_ */
  
  R_Ring_Vector c2 = Scale(c1, (*pk)[c.second].Get_q(), (*pk)[c.second - 1].Get_q(), my_p);
  
#ifdef _FHE_DEBUG_MODE_
  /**** debug code ****/
  // noise check
  if (sk != NULL) {
    s_pp.Decrease_Modul(c2.Get_q());
    B = c2.Dot_Product(s_pp).Get_Norm();
    // assert(B * 2 < c2.Get_q());
  }
  /**** endof debug code ****/
#endif /* _FHE_DEBUG_MODE_ */
  
  // assert(c2.Get_q() == (*pk)[L + c.second].Get_q());
  R_Ring_Vector c3 = Switch_Key((*pk)[L + c.second], c2);

#ifdef _CHECK_NOISE_
  ThNoise = ThNoise + (int)(2 * sqrt((double)c.first.Get_d()) * 2) * c.first.Get_Field_Expansion() * (((*pk)[c.second].Get_Noof_Columns() - 1) * (*pk)[c.second].Get_Noof_Columns()) / 2 * NumBits(c.first.Get_q()) * NumBits(c.first.Get_q());
  // assert(ThNoise * 2 < c2.Get_q());
#endif /* _CHECK_NOISE_ */

#ifdef _FHE_DEBUG_MODE_
  /**** debug code ****/
  // noise check
  if (sk != NULL) {
    ZZ new_B = c3.Dot_Product((*sk)[c.second - 1]).Get_Norm();
    // 2 multiplier before second summand stands for B_{\xi}
    ZZ expected_B = B + (int)(2 * sqrt((double)c.first.Get_d()) * 2) * c.first.Get_Field_Expansion() * (((*sk)[c.second].Get_Dimension() - 1) * (*sk)[c.second].Get_Dimension()) / 2 * NumBits(c.first.Get_q()) * NumBits(c.first.Get_q());
    if (expected_B * 2 >= c3.Get_q()) {
      std::cout << "expected_B = " << expected_B << std::endl;
      std::cout << "expected_B * 2 = " << expected_B * 2 << std::endl;
      std::cout << "c3.Get_q() = " << c3.Get_q() << std::endl;
      // assert(expected_B * 2 < c3.Get_q());
    }
    // assert(new_B * 2 < c2.Get_q());
  }
  /**** endof debug code ****/
#endif /* _FHE_DEBUG_MODE_ */

  c_cipher.my_cipher = Pair<R_Ring_Vector, int>(c3, c.second - 1); // updating cipher
  c_cipher.ThNoise = ThNoise; // updating theoretical noise bound
}
  
void FHE_Cipher_Text::Update_To_Same_Level(FHE_Cipher_Text &c1_cipher, FHE_Cipher_Text &c2_cipher) const {
  Pair<R_Ring_Vector, int> &c1 = c1_cipher.my_cipher, &c2 = c2_cipher.my_cipher;
  if (c1.second == c2.second) {
    return;
  }
  if (c2.second < c1.second) {
    Update_To_Same_Level(c2_cipher, c1_cipher);
    return;
  }    
  assert(c1.second < c2.second);

  int L = ((*my_pk).size() - 1) / 2;
  while (c1.second != c2.second) {
    // assert(c2.second - 1 >= 0);
    // switch c2 to next level
    // c2 goes from level c2.second to level c2.second - 1
    R_Ring_Vector res_c(c2.first.Get_q(), c2.first.Get_d(), (c2.first.Get_Dimension() * (c2.first.Get_Dimension() + 1)) / 2);
    for (int i = 0; i < c2.first.Get_Dimension(); i++) {
      res_c[i] = c2.first[i];
    }
    Pair<R_Ring_Vector, int> result(res_c, c2.second);
    // assert(res_c.Get_q() == (*my_pk)[c2.second].Get_q());
    c2_cipher.my_cipher = result; // copying here
    Refresh(c2_cipher);
    // assert(c2_cipher.my_cipher.first.Get_q() == (*my_pk)[c2.second].Get_q());
    // assert(c2.first.Get_q() == (*my_pk)[c2.second].Get_q());
  }
}

FHE_Cipher_Text FHE_Cipher_Text::Mult(const FHE_Cipher_Text &c_cipher, ZZ n) const {
  const Pair<R_Ring_Vector, int> &c = c_cipher.my_cipher;
  ZZ ThNoise = c_cipher.ThNoise;
  FHE_Secret_Key_Type *sk = my_sk;
  FHE_Public_Key_Type *pk = my_pk;

#ifdef _CHECK_NOISE_
  ThNoise = c.first.Get_Field_Expansion() * ThNoise * n;
  // assert(ThNoise * 2 < c.first.Get_q());
#endif /* _CHECK_NOISE_ */
  
#ifdef _FHE_DEBUG_MODE_
  /**** debug code ****/
  // noise check
  ZZ B;
  if (sk != NULL) {
    B = c.first.Dot_Product((*sk)[c.second]).Get_Norm();
    // assert(c.first.Get_Field_Expansion() * B * n * 2 < c.first.Get_q());
  }
  /**** endof debug code ****/
#endif /* _FHE_DEBUG_MODE_ */

  int L = ((*my_pk).size() - 1) / 2;
  R_Ring_Vector res_c(c.first.Get_q(), c.first.Get_d(), (c.first.Get_Dimension() * (c.first.Get_Dimension() + 1)) / 2);
  for (int i = 0; i < c.first.Get_Dimension(); i++) {
    res_c[i] = c.first[i] * n;
  }
  
#ifdef _FHE_DEBUG_MODE_
  /**** debug code ****/
  if (sk != NULL) {
    B = res_c.Dot_Product((*sk)[c.second].Tensor_Product((*sk)[c.second])).Get_Norm();
    // assert(B * 2 < c.first.Get_q());
  }
  /**** endof debug code ****/
#endif /* _FHE_DEBUG_MODE_ */

  FHE_Cipher_Text result = FHE_Cipher_Text(Pair<R_Ring_Vector, int>(res_c, c.second), my_pk, my_p, ThNoise, my_sk);
  Refresh(result);
  return result;
}

// here should be copying, i.e. we do not want to change given ciphers
FHE_Cipher_Text FHE_Cipher_Text::Add(FHE_Cipher_Text &c1_cipher, FHE_Cipher_Text &c2_cipher, bool sign) const { // sign = true for addition, sign = false for substraction
  Pair<R_Ring_Vector, int> &c1 = c1_cipher.my_cipher;
  Pair<R_Ring_Vector, int> &c2 = c2_cipher.my_cipher;
  FHE_Secret_Key_Type *sk = my_sk;

  int L = ((*my_pk).size() - 1) / 2;
  Update_To_Same_Level(c1_cipher, c2_cipher);
  bool need_refresh = false; // TODO: deduce this value from theoretical estimations
#ifndef _CHECK_NOISE
  ZZ ThNoise = ZZ::zero();
#else
  ZZ ThNoise = c1_cipher.ThNoise + c2_cipher.ThNoise;

  if (ThNoise * 2 >= c1.first.Get_q()) {
    std::cout << "c1.ThNoise = " << c1_cipher.ThNoise << std::endl;
    std::cout << "c2.ThNoise = " << c2_cipher.ThNoise << std::endl;
    std::cout << "ThNoise = " << ThNoise << std::endl;
    std::cout << c1.first.Get_q() << std::endl;
  }
  // assert(ThNoise * 2 < c1.first.Get_q()); // TODO: to think how to handle this event
#endif /* _CHECK_NOISE_ */

#ifdef _FHE_DEBUG_MODE_
  /**** debug code ****/
  ZZ B1, B2;
  if (sk != NULL) {
    B1 = c1.first.Dot_Product((*sk)[c1.second]).Get_Norm();
    B2 = c2.first.Dot_Product((*sk)[c1.second]).Get_Norm();
    if ((B1 + B2) * 2 >= c1.first.Get_q()) {
      std::cout << "B1 = " << B1 << std::endl;
      std::cout << "B2 = " << B2 << std::endl;
      std::cout << "(B1 + B2) * 2 = " << (B1 + B2) * 2 << std::endl;
      std::cout << "(B1 + B2) = " << (B1 + B2) << std::endl;
      std::cout << "c1.first.Get_q() = " << c1.first.Get_q() << std::endl;
      need_refresh = true;
    }
    //    // assert((B1 + B2) * 2 < c1.first.Get_q());
  }
  /**** endof debug code ****/
#endif _FHE_DEBUG_MODE_

  R_Ring_Vector c3_p = sign ? (c1.first + c2.first) : (c1.first - c2.first);
  
  if (need_refresh) {
    // assert(c1.first.Get_Dimension() == c2.first.Get_Dimension());
    R_Ring_Vector c3(c3_p.Get_q(), c3_p.Get_d(), (c1.first.Get_Dimension() * (c1.first.Get_Dimension() + 1)) / 2);
    for (int i = 0; i < c1.first.Get_Dimension(); i++) {
      c3[i] = c3_p[i];
    }
    // assert(c1.second >= 1);
    // going from c1.second level to c1.second - 1 level
    FHE_Cipher_Text result = FHE_Cipher_Text(Pair<R_Ring_Vector, int>(c3, c1.second), my_pk, my_p, ThNoise, my_sk);
    Refresh(result);

#ifdef _FHE_DEBUG_MODE_
    /**** debug code ****/
    if (sk != NULL) {
      ZZ B = result.my_cipher.first.Dot_Product((*sk)[result.my_cipher.second].Tensor_Product((*sk)[result.my_cipher.second])).Get_Norm();
      // assert(B * 2 < result.my_cipher.first.Get_q());
    }
    /**** endof debug code ****/
#endif _FHE_DEBUG_MODE_
    return result;
  }
  return FHE_Cipher_Text(Pair<R_Ring_Vector, int>(c3_p, c1.second), my_pk, my_p, ThNoise, my_sk);
}

FHE_Cipher_Text FHE_Cipher_Text::Mult(FHE_Cipher_Text &c1_cipher, FHE_Cipher_Text &c2_cipher) const {
  FHE_Secret_Key_Type *sk = my_sk;
  Pair<R_Ring_Vector, int> &c1 = c1_cipher.my_cipher, &c2 = c2_cipher.my_cipher;
  
  int L = ((*my_pk).size() - 1) / 2;
  Update_To_Same_Level(c1_cipher, c2_cipher);

#ifndef _CHECK_NOISE_
  ZZ ThNoise = ZZ::zero();
#else
  ZZ ThNoise = c1_cipher.ThNoise * c2_cipher.ThNoise * c1.first.Get_Field_Expansion();
  // assert(ThNoise * 2 < c1.first.Get_q());
#endif /* _CHECK_NOISE_ */

#ifdef _FHE_DEBUG_MODE_
  /**** debug code ****/
  ZZ B1, B2;
  if (sk != NULL) {
    B1 = c1.first.Dot_Product((*sk)[c1.second]).Get_Norm();
    B2 = c2.first.Dot_Product((*sk)[c1.second]).Get_Norm();
    // assert(c1.first.Get_Field_Expansion() * B1 * B2 * 2 < c1.first.Get_q());
  }
  /**** endof debug code ****/
#endif _FHE_DEBUG_MODE_

  // assert(c1.second == c2.second);
  // assert(c1.first.Get_Dimension() == c2.first.Get_Dimension());
  // assert(c1.first.Get_q() == c2.first.Get_q());
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
  // assert(index == new_dimension);
  // assert(c1.second >= 1);

#ifdef _FHE_DEBUG_MODE_
  /**** debug code ****/
  if (sk != NULL) {
    B1 = c3.Dot_Product((*sk)[c1.second].Tensor_Product((*sk)[c1.second])).Get_Norm();
    // assert(B1 * 2 < c3.Get_q());
  }
  /**** endof debug code ****/
#endif _FHE_DEBUG_MODE_

  // going from c1.second level to c1.second - 1 level
  FHE_Cipher_Text result = FHE_Cipher_Text(Pair<R_Ring_Vector, int>(c3, c1.second), my_pk, my_p, ThNoise, my_sk);
  Refresh(result);
  return result;
}
