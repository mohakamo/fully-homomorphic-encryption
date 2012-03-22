#include "SI_HE_Cipher_Text.h"

SI_HE_Cipher_Text SI_HE_Cipher_Text::Add(SI_HE_Cipher_Text &c1,
					 SI_HE_Cipher_Text &c2,
						bool sign) { // sign = true for addition, sign = false for substraction
  // assuming that c1_cipher and c2_cipher are from the same encryption scheme,
  // meaning that pk, sk and evalk are the same for them

  Update_To_Same_Level(c1, c2);
  R_Ring_Vector sum = sign ? (c1.my_cipher + c2.my_cipher) : (c1.my_cipher - c2.my_cipher);
  return SI_HE_Cipher_Text(sum, c1.my_level, c1.my_pk, c1.my_eval, c1.my_p, ZZ::zero(), c1.my_sk);
}

SI_HE_Cipher_Text SI_HE_Cipher_Text::Mult(const SI_HE_Cipher_Text &c, ZZ n) {
  R_Ring_Vector new_cipher = c.my_cipher * n;

  return SI_HE_Cipher_Text(new_cipher, c.my_level, c.my_pk, c.my_eval, c.my_p, ZZ::zero(), c.my_sk);
}

SI_HE_Cipher_Text SI_HE_Cipher_Text::Mult(SI_HE_Cipher_Text &c1, SI_HE_Cipher_Text &c2) {
  Update_To_Same_Level(c1, c2);

  R_Ring_Vector c1_powers, c2_powers;
  Powersof2(c1.my_cipher, c1.my_cipher.Get_q(), c1_powers);
  Powersof2(c2.my_cipher, c2.my_cipher.Get_q(), c2_powers);

  // TODO: to compute all of this once - not on every mult
  ZZ q = c1.my_cipher.Get_q();
  ZZ q_half = q / 2;
  ZZ q_quater = q / 4;
  int tens_dimension = (c1_powers.Get_Dimension() * (c1_powers.Get_Dimension() + 1)) / 2;
  R_Ring_Vector q_quater_v(q, c1_powers.Get_d(), tens_dimension);
  q_quater_v[0] = q_quater;
  for (int i = 0; i < tens_dimension; i++) {
    q_quater_v[i] = q_quater;
  }

  /*  if (c1_powers.Get_q() != q_quater_v.Get_q()) {
    std::cout << "c1_powers.Get_q() =  " << c1_powers.Get_q() << std::endl;
    std::cout << "q_quater_v.Get_q() =  " << q_quater_v.Get_q() << std::endl;
    }*/

  c1_powers = (c1_powers + q_quater_v) / q_half;
  R_Ring_Vector tens = c1_powers.Tensor_Product(c2_powers);
  R_Ring_Vector c_mult = tens;//(tens + q_quater_v) / q_half;

  R_Ring_Vector new_cipher = Switch_Key((*c1.my_eval)[c1.my_level], c_mult);
  return SI_HE_Cipher_Text(new_cipher,
			   c1.my_level + 1, c1.my_pk, c1.my_eval, c1.my_p, ZZ::zero(), c1.my_sk);
}

void SI_HE_Cipher_Text::Raise_Level() {
  int new_dimension = (my_cipher.Get_Dimension() * (my_cipher.Get_Dimension() + 1)) / 2;
  R_Ring_Vector expanded_cipher(my_cipher.Get_q(), my_cipher.Get_d(), new_dimension);
  for (int i = 0; i < my_cipher.Get_Dimension(); i++) {
    expanded_cipher[i] = my_cipher[i];
  }
  my_cipher = Switch_Key((*my_eval)[my_level], expanded_cipher);
  my_level++;
}

// The result is modul c1.Get_q(), to switch to other modul, use Scale
R_Ring_Vector SI_HE_Cipher_Text::Switch_Key(const R_Ring_Matrix &A, const R_Ring_Vector &c) {
  R_Ring_Vector res;
  Bit_Decomposition(c, A.Get_q(), res);
  return res * A;
}

void SI_HE_Cipher_Text::Update_To_Same_Level(SI_HE_Cipher_Text &c1, SI_HE_Cipher_Text &c2) {
  if (c1.my_level == c2.my_level) {
    return;
  }
  if (c1.my_level > c2.my_level) {
    Update_To_Same_Level(c2, c1);
    return;
  }
  R_Ring_Vector one_v(c1.my_cipher.Get_q(), c1.my_cipher.Get_d(), c1.my_cipher.Get_Dimension());
  one_v[0] = 1;
  R_Ring_Vector one_v_powers;
  Powersof2(one_v, c1.my_cipher.Get_q(), one_v_powers);
  while (c1.my_level < c2.my_level) {
    R_Ring_Vector c_new = c1.my_cipher.Tensor_Product(one_v_powers);
    c1.my_cipher = Switch_Key((*c1.my_eval)[c1.my_level], c1.my_cipher);
    c1.my_level++;
  }
}

ostream& operator<<(ostream& s, const SI_HE_Cipher_Text& a) {
  s << "(" << a.my_cipher << ", " << a.my_level << ")";
  return s;
}
