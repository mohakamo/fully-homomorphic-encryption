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
  ZZ q_half = q / c1.my_p;
  ZZ q_quater = q / (c1.my_p * 2);
  int tens_dimension = c1_powers.Get_Dimension();

  /*  if (c1_powers.Get_q() != q_quater_v.Get_q()) {
    std::cout << "c1_powers.Get_q() =  " << c1_powers.Get_q() << std::endl;
    std::cout << "q_quater_v.Get_q() =  " << q_quater_v.Get_q() << std::endl;
    }*/

  // c1_powers = (c1_powers + q_quater_v) / q_half;
  // NB!
  // get tensor product without reduction (second parameter - false)
  R_Ring_Vector tens = c1_powers.Tensor_Product2(c2_powers, false);
  // in division operator reduction will be made
  R_Ring_Vector q_quater_v(q, c1_powers.Get_d(), tens.Get_Dimension());
  q_quater_v[0] = q_quater;
  for (int i = 0; i < tens.Get_Dimension(); i++) {
    q_quater_v[i] = q_quater;
  }

  R_Ring_Vector c_mult = R_Ring_Vector::sum(tens, q_quater_v, false) / q_half;

  /*** DEBUG CODE ***/
  R_Ring_Vector sk_prime = (*c1.my_sk)[c1.my_level];
  R_Ring_Vector sk(sk_prime.Get_q(), sk_prime.Get_d(), sk_prime.Get_Dimension() + 1);
  for (int i = 1; i < sk.Get_Dimension(); i++) {
    sk[i] = sk_prime[i - 1];
  }
  sk[0] = 1;
  
  R_Ring_Vector sk_bd;
  Bit_Decomposition(sk, sk.Get_q(), sk_bd);
  sk_bd = sk_bd.Tensor_Product2(sk_bd);

  //  std::cout << "noise1 = " << c1.my_cipher.Dot_Product(sk) << "\n";
  //  std::cout << "noise2 = " << c2.my_cipher.Dot_Product(sk) << "\n";
  //  std::cout << "noise after tensoring = " << tens.Dot_Product(sk_bd) << "\n";
  //  std::cout << "noise after division = " << c_mult.Dot_Product(sk_bd) << "\n";
  /*** END OF DEBUG CODE ***/

  R_Ring_Vector new_cipher = Switch_Key((*c1.my_eval)[c1.my_level], c_mult);
  /*** DEBUG CODE ***/
  R_Ring_Vector new_sk_prime = (*c1.my_sk)[c1.my_level + 1];
  R_Ring_Vector new_sk(new_sk_prime.Get_q(), new_sk_prime.Get_d(), new_sk_prime.Get_Dimension() + 1);
  new_sk[0] = 1;
  for (int i = 0; i < new_sk_prime.Get_Dimension(); i++) {
    new_sk[i + 1] = new_sk_prime[i];
  }
  //  std::cout << "noise after switch keys = " << new_cipher.Dot_Product(new_sk) << "\n";
  /*** END OF DEBUG CODE ***/
  return SI_HE_Cipher_Text(new_cipher,
			   c1.my_level + 1, c1.my_pk, c1.my_eval, c1.my_p, ZZ::zero(), c1.my_sk);
}

void SI_HE_Cipher_Text::Raise_Level() {
  R_Ring_Vector powers_c;
  //  assert(my_cipher.Get_Dimension() == 
  Powersof2(my_cipher, my_cipher.Get_q(), powers_c);
  R_Ring_Vector uno(my_cipher.Get_q(), my_cipher.Get_d(), powers_c.Get_Dimension());
  uno[0] = 1;
  uno = powers_c.Tensor_Product2(uno);
  
  /*  ZZ q = my_cipher.Get_q();
  assert(uno.Get_Dimension() == (my_cipher.Get_Dimension() * NumBits(q)) *
	 my_cipher.Get_Dimension() * NumBits(q));
  std::cout << "cipher_length = " << my_cipher.Get_Dimension() << "\n";
  std::cout << "expected = " << NumBits(q) * NumBits(q) << "\n";
  */

  my_cipher = Switch_Key((*my_eval)[my_level], uno);
  my_level++;
}

R_Ring_Vector SI_HE_Cipher_Text::Switch_Key(const R_Ring_Matrix &A, const R_Ring_Vector &c) {
  R_Ring_Vector res;
  Bit_Decomposition(c, A.Get_q(), res);
  if (res.Get_Dimension() != A.Get_Noof_Rows()) {
    ZZ q = c.Get_q();
    int n = A.Get_Noof_Columns() - 1;
    std::cout << "q = " << q << ", n = " << n << "\n";
    int expected_noof_columns = (n + 1) * (n + 1) * NumBits(q) * NumBits(q) * NumBits(q);
    std::cout << "expected_noof_columns = " << expected_noof_columns << "\n";
    std::cout << "res.Get_Dimension() = " << res.Get_Dimension() << "\n";
    std::cout << "A:  = " << A.Get_Noof_Rows() << " x " << A.Get_Noof_Columns() << "\n";
    assert(res.Get_Dimension() == A.Get_Noof_Rows());
  }
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
  //  R_Ring_Vector one_v(c1.my_cipher.Get_q(), c1.my_cipher.Get_d(), c1.my_cipher.Get_Dimension());
  //  one_v[0] = 1;
  //  R_Ring_Vector one_v_powers;
  //  Powersof2(one_v, c1.my_cipher.Get_q(), one_v_powers);
  while (c1.my_level < c2.my_level) {
    c1.Raise_Level();
    //    R_Ring_Vector c_new = c1.my_cipher.Tensor_Product2(one_v_powers);
    //    c1.my_cipher = Switch_Key((*c1.my_eval)[c1.my_level], c1.my_cipher);
    //    c1.my_level++;
  }
}

ostream& operator<<(ostream& s, const SI_HE_Cipher_Text& a) {
  s << "(" << a.my_cipher << ", " << a.my_level << ")";
  return s;
}
