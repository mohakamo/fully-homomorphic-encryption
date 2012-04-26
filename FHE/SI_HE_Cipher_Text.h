/**
 *  SI_HE_Cipher_Text.h
 *  Scale Invariant Homomorphic Encryption Cipher Text Class
 *  Created by valerini on 1/5/12.
 */
#ifndef _SI_HE_CIPHER_TEXT_H_
#define _SI_HE_CIPHER_TEXT_H_

#include "Assert.h"
#include "Regev.h"
#include "Pair.h"
#include "time.h"

typedef R_Ring_Matrix SI_HE_Public_Key_Type;
typedef std::vector<R_Ring_Vector> SI_HE_Secret_Key_Type;
typedef std::vector<R_Ring_Matrix> SI_HE_Evaluation_Key_Type;

class SI_HE_Cipher_Text {
  R_Ring_Vector my_cipher;
  int my_level;
  SI_HE_Public_Key_Type *my_pk;
  SI_HE_Secret_Key_Type *my_sk;
  SI_HE_Evaluation_Key_Type *my_eval;
  ZZ my_p; // message modul
  ZZ ThNoise;
  
 private:
  /**
   * Adding two ciphertexts together and returning the result
   * @param c1_cipher
   * @param c2_cipher
   * @param sign true for + and false for -
   * @return new ciphertext that is the sum of c1_cipher and c2_cipher
   **/
  static SI_HE_Cipher_Text Add(SI_HE_Cipher_Text &c1, SI_HE_Cipher_Text &c2, bool sign = true);
  static SI_HE_Cipher_Text Mult(SI_HE_Cipher_Text &c1, SI_HE_Cipher_Text &c2);
  static SI_HE_Cipher_Text Mult(const SI_HE_Cipher_Text &c1, ZZ n);
 public:
  static void Update_To_Same_Level(SI_HE_Cipher_Text &c1, SI_HE_Cipher_Text &c2);
  void Raise_Level();
 private:
  static R_Ring_Vector Switch_Key(const R_Ring_Matrix &A, const R_Ring_Vector &c);

 public:
 SI_HE_Cipher_Text(R_Ring_Vector &cipher,
		   int level,
		   SI_HE_Public_Key_Type *pk,
		   SI_HE_Evaluation_Key_Type *eval,
		   ZZ p,
		   ZZ noise,
		   SI_HE_Secret_Key_Type *sk = NULL) :
  my_cipher(cipher),
    my_level(level),
    my_pk(pk),
    my_p(p),
    my_sk(sk),
    ThNoise(noise),
    my_eval(eval) {}

 SI_HE_Cipher_Text(const SI_HE_Cipher_Text &c) :
  my_cipher(c.my_cipher),
  my_level(c.my_level),
    my_pk(c.my_pk),
    my_p(c.my_p),
    my_sk(c.my_sk),
    ThNoise(c.ThNoise),
    my_eval(c.my_eval) {}

 SI_HE_Cipher_Text() : my_cipher(R_Ring_Vector()), my_level(0), my_pk(NULL), my_sk(NULL), my_eval(NULL) {
    my_p = ZZ::zero();
    ThNoise = ZZ::zero();
  }

 SI_HE_Cipher_Text& operator =(const SI_HE_Cipher_Text &c) {
   my_cipher = c.my_cipher;
   my_level = c.my_level;
   my_pk = c.my_pk;
   my_p = c.my_p;
   my_sk = c.my_sk;
   my_eval = c.my_eval;
   ThNoise = c.ThNoise;
   return *this;
 }
  
  SI_HE_Cipher_Text operator +(SI_HE_Cipher_Text &c) {
    return Add(*this, c, true);
  }

  SI_HE_Cipher_Text operator -(SI_HE_Cipher_Text &c) {
    return Add(*this, c, false);
  }
  
  SI_HE_Cipher_Text operator -() const { // assuming the field operations keep numbers in range [-(q - 1) / 2; (q - 1) / 2]
    SI_HE_Cipher_Text result = *this;
    for (int i = 0; i < result.my_cipher.Get_Dimension(); i++) {
      result.my_cipher[i] = -result.my_cipher[i];
    }
    return result;
  }

  SI_HE_Cipher_Text operator *(ZZ n) const {
    return Mult(*this, n);
  }
  
  SI_HE_Cipher_Text operator *(SI_HE_Cipher_Text &c) {
    return Mult(*this, c);
  }

  R_Ring_Number Decrypt(const Regev_Params &params, const SI_HE_Secret_Key_Type &sk) const {
    int j = my_level;
    return Regev::Decrypt(params, sk[j], my_cipher);
  }

  void print(void) {
    std::cout << "(";
    my_cipher.print();
    std::cout << ", " << my_level << ")";
  }
  // debugging functions
  R_Ring_Number Get_Noise(Regev_Params &params, SI_HE_Secret_Key_Type &sk) const {
    int j = my_level;
    return Regev::Get_Noise(params, sk[j], my_cipher);
  }

  Pair<R_Ring_Vector, int> Get_Cipher() {
    return Pair<R_Ring_Vector, int>(my_cipher, my_level);
  }

  R_Ring_Vector Get_Cipher_Copy() {
    R_Ring_Vector cipher_clone = my_cipher;
    return cipher_clone;
  }

  void Add_Secret_Key_Info(SI_HE_Secret_Key_Type *sk) {
    my_sk = sk;
  }
 public:
	long long getSizeInBytes(void) {
		return my_cipher.Get_Dimension() * NumBits(my_cipher[0].Get_q()) * my_cipher[0].Get_d() / 8;
	}

  friend ostream& operator<<(ostream& s, const SI_HE_Cipher_Text& a);
};

ostream& operator<<(ostream& s, const SI_HE_Cipher_Text& a);

#endif /* _SI_HE_CIPHER_TEXT_H_ */
