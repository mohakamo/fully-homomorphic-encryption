/**
 *  SI_HE.h
 *  Scale Invariant Homomorphic Encryption Scheme
 *  Created by valerini on 1/5/12.
 */
#ifndef _SI_HE_H_
#define _SI_HE_H_

#include "SI_HE_Cipher_Text.h"
#define SQ(x) ((x) * (x))

class SI_HE {
 private:  
  R_Ring_Matrix Switch_Key_Gen(const R_Ring_Vector &s_from, const R_Ring_Vector &s_to, Regev_Params &params) const {
    assert(s_from.Get_Dimension() == SQ((params.n + 1) * NumBits(params.q)));
    int N = s_from.Get_Dimension() * NumBits(params.q);
    R_Ring_Matrix result(params.q, params.d, N, s_to.Get_Dimension() + 1);
    R_Ring_Matrix A_prime = R_Ring_Matrix::Uniform_Rand(params.q, params.d, N, params.n);
    R_Ring_Vector ksi_noise(params.q, params.d, N);
    for (int i = 0; i < N; i++) {
      ksi_noise[i] = params.ksi();
    }

    R_Ring_Vector s_from_powers;
    Powersof2(s_from, params.q, s_from_powers);

    result.Set_Column(0, A_prime * s_to + ksi_noise + s_from_powers);
    result.Set_Block(0, 1, -A_prime);
    return result;
  }
  
  Regev E;
  int my_L;
 public:

  /**
   * Setting up parameters
   * @param L			maximum depth of the circuit that the scheme can evaluate
   * @param b			b == LWE_Based or b == RLWE_Based
   * @param p                   module for message representation (should be coprime with all the modules in the ladder)
   * @return			set of parameters for each level of the circuit
   **/
  Regev_Params Setup(int L, GLWE_Type b, ZZ p = ZZ(INIT_VAL, 2), int n = -1) {
    my_L = L;
    Regev_Params params = E.Setup(b, p, n);
    if (NumBits(params.q) * NumBits(params.q) * (params.n + 1) >= params.q / 4) {
      std::cout << "log(q) * log(q) * (n + 1) = " << NumBits(params.q) * NumBits(params.q) * (params.n + 1) << std::endl;
      std::cout << "q / 4 = " << params.q / 4 << std::endl;
      //    assert(NumBits(params.q) * NumBits(params.q) * (params.n + 1) < params.q / 4);
    }
    return params;
  }
  
  SI_HE_Secret_Key_Type Secret_Key_Gen(Regev_Params &params) const {
    SI_HE_Secret_Key_Type sk(my_L + 1);
    for (int i = 0; i <= my_L; i++) {
      sk[i] = E.Secret_Key_Gen(params);
    }
    return sk; 
  }

  SI_HE_Public_Key_Type Public_Key_Gen(Regev_Params &params, SI_HE_Secret_Key_Type &sk) const {
    return E.Public_Key_Gen(params, sk[0]);
  }

  SI_HE_Evaluation_Key_Type Evaluation_Key_Gen(Regev_Params &params,
					       SI_HE_Secret_Key_Type &sk,
					       SI_HE_Public_Key_Type &pk) const {
    std::vector<R_Ring_Vector> sk_decomposed(my_L);
    SI_HE_Evaluation_Key_Type evalk(my_L);
    assert(sk.size() == my_L + 1);
    for (int i = 1; i <= my_L; i++) {
      /* Computing sk_2 := BitDecomp((1, sk_{i - 1}) \tensor_prod BitDecomp((1, sk_{i - 1}) */
      R_Ring_Vector sk1(sk[i - 1].Get_q(), sk[i - 1].Get_d(), sk[i - 1].Get_Dimension() + 1);
      sk1[0] = 1;
      for (int j = 1; j <= sk[i - 1].Get_Dimension(); j++) {
	sk1[j] = sk[i - 1][j - 1];
      }
      R_Ring_Vector sk2;
      Bit_Decomposition(sk1, params.q, sk2);
      sk2 = sk2.Tensor_Product2(sk2);

      evalk[i - 1] = Switch_Key_Gen(sk2, sk[i], params);
    }
    return evalk;
  }

  // For params should be a zero level of parameters generated by Setup, i.e. Setup(...)[0]
  SI_HE_Cipher_Text Encrypt(Regev_Params &params, SI_HE_Public_Key_Type *pk, R_Ring_Number m,
			    SI_HE_Evaluation_Key_Type *evalk) {
    ZZ ThNoise = ZZ::zero();
    R_Ring_Vector vec = E.Encrypt(params, *pk, m);
    return SI_HE_Cipher_Text(vec, 0, pk, evalk, params.p, ThNoise);
  }
  
  // For params should be a zero level of parameters generated by Setup, i.e. Setup(...)[0]
  SI_HE_Cipher_Text Encrypt(Regev_Params &params, SI_HE_Public_Key_Type *pk, ZZ m) {
    R_Ring_Number mp(params.p, params.d);
    mp = m;
    return Encrypt(params, pk, mp);
  }

  friend class SI_HE_Cipher_Text;
  friend class Regev;
};

#endif /* _SI_HE_H_ */
