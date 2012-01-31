/*
 *  FHE.h
 *  Created by valerini on 1/5/12.
 */
#ifndef _FHE_H_
#define _FHE_H_

#include "GLWE.h"
#include "Pair.h"

typedef std::vector<R_Ring_Matrix> FHE_Public_Key_Type;
typedef std::vector<R_Ring_Vector> FHE_Secret_Key_Type;
typedef std::vector<GLWE_Params> FHE_Params;
typedef Pair<R_Ring_Vector, int> FHE_Ciphertext_Type;

class FHE_Cipher_Text {
  Pair<R_Ring_Vector, int> my_cipher;
  FHE_Public_Key_Type *my_pk;
  int my_p;

  R_Ring_Vector Switch_Key(R_Ring_Matrix A, R_Ring_Vector c1);
 public: // for testing purposes
  static R_Ring_Vector Scale(R_Ring_Vector &x, int q, int p, int r);
 private:
  void Update_To_Same_Level(Pair<R_Ring_Vector, int> &c1, Pair<R_Ring_Vector, int> &c2);
  Pair<R_Ring_Vector, int> Add(Pair<R_Ring_Vector, int> &c1, Pair<R_Ring_Vector, int> &c2, bool sign = true);
  Pair<R_Ring_Vector, int> Mult(Pair<R_Ring_Vector, int> &c1, Pair<R_Ring_Vector, int> &c2);
  Pair<R_Ring_Vector, int> Mult(Pair<R_Ring_Vector, int> &c1, int n);
 public:
 FHE_Cipher_Text(const Pair<R_Ring_Vector, int> &cipher, FHE_Public_Key_Type *pk, int p = 2) : my_cipher(Pair<R_Ring_Vector, int>(R_Ring_Vector(cipher.first), cipher.second)), my_pk(pk), my_p(p) {}
 FHE_Cipher_Text(const FHE_Cipher_Text &c) : my_cipher(R_Ring_Vector(c.my_cipher.first), c.my_cipher.second), my_pk(c.my_pk), my_p(c.my_p) {}
 FHE_Cipher_Text() : my_cipher(R_Ring_Vector(), 0), my_pk(NULL), my_p(0) {}

  /*  FHE_Cipher_Text& operator =(const FHE_Cipher_Text &c) {
    my_cipher = Pair<R_Ring_Vector, int>(R_Ring_Vector(c.my_cipher.first), c.my_cipher.second);
    my_pk = c.my_pk;
    return *this;
    }*/
  
  FHE_Cipher_Text operator +(FHE_Cipher_Text &c) {
    assert(my_pk == c.my_pk); // addresses comparison
    return FHE_Cipher_Text(Add(my_cipher, c.my_cipher), my_pk, my_p);
  }

  FHE_Cipher_Text operator -(FHE_Cipher_Text &c) {
    assert(my_pk == c.my_pk);
    return FHE_Cipher_Text(Add(my_cipher, c.my_cipher, false), my_pk, my_p);
  }

  FHE_Cipher_Text operator *(int n) {
    return FHE_Cipher_Text(Mult(my_cipher, n), my_pk, my_p);
  }
  
  FHE_Cipher_Text operator *(FHE_Cipher_Text &c) {
    assert(my_pk == c.my_pk); // addresses, comparison
    return FHE_Cipher_Text(Mult(my_cipher, c.my_cipher), my_pk, my_p);
  }

  Pair<R_Ring_Vector, int> Copy_Cipher(void) const {
    R_Ring_Vector new_vector = my_cipher.first;
    return Pair<R_Ring_Vector, int>(new_vector, my_cipher.second);
  }

  R_Ring_Number Decrypt(FHE_Params &params, FHE_Secret_Key_Type &sk) {
    int j = my_cipher.second;
    return GLWE::Decrypt(params[j], sk[j], my_cipher.first);
  }

  void Refresh(Pair<R_Ring_Vector, int> &c, FHE_Public_Key_Type *pk);
  void print(void) {
    std::cout << "(";
    my_cipher.first.print();
    std::cout << ", " << my_cipher.second << ")";
  }
  // debugging functions
  R_Ring_Number Get_Noise(FHE_Params &params, FHE_Secret_Key_Type &sk) {
    int j = my_cipher.second;
    return GLWE::Get_Noise(params[j], sk[j], my_cipher.first);
  }
  Pair<R_Ring_Vector, int> Get_Cipher() {
    return my_cipher;
  }
};

class FHE {
 public: // for testing purposes
  static R_Ring_Vector Bit_Decomposition(R_Ring_Vector x, int q) {
    assert(q == x.Get_q());
    int noof_vectors = ceil(log2(q));
    R_Ring_Vector res_r(q, x.Get_d(), noof_vectors * x.Get_Dimension());
    for (int i = 0; i < x.Get_Dimension(); i++) {
      for (int p = 0; p < noof_vectors; p++) {
	for (int j = 0; j < x.Get_d(); j++) {
	  /* if (x[i][j] < 0) {
	    std::cout << "x = "; x.print(); std::cout << std::endl;
	    exit(1);
	    } */
	  int n = R_Ring_Number::Clamp(x[i][j], q);
	  n = n + (n < 0 ? q : 0);
	  //	  assert(n >= 0 && n < q);
	  res_r[p + i * noof_vectors][j] = (n >> p) & 1; 
	}
      }
    }
    return res_r;
  }
  
  static R_Ring_Vector Powersof2(R_Ring_Vector x, int q) {
    int noof_vectors = ceil(log2(q));
    R_Ring_Vector res_r(q, x.Get_d(), noof_vectors * x.Get_Dimension());
    int index;
    
    for (int i = 0; i < x.Get_Dimension(); i++) {
      res_r[i * noof_vectors] = x[i];
      for (int p = 1; p < noof_vectors; p++) {
	index = p + i * noof_vectors;
	res_r[index] = res_r[index - 1] * 2; // not to have overflowing
      }
    }
    return res_r;
  }
 private:  
  R_Ring_Matrix Switch_Key_Gen(R_Ring_Vector s1, R_Ring_Vector s2, GLWE_Params params2) {
    // s1 should be in the bit representation, i.e. bounded moduli 2
    // assert(s1.Get_q() == 2);
    int N = s1.Get_Dimension() * ceil(log2(s2.Get_q()));
    params2.N = N;
    R_Ring_Matrix A = E.Public_Key_Gen(params2, s2);
    assert(A.Get_Noof_Columns() == s2.Get_Dimension());
    assert(A.Get_Noof_Rows() == N);
    return A.Add_To_Column(0, Powersof2(s1, s2.Get_q()));
  }
  
  int Choose_mu(int lambda, int L) {
    // TODO: to be implemented
    return 6;
  }

  GLWE E;
  int my_L;
 public:
  /**
   * Setting up parameters
   * @param lambda	security parameter, 100 is a good value for it
   * @param L			maximum depth of the circuit that the scheme can evaluate
   * @param b			b == LWE_Based or b == RLWE_Based
   * @param p                   module for message representation (should be coprime with all the modules in the ladder)
   * @return			set of parameters for each level of the circuit
   **/
  FHE_Params Setup(int lambda, int L, GLWE_Type b, int p = 2) {
    my_L = L;
    // initial modul length
    int q_size = 5;
    int d = E.Choose_d(lambda, 0, b);
    int n = E.Choose_n(lambda, 0, b);
    int mu = ceil(4 * sqrt((double)d));
    
    while (q_size < 8 * sizeof(long long) && ((1 << q_size) - 2) / q_size < 4 * p * d * (2 * n + 1)) {
      q_size++;
    }
    // q_size++;
    q_size += 4;
    //    std::cout << "Chosen initial modul size = " << q_size << std::endl;
    
    std::vector<GLWE_Params> params(L + 1);
    for (int i = 0; i <= L; i++) {
      params[i] = E.Setup(lambda, q_size + i * mu, b, p); // modules q are increasing, so while doing noise cleaning we need to switch from bigger module to smaller, i.e. from j to j - 1
    }
    for (int i = 0; i < L; i++) {
      params[i].d = params[L].d;
      // params[i].noise = params[L].noise;
      params[L - i].B = params[0].B;
    }
    // std::cout << "B = " << params[0].B << std::endl;
    return params;
  }
  
  /**
   * Key generation function
   * @return Pair<FHE_Secret_Key, FHE_Public_Key>
   */
  Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> Key_Gen(std::vector<GLWE_Params> params) {
    std::vector<R_Ring_Vector> sk(my_L + 1); // s_0, ..., s_L
    std::vector<R_Ring_Matrix> pk(my_L + 1 + my_L); // first L + 1 slots are devoted to A_0, ..., A_L and next L slots are devoted to tau_(1 -> 0), ..., tau_(L -> L - 1)
    for (int i = my_L; i >= 0; i--) {
      R_Ring_Vector s_i = E.Secret_Key_Gen(params[i]);
      /* for (int j = 0; j < s_i.Get_Dimension(); j++) {
	for (int jj = 0; jj < s_i.Get_d(); jj++) {
	  if (s_i[j][jj] < 0) {
	    std::cout << "s_" << i << " = "; s_i.print(); std::cout << std::endl;
	    exit(1);
	  }
	}
	} */

      sk[i] = s_i;
      R_Ring_Matrix A_i = E.Public_Key_Gen(params[i], s_i);
      pk[i] = A_i;
      if (i != my_L) {
	/* for (int j = 0; j < sk[i + 1].Get_Dimension(); j++) {
	  for (int jj = 0; jj < sk[i + 1].Get_d(); jj++) {
	    if (sk[i + 1][j][jj] < 0) {
	      std::cout << "sk[i + 1]" << i << " = "; sk[i + 1].print(); std::cout << std::endl;
	      exit(1);
	    }
	  }
	  } */
	
	R_Ring_Vector s_i_prime = sk[i + 1].Tensor_Product(sk[i + 1]); // from R_{q_j}^{(n_j + 1, 2)} space, i.e. there are n_j * (n_j + 1) element
	assert(s_i_prime.Get_Dimension() == (sk[i+1].Get_Dimension() * (sk[i+1].Get_Dimension() + 1)) / 2);

	/* for (int j = 0; j < s_i_prime.Get_Dimension(); j++) {
	  for (int jj = 0; jj < s_i_prime.Get_d(); jj++) {
	    if (s_i_prime[j][jj] < 0) {
	      std::cout << "sk[i + 1] = "; sk[i + 1].print(); std::cout << std::endl;
	      std::cout << "s_i_prime" << " = "; s_i_prime.print(); std::cout << std::endl;
	      R_Ring_Vector s_i_prime2 = sk[i + 1].Tensor_Product(sk[i + 1]); // from R_{q_j}^{(n_j + 1, 2)} space, i.e. there are n_j * (n_j + 1) element
	      exit(1);
	    }
	  }
	  } */

	R_Ring_Vector s_i_prime_prime = Bit_Decomposition(s_i_prime, params[i + 1].q);
	assert(s_i_prime_prime.Get_Dimension() == s_i_prime.Get_Dimension() * ceil(log2(params[i + 1].q)));
	R_Ring_Matrix tau_i = Switch_Key_Gen(s_i_prime_prime, s_i, params[i]); // give the bigger modul
	assert(tau_i.Get_Noof_Columns() == s_i.Get_Dimension());
	assert(s_i.Get_Dimension() == params[i].n + 1);
	assert(tau_i.Get_Noof_Rows() == s_i_prime_prime.Get_Dimension() * ceil(log2(params[i].q)));
	pk[my_L + i + 1] = tau_i;
      }
    }
    return Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type>(sk, pk);
  }
  
  // For params should be a zero level of parameters generated by Setup, i.e. Setup(...)[0]
  /*
    Pair<R_Ring_Vector, int> Encrypt(FHE_Params params, FHE_Public_Key_Type pk, R_Ring_Number m) {
    return Pair<R_Ring_Vector, int> (E.Encrypt(params[L], pk[L], m), L); // pk[L] == A_L
    }
  */
  FHE_Cipher_Text Encrypt(FHE_Params params, FHE_Public_Key_Type *pk, R_Ring_Number m) {
    return FHE_Cipher_Text(Pair<R_Ring_Vector, int> (E.Encrypt(params[my_L], (*pk)[my_L], m), my_L), pk, params[0].p);
  }
  
  friend class FHE_Cipher_Text;
  friend class GLWE;
};

#endif /* _FHE_H_ */
