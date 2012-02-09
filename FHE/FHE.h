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
  FHE_Secret_Key_Type *my_sk;
  int my_p;

  R_Ring_Vector Switch_Key(R_Ring_Matrix A, R_Ring_Vector c1);
 public: // for testing purposes
  static R_Ring_Vector Scale(R_Ring_Vector &x, ZZ q, ZZ p, int r);
  void Update_To_Same_Level(Pair<R_Ring_Vector, int> &c1, Pair<R_Ring_Vector, int> &c2);
 private:
  Pair<R_Ring_Vector, int> Add(Pair<R_Ring_Vector, int> c1, Pair<R_Ring_Vector, int> c2, bool sign = true, FHE_Secret_Key_Type *sk = NULL);
  Pair<R_Ring_Vector, int> Mult(Pair<R_Ring_Vector, int> c1, Pair<R_Ring_Vector, int> c2, FHE_Secret_Key_Type *sk = NULL);
  Pair<R_Ring_Vector, int> Mult(Pair<R_Ring_Vector, int> c1, ZZ n, FHE_Secret_Key_Type *sk = NULL);
 public:
 FHE_Cipher_Text(const Pair<R_Ring_Vector, int> &cipher, FHE_Public_Key_Type *pk, int p = 2, FHE_Secret_Key_Type *sk = NULL) :
  my_cipher(Pair<R_Ring_Vector, int>(R_Ring_Vector(cipher.first), cipher.second)), my_pk(pk), my_p(p), my_sk(sk) {}
 FHE_Cipher_Text(const FHE_Cipher_Text &c) : my_cipher(R_Ring_Vector(c.my_cipher.first), c.my_cipher.second), my_pk(c.my_pk), my_p(c.my_p), my_sk(c.my_sk) {}
 FHE_Cipher_Text() : my_cipher(R_Ring_Vector(), 0), my_pk(NULL), my_p(0), my_sk(NULL) {}

 FHE_Cipher_Text& operator =(const FHE_Cipher_Text &c) {
    my_cipher = Pair<R_Ring_Vector, int>(R_Ring_Vector(c.my_cipher.first), c.my_cipher.second);
    my_pk = c.my_pk;
    my_p = c.my_p;
    my_sk = c.my_sk;
    return *this;
 }
  
  FHE_Cipher_Text operator +(FHE_Cipher_Text &c) {
    assert(my_pk == c.my_pk); // addresses comparison
    return FHE_Cipher_Text(Add(my_cipher, c.my_cipher, true, my_sk), my_pk, my_p, my_sk);
  }

  FHE_Cipher_Text operator -(FHE_Cipher_Text &c) {
    assert(my_pk == c.my_pk);
    return FHE_Cipher_Text(Add(my_cipher, c.my_cipher, false, my_sk), my_pk, my_p, my_sk);
  }
  
  FHE_Cipher_Text operator -() { // assuming the field operations keep numbers in range [-(q - 1) / 2; (q - 1) / 2]
    FHE_Cipher_Text result = *this;
    for (int i = 0; i < result.my_cipher.first.Get_Dimension(); i++) {
      result.my_cipher.first[i] = -result.my_cipher.first[i];
    }
    return result;
  }

  FHE_Cipher_Text operator *(ZZ n) {
    return FHE_Cipher_Text(Mult(my_cipher, n, my_sk), my_pk, my_p, my_sk);
  }
  
  FHE_Cipher_Text operator *(FHE_Cipher_Text &c) {
    assert(my_pk == c.my_pk); // addresses, comparison
    return FHE_Cipher_Text(Mult(my_cipher, c.my_cipher, my_sk), my_pk, my_p, my_sk);
  }

  Pair<R_Ring_Vector, int> Copy_Cipher(void) const {
    R_Ring_Vector new_vector = my_cipher.first;
    return Pair<R_Ring_Vector, int>(new_vector, my_cipher.second);
  }

  R_Ring_Number Decrypt(FHE_Params &params, FHE_Secret_Key_Type &sk) {
    int j = my_cipher.second;
    return GLWE::Decrypt(params[j], sk[j], my_cipher.first);
  }

  void Refresh(Pair<R_Ring_Vector, int> &c, FHE_Public_Key_Type *pk, FHE_Secret_Key_Type *sk = NULL);
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
  Pair<R_Ring_Vector, int>& Get_Cipher() {
    return my_cipher;
  }

  void Add_Secret_Key_Info(FHE_Secret_Key_Type *sk) {
    my_sk = sk;
  }
};

class FHE {
 public: // for testing purposes
  static R_Ring_Vector Bit_Decomposition(R_Ring_Vector x, ZZ q) {
    assert(q == x.Get_q());
    int noof_vectors = NumBits(q);
    R_Ring_Vector res_r(q, x.Get_d(), noof_vectors * x.Get_Dimension());
    for (int i = 0; i < x.Get_Dimension(); i++) {
      for (int p = 0; p < noof_vectors; p++) {
	for (int j = 0; j < x.Get_d(); j++) {
	  /* if (x[i][j] < 0) {
	    std::cout << "x = "; x.print(); std::cout << std::endl;
	    exit(1);
	    } */
	  ZZ n = R_Ring_Number::Clamp(x[i][j], q);
	  n = n + (n < 0 ? q : ZZ::zero());
	  //	  assert(n >= 0 && n < q);
	  res_r[p + i * noof_vectors][j] = (n >> p) & 1; 
	}
      }
    }
    return res_r;
  }
  
  static R_Ring_Vector Powersof2(R_Ring_Vector x, ZZ q) {
    int noof_vectors = NumBits(q);
    R_Ring_Vector res_r(q, x.Get_d(), noof_vectors * x.Get_Dimension());
    int index;
    
    for (int i = 0; i < x.Get_Dimension(); i++) {
      res_r[i * noof_vectors] = x[i];
      for (int p = 1; p < noof_vectors; p++) {
	index = p + i * noof_vectors;
	res_r[index] = res_r[index - 1] * ZZ(INIT_VAL, 2); // not to have overflowing
      }
    }
    return res_r;
  }
 private:  
  R_Ring_Matrix Switch_Key_Gen(R_Ring_Vector s1, R_Ring_Vector s2, GLWE_Params params2) {
    // s1 should be in the bit representation, i.e. bounded moduli 2
    // assert(s1.Get_q() == 2);
    int N = s1.Get_Dimension() * NumBits(s2.Get_q());
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

  ZZ Choose_Noise_Bound(ZZ p, FHE_Params &params) { // suppose that all n are the same (for LWE scheme)
    ZZ max_fraq_num, max_fraq_den;
    max_fraq_num = max_fraq_den = -1;
    int d = params[params.size() - 1].d;
    int field_expansion = sqrt(d);
    int L = params.size() - 1;
    int n = params[L].n;
    for (int i = 1; i < params.size(); i++) {
      ZZ fraq_num, fraq_den;
      fraq_num = params[i].q;
      fraq_den = params[i - 1].q;
      if (max_fraq_num == -1 || fraq_num * max_fraq_den > max_fraq_num * fraq_den) {
	max_fraq_num = fraq_num;
	max_fraq_den = fraq_den;
      }
    }
    //    assert(max_fraq > 0);
    if (max_fraq_num <= 0) {return ZZ(INIT_VAL, -1);}
    ZZ noise_UB = (max_fraq_num / max_fraq_den / 2 / field_expansion - 1) / p / d / (2 * params[L].n + 1) / NumBits(params[L].q); // noise upper bound
    if (noise_UB <= 1) {return ZZ(INIT_VAL, -1);}
    ZZ noise_LB_temp;
    noise_LB_temp = (2 * d * (int)field_expansion * (n + 1) * n / 2 * NumBits(params[L].q) - 1) / d; // noise lower bound without denominator
    ZZ noise_LB = ZZ::zero(); // noise lower bound

    // computing denominators and trying to find lower bound
    for (int j = 1; j < params.size(); j++) {
      ZZ noise_LB_den = p * (2 * params[L].n + 1) * NumBits(params[L].q) - 4 * (params[L].n + 1) * params[L].n / 2 * NumBits(params[j].q) * NumBits(params[j - 1].q);
      ZZ potential_LB = noise_LB_temp / noise_LB_den;
      if (noise_LB_den > 0 && potential_LB > noise_LB) {
	noise_LB = potential_LB;
      } else if (noise_LB_den < 0 && -potential_LB < noise_UB) {
	noise_UB = -potential_LB;
      }
    }
    std::cout << "Theoretical noise bound = [" << noise_LB << ", " << noise_UB << "]" << std::endl;
    if (noise_LB > noise_UB || noise_UB <= 1) {
      return ZZ(INIT_VAL, -1);
    }
    return noise_LB + 1 > 1 ? noise_LB + 1 : ZZ(INIT_VAL, 2); // trying to make noise as small as possible
  }

  GLWE E;
  int my_L;
 public:
  void Print_Possible_Parameters(int L, GLWE_Type type, ZZ modul) {
    /*
    // inifinite cycle inside: interrupt manually
    bool found_noise_bound = false;
    ZZ x, y, q0, mu, q0_l, mu_l;
    q0_l = NumBits(modul) + 1;
    mu_l = 3;
    
    for (y = 0; ; y++) {
      for (x = 0; x <= y; x++) {
	mu = x + mu_l;
	q0 = y - x + q0_l;

	std::vector<GLWE_Params> try_params(L + 1);
	for (int i = 0; i <= L; i++) {
	  ZZ modul_size;
	  modul_size = q0 + i * mu;
	  try {
	    try_params[i] = E.Setup(100, modul_size, type, modul); // modules q are increasing, so while doing noise cleaning we need to switch from bigger module to smaller, i.e. from j to j - 1
	  } catch (bool b) {
	    std::cout << "Exception caught" << std::endl;
	    continue;
	  }
	}

	if (Choose_Noise_Bound(modul, try_params) > 0) {
	  found_noise_bound = true;
	  std::cout << "Parameters found, q0 = " << q0 << ", mu = " << mu << std::endl << std::endl;
	}
      }
      }*/
  }

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
    
    while (q_size < 8 * sizeof(ZZ) && ((1 << q_size) - 2) / q_size < 2 * p * d * (2 * n + 1)) {
      q_size++;
    }
    q_size += 2;
    //    q_size += 4;
    mu += 9;
    q_size = 10;
    mu = 12;

    /*** Choosing modules step ***/
    /*
    bool found_noise_bound = false;
    for (int q0 = 5; q0 < 60; q0++) {
      //      std::cout << "q0 = " << q0 << std::endl;
      for (int mu = 3; L * mu + q0 < 61 && mu < q0; mu++) {

	std::vector<GLWE_Params> try_params(L + 1);
	for (int i = 0; i <= L; i++) {
	  int modul_size = q_size + i * mu;
	  try {
	    try_params[i] = E.Setup(lambda, modul_size, b, p); // modules q are increasing, so while doing noise cleaning we need to switch from bigger module to smaller, i.e. from j to j - 1
	  } catch (bool b) {
	    std::cout << "Exception caught" << std::endl;
	    continue;
	  }
	}

	if (Choose_Noise_Bound(ZZ(INIT_VAL, p), try_params) > 0) {
	  found_noise_bound = true;
	  std::cout << "Parameters found, q0 = " << q0 << ", mu = " << mu << std::endl;
	}
      }
    }
    if (!found_noise_bound) {
      std::cout << "Noise bound not found" << std::endl;
      std::cout << "Message modul = " << p << std::endl;
      assert(false);
    }*/
    q_size = 23;
    mu = 16;
    int B = 5;
    /*
    //    mu = 9;
    int mu_low_bound = mu;
    while ((q_size + (L - 1) * mu) * 2 >= sizeof(ZZ) * 8 - 1 && mu >= mu_low_bound) {
      mu--;
      std::cout << "!";
    }
    */
    
    std::vector<GLWE_Params> params(L + 1);
    for (int i = 0; i <= L; i++) {
      int modul_size = q_size + i * mu;
      //      if (modul_size * 2 >= sizeof(ZZ) * 8 - 1) {
      //	std::cout << "ALARM: moduls_size = " << modul_size << " > " << (sizeof(ZZ) * 8 - 1) / 2 << std::endl;
      //	assert(modul_size * 2 < sizeof(ZZ) * 8 - 1); // to make addition, we need module twice smaller
      //      }
      params[i] = E.Setup(lambda, modul_size, b, p); // modules q are increasing, so while doing noise cleaning we need to switch from bigger module to smaller, i.e. from j to j - 1
    }
    /*
    ZZ min_fraq = ceil(params[L].q / (double)params[L - 1].q);
    for (int i = L - 1; i >= 1; i--) {
      if (ceil(params[i].q / (double)params[i - 1].q) < min_fraq) {
	min_fraq = ceil(params[i].q / (double)params[i - 1].q);
      }
    }
    ZZ B = 1 + p * params[L].d * params[L].N * params[L].B;
    ZZ B_upper_bound = min_fraq / 2 / sqrt(1.0 * d);
    ZZ B_lower_bound = 2 * (sqrt(1.0 * d) * ((params[L].n * (params[L].n + 1)) / 2) * ceil(log2(params[L].q)) * (params[L].d + 2 * sqrt(1.0 * params[L].d) * params[L].B));
    ZZ Bx = 2;
    if (B < B_lower_bound) {
      std::cout << "Adjusting noise " << std::endl;
      std::cout << "params[L].d = " << params[L].d << ", params[L].N = " << params[L].N << ", params[L].B = " << params[L].B << std::endl;
      while (B < B_lower_bound && Bx < params[0].q) {
	std::cout << "B = " << B << " not in [" << B_lower_bound << ", " << B_upper_bound << "], Bx = " << Bx << std::endl;
	B_lower_bound = 2 * (sqrt(1.0 * d) * ((params[L].n * (params[L].n + 1)) / 2) * ceil(log2(params[L].q)) * (params[L].d + 2 * sqrt(1.0 * params[L].d) * Bx));
	B = 1 + p * params[L].d * params[L].N * Bx;
	if (B_lower_bound < B_upper_bound && B < B_upper_bound && B > B_lower_bound) {
	  break;
	}
	if (B_lower_bound > B_upper_bound) {
	  break;
	}
	Bx++;
      }
    }
    
    if (B > B_upper_bound || B < B_lower_bound) {
      std::cout << "B = " << B << " not int [" << B_lower_bound << ", " << B_upper_bound << "]" << std::endl;
      std::cout << "min_fraq = " << min_fraq << std::endl;
      std::cout << "n = " << params[L].n << ", log(q) = " << ceil(log2(params[L].q)) << ", d = " << params[L].d << ", B = " << params[0].B << std::endl;
      exit(1);
    }
    //    int B = params[0].q - 2);
    */
			     
    for (int i = 0; i < L; i++) {
      params[i].d = params[L].d;
      // params[i].noise = params[L].noise;
      params[L - i].B = B;//params[0].B;
    }
    // print noise bounds parameters and assert if the interval is empty
    Choose_Noise_Bound(ZZ(INIT_VAL, p), params);
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
	assert(s_i_prime_prime.Get_Dimension() == s_i_prime.Get_Dimension() * NumBits(params[i + 1].q));
	R_Ring_Matrix tau_i = Switch_Key_Gen(s_i_prime_prime, s_i, params[i]); // give the bigger modul
	assert(tau_i.Get_Noof_Columns() == s_i.Get_Dimension());
	assert(s_i.Get_Dimension() == params[i].n + 1);
	assert(tau_i.Get_Noof_Rows() == s_i_prime_prime.Get_Dimension() * NumBits(params[i].q));
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
