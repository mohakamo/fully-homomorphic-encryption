/**
 *  FHE.h
 *  Created by valerini on 1/5/12.
 */
#ifndef _FHE_H_
#define _FHE_H_

#include "Assert.h"
#include "GLWE.h"
#include "Pair.h"
#include "time.h"

typedef std::vector<R_Ring_Matrix> FHE_Public_Key_Type;
typedef std::vector<R_Ring_Vector> FHE_Secret_Key_Type;
typedef std::vector<GLWE_Params> FHE_Params;
typedef Pair<R_Ring_Vector, int> FHE_Ciphertext_Type;

class FHE_Cipher_Text {
  Pair<R_Ring_Vector, int> my_cipher;
  FHE_Public_Key_Type *my_pk;
  FHE_Secret_Key_Type *my_sk;
  ZZ my_p;

  R_Ring_Vector Switch_Key(const R_Ring_Matrix &A, const R_Ring_Vector &c1) const;
 public: // for testing purposes
  ZZ ThNoise;
  static R_Ring_Vector Scale(R_Ring_Vector &x, ZZ q, ZZ p, ZZ r);
  void Update_To_Same_Level(FHE_Cipher_Text &c1, FHE_Cipher_Text &c2) const;
 private:
  FHE_Cipher_Text Add(FHE_Cipher_Text &c1, FHE_Cipher_Text &c2, bool sign = true) const;
  FHE_Cipher_Text Mult(FHE_Cipher_Text &c1, FHE_Cipher_Text &c2) const;
  FHE_Cipher_Text Mult(const FHE_Cipher_Text &c1, ZZ n) const;
 public:
 FHE_Cipher_Text(const Pair<R_Ring_Vector, int> &cipher, FHE_Public_Key_Type *pk, ZZ p, ZZ noise, FHE_Secret_Key_Type *sk = NULL) :
  my_cipher(Pair<R_Ring_Vector, int>(R_Ring_Vector(cipher.first), cipher.second)), my_pk(pk), my_p(p), my_sk(sk), ThNoise(noise) {}
 FHE_Cipher_Text(const FHE_Cipher_Text &c) : my_cipher(R_Ring_Vector(c.my_cipher.first), c.my_cipher.second), my_pk(c.my_pk), my_p(c.my_p), my_sk(c.my_sk), ThNoise(c.ThNoise) {}
 FHE_Cipher_Text() : my_cipher(R_Ring_Vector(), 0), my_pk(NULL), my_sk(NULL) {
    my_p = ZZ::zero();
    ThNoise = ZZ::zero();
  }

 FHE_Cipher_Text& operator =(const FHE_Cipher_Text &c) {
    my_cipher = Pair<R_Ring_Vector, int>(R_Ring_Vector(c.my_cipher.first), c.my_cipher.second);
    my_pk = c.my_pk;
    my_p = c.my_p;
    my_sk = c.my_sk;
    ThNoise = c.ThNoise;
    return *this;
 }
  
  FHE_Cipher_Text operator +(FHE_Cipher_Text &c) {
    // assert(my_pk == c.my_pk); // addresses comparison
    // assert(ThNoise != 0);
    // assert(c.ThNoise != 0);
    return Add(*this, c, true);
  }

  FHE_Cipher_Text operator -(FHE_Cipher_Text &c) {
    // assert(my_pk == c.my_pk);
    // assert(ThNoise != 0);
    // assert(c.ThNoise != 0);
    return Add(*this, c, false);
  }
  
  FHE_Cipher_Text operator -() const { // assuming the field operations keep numbers in range [-(q - 1) / 2; (q - 1) / 2]
    FHE_Cipher_Text result = *this;
    for (int i = 0; i < result.my_cipher.first.Get_Dimension(); i++) {
      result.my_cipher.first[i] = -result.my_cipher.first[i];
    }
    return result;
  }

  FHE_Cipher_Text operator *(ZZ n) const {
    // assert(ThNoise != 0);
    return Mult(*this, n);
  }
  
  FHE_Cipher_Text operator *(FHE_Cipher_Text &c) {
    // assert(my_pk == c.my_pk); // addresses, comparison
    // assert(ThNoise != 0);
    // assert(c.ThNoise != 0);
    return Mult(*this, c);
  }

  Pair<R_Ring_Vector, int> Copy_Cipher(void) const {
    R_Ring_Vector new_vector = my_cipher.first;
    return Pair<R_Ring_Vector, int>(new_vector, my_cipher.second);
  }

  R_Ring_Number Decrypt(const FHE_Params &params, const FHE_Secret_Key_Type &sk) const {
    int j = my_cipher.second;
    return GLWE::Decrypt(params[j], sk[j], my_cipher.first);
  }

  void Refresh(FHE_Cipher_Text &c) const;
  void print(void) {
    std::cout << "(";
    my_cipher.first.print();
    std::cout << ", " << my_cipher.second << ")";
  }
  // debugging functions
  R_Ring_Number Get_Noise(FHE_Params &params, FHE_Secret_Key_Type &sk) const {
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
  static void Bit_Decomposition(const R_Ring_Vector &x, ZZ q, R_Ring_Vector &result) {
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
  
  static void Powersof2(const R_Ring_Vector &x, ZZ q, R_Ring_Vector &result) {
    static ZZ TWO = ZZ(INIT_VAL, 2);
    int noof_vectors = NumBits(q);
    if (result.Get_Dimension() != noof_vectors * x.Get_Dimension() ||
	result.Get_d() != x.Get_d() ||
	result.Get_q() != q) {
      result = R_Ring_Vector(q, x.Get_d(), noof_vectors * x.Get_Dimension());
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
 private:  
  R_Ring_Matrix Switch_Key_Gen(const R_Ring_Vector &s1, const R_Ring_Vector &s2, GLWE_Params &params2) {
    // s1 should be in the bit representation, i.e. bounded moduli 2
    // // assert(s1.Get_q() == 2);
    int N = s1.Get_Dimension() * NumBits(s2.Get_q());
    int Old_N = params2.N;
    params2.N = N;
    R_Ring_Matrix A = E.Public_Key_Gen(params2, s2);
    params2.N = Old_N;
    // assert(A.Get_Noof_Columns() == s2.Get_Dimension());
    // assert(A.Get_Noof_Rows() == N);
    R_Ring_Vector column;
    Powersof2(s1, s2.Get_q(), column);
    return A.Add_To_Column(0, column);
  }
  
  int Choose_mu(int lambda, int L) {
    // mu is chosen separately, this value is not used
    return 6;
  }

  ZZ Choose_Noise_Bound(ZZ p, FHE_Params &params, bool printInfo = false) { // suppose that all n are the same (for LWE scheme)
    int d = params[params.size() - 1].d;
    int field_expansion = sqrt(d);
    int L = params.size() - 1;
    int n = params[L].n;
    ZZ noise_UB_num = params[1].q;
    ZZ noise_UB_den = (params[0].q * 2 * field_expansion * p * d * (2 * params[L].n + 1) * NumBits(params[L].q)); // noise upper bound
    if (printInfo) {
      std::cout << field_expansion  << ", " << p << ", " << d << ", " << (2 * params[L].n + 1) << ", " << NumBits(params[L].q) << ", L = " << L << std::endl;
      //      std::cout << "noise_UB = " << noise_UB << std::endl;
    }
    if (noise_UB_num <= noise_UB_den) {
      //      std::cout << "noise_UB = " << noise_UB << std::endl;
      return ZZ(INIT_VAL, -1);
    }
    ZZ noise_LB_num = ZZ::zero(), noise_LB_den = ZZ(INIT_VAL, 1); // noise lower bound

    // computing denominators and trying to find lower bound
    ZZ a = d * p * (2 * n + 1) * NumBits(params[L].q);
    ZZ b = ZZ(INIT_VAL, 2 * d * (n + 1) * n);
    ZZ c = ZZ(INIT_VAL, d * (int)field_expansion * (n + 1) * n);
    for (int j = 1; j < params.size(); j++) {
      ZZ noise_LB_num_tmp = c * NumBits(params[j].q) - 1; // noise lower bound without denominator
      ZZ noise_LB_den_tmp = a - b * NumBits(params[j].q) * NumBits(params[j - 1].q);
      if (noise_LB_den_tmp == 0) {
	continue;
      }
      //      ZZ potential_LB = noise_LB_temp / noise_LB_den;
      if (noise_LB_den_tmp > 0 && noise_LB_num_tmp * noise_LB_den > noise_LB_num * noise_LB_den_tmp) {
	noise_LB_num = noise_LB_num_tmp;
	noise_LB_den = noise_LB_den_tmp;
	//	if (printInfo) 	std::cout << "noise_LB = " << noise_LB << std::endl;
      } else if (noise_LB_den_tmp < 0 && -noise_LB_num_tmp * noise_UB_den < noise_UB_num * noise_LB_den_tmp) {
	noise_UB_num = noise_LB_num_tmp;
	noise_UB_den = -noise_LB_den_tmp;

	//	if (printInfo) {
	//	  std::cout << "noise_UB = " << noise_UB << std::endl;
	//	  std::cout << "noise_LB_temp = " << noise_LB_temp << std::endl;
	//	  std::cout << "noise_LB_den = " << noise_LB_den << std::endl;
	//	  std::cout << "p = " << p << std::endl;
	//	}
      }
      if (noise_LB_num * noise_UB_den > noise_UB_num * noise_LB_den || noise_UB_num <= noise_UB_den) {
	return ZZ(INIT_VAL, -1);
      }
    }
    //    std::cout << "Theoretical noise bound = [" << noise_LB << ", " << noise_UB << "]" << std::endl;
    return noise_UB_num / noise_UB_den;
  }

  GLWE E;
  int my_L;
 public:
  void Print_Possible_Parameters(int L, GLWE_Type type, ZZ modul, FHE_Params &params, ZZ &B,
				 bool justPrint = true) {
    std::vector<int> brokenModules;
    // inifinite cycle inside: interrupt manually
    bool found_noise_bound = false;
    int x, y, q0, mu, q0_l, mu_l;
    static int x0 = 0, y0 = 0;
    ZZ max_noise = ZZ::zero();
    q0_l = NumBits(modul) + 1;// * 2; // Why do I did |* 2 here???
    mu_l = NumBits(modul);
    if (justPrint) std::cout << "q0_l = " << q0_l << std::endl;
    //    x0 = 22, y0 = 22;
    std::vector<int> modules(L + 1);
    ZZ B_desired = ZZ(INIT_VAL, 10);
    
    for (y = y0; ; y++) {
      //      std::cout << y << std::endl;
      for (y == y0 ? x = x0 : x = 0; x <= y; x++) {
	//	std::cout << "  " << x << std::endl;
	mu = x + mu_l;
	q0 = y - x + q0_l;
	int q = q0;
	for (int i = 0; i <= L; i++) {
	  bool modul_found = false;
	  while (!modul_found) {
	    if (std::find(brokenModules.begin(), brokenModules.end(), q) != brokenModules.end()) {
	      q++;
	    } else {
	      try {
		GLWE::Choose_q(q, modul);
		modules[i] = q;
		modul_found = true;
	      } catch (bool t) {
		brokenModules.push_back(q);
		q++;
	      }
	    }
	  }
	  q += mu;
	}
	//	for (int nd = 2; nd < (q0 + L * mu - 1000) * 100; nd++)
	//	std::cout << "(" << q0 << ", " << mu << "), (" << x << ", " << y << ")" << std::endl;

	for (int i = 0; i <= L; i++) {
	  params[i] = E.Setup(100, modules[i], type, modul); // modules q are increasing, so while doing noise cleaning we need to switch from bigger module to smaller, i.e. from j to j - 1
	}

	ZZ Initial_noise = 1 + modul * params[0].d * (2 * params[0].n + 1) * NumBits(params[L].q) * B_desired;
	if (Initial_noise * 2 > params[0].q) {
	  continue;
	}

	ZZ noise_bound;
	noise_bound = Choose_Noise_Bound(modul, params);

	if (noise_bound > B_desired) {
	  found_noise_bound = true;
	  std::cout << "(x, y) = (" << x << ", " << y << ")" << std::endl;
	  std::cout << "Theoretical upper noise bound = [" <<  noise_bound << "]" << std::endl;
	  std::cout << "Parameters found, q0 = " << q0 << ", mu = " << mu << std::endl << std::endl;
	  
	  x0 = x;
	  y0 = y;

	  B = B_desired;
	  return;
	}
      }
    }
  }

  /**
   * Setting up parameters
   * @param lambda	security parameter, 100 is a good value for it
   * @param L			maximum depth of the circuit that the scheme can evaluate
   * @param b			b == LWE_Based or b == RLWE_Based
   * @param p                   module for message representation (should be coprime with all the modules in the ladder)
   * @return			set of parameters for each level of the circuit
   **/
  FHE_Params Setup(int lambda, int L, GLWE_Type b, ZZ p = ZZ(INIT_VAL, 2)) {
    my_L = L;
    // initial modul length
    clock_t start = clock();
    ZZ B;

    std::cout << "|p| = " << NumBits(p) << std::endl;

    std::vector<GLWE_Params> params(L + 1);
    Print_Possible_Parameters(L, b, p, params, B, false);

    std::cout << "Time for finding parameters = " << (clock() - start) / (double)CLOCKS_PER_SEC << std::endl;
    
    for (int i = 0; i <= L; i++) {
      std::cout << "|q" << i << "| = " << NumBits(params[i].q) << " ";
      params[i].d = params[L].d;
      params[L - i].B = B;
    }
    std::cout << std::endl;
    // print noise bounds parameters and ASSERT if the interval is empty
    // Choose_Noise_Bound(ZZ(INIT_VAL, p), params);

    return params;
  }
  
  /**
   * Key generation function
   * @return Pair<FHE_Secret_Key, FHE_Public_Key>
   */
  Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> Key_Gen(std::vector<GLWE_Params> params) {
    clock_t start = clock();
    std::vector<R_Ring_Vector> sk(my_L + 1); // s_0, ..., s_L
    std::vector<R_Ring_Matrix> pk(my_L + 1 + my_L); // first L + 1 slots are devoted to A_0, ..., A_L and next L slots are devoted to tau_(1 -> 0), ..., tau_(L -> L - 1)
    for (int i = my_L; i >= 0; i--) {
      R_Ring_Vector s_i = E.Secret_Key_Gen(params[i]);

      sk[i] = s_i;
      R_Ring_Matrix A_i = E.Public_Key_Gen(params[i], s_i);
      pk[i] = A_i;
      if (i != my_L) {
	R_Ring_Vector s_i_prime = sk[i + 1].Tensor_Product(sk[i + 1]); // from R_{q_j}^{(n_j + 1, 2)} space, i.e. there are n_j * (n_j + 1) element
	// assert(s_i_prime.Get_Dimension() == (sk[i+1].Get_Dimension() * (sk[i+1].Get_Dimension() + 1)) / 2);

	R_Ring_Vector s_i_prime_prime;
	Bit_Decomposition(s_i_prime, params[i + 1].q, s_i_prime_prime);
	// assert(s_i_prime_prime.Get_Dimension() == s_i_prime.Get_Dimension() * NumBits(params[i + 1].q));
	R_Ring_Matrix tau_i = Switch_Key_Gen(s_i_prime_prime, s_i, params[i]); // give the bigger modul
	// assert(tau_i.Get_Noof_Columns() == s_i.Get_Dimension());
	// assert(s_i.Get_Dimension() == params[i].n + 1);
	// assert(tau_i.Get_Noof_Rows() == s_i_prime_prime.Get_Dimension() * NumBits(params[i].q));
	pk[my_L + i + 1] = tau_i;
      }
    }
    std::cout << "Time for key generation = " << (clock() - start) / (double)CLOCKS_PER_SEC << std::endl;
    return Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type>(sk, pk);
  }
  
  // For params should be a zero level of parameters generated by Setup, i.e. Setup(...)[0]
  FHE_Cipher_Text Encrypt(FHE_Params &params, FHE_Public_Key_Type *pk, R_Ring_Number m) {
    ZZ ThNoise;
    ThNoise = 1 + params[my_L].p * params[my_L].d * sqrt(params[my_L].d) * params[my_L].N * params[my_L].B;
    return FHE_Cipher_Text(Pair<R_Ring_Vector, int> (E.Encrypt(params[my_L], (*pk)[my_L], m), my_L), pk, params[0].p, ThNoise);
  }
  
  // For params should be a zero level of parameters generated by Setup, i.e. Setup(...)[0]
  FHE_Cipher_Text Encrypt(FHE_Params &params, FHE_Public_Key_Type *pk, ZZ m) {
    R_Ring_Number mp(params[0].p, params[0].d);
    mp = m;
    return Encrypt(params, pk, mp);
  }

  friend class FHE_Cipher_Text;
  friend class GLWE;
};

#endif /* _FHE_H_ */
