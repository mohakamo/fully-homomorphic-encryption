/*
 *  GLWE.h
 *  Created by valerini on 1/5/12.
 */
#ifndef _GLWE_H_
#define _GLWE_H_
#include "NormDistr.h"
#include "R_Ring_Matrix.h"
#include <math.h>
#include <map>

/*** Params for GLWE_Encrption_Scheme ***/
class GLWE_Params {
public:
  GLWE_Params() {
    d = n = N = 0;
    p = B = q = ZZ::zero();
    noise = NULL;
  }
  GLWE_Params(ZZ q, int d, int n, int N, ZZ B, ZZ p,  R_Ring_Number (*noise)(ZZ q, int d, ZZ B, ZZ p)) {
    this->q = q;
    this->d = d;
    this->n = n;
    this->N = N;
    this->B = B;
    this->p = p;
    this->noise = noise;
  }
  ZZ B, q, p;
  int d, n, N;
  /*** Noise distribution ***/
  R_Ring_Number (*noise)(ZZ q, int d, ZZ B, ZZ p);
  
  R_Ring_Number ksi() const {
    return noise(q, d, B, p);
  }
  
  GLWE_Params(const GLWE_Params &r) {
    q = r.q;
    d = r.d;
    n = r.n;
    N = r.N;
    B = r.B;
    p = r.p;
    noise = r.noise;
  }

  void print() {
    std::cout << "(q=" << q << ", ";
    std::cout << "d=" << d << ", ";
    std::cout << "n=" << n << ", ";
    std::cout << "N=" << N << ", ";
    std::cout << "B=" << B << ",";
    std::cout << "p=" << p << ")";
  }
};

enum GLWE_Type {LWE_Based = 1, RLWE_Based};
typedef R_Ring_Vector GLWE_Secret_Key_Type;
typedef R_Ring_Matrix GLWE_Public_Key_Type;
typedef R_Ring_Vector GLWE_Ciphertext_Type;

/*** Basic GLWE_Encryption_Scheme ***/
class GLWE {
 public:
  // q is chosen to be prime of length mu such that q = 1 (mod r), where r is the message module
  static ZZ Choose_q(int mu, ZZ r) {
    static std::map<int, ZZ> primes;
    static ZZ my_r = r;

    if (my_r != r) {
      my_r = r;
      primes.clear();
    }

    std::map<int, ZZ>::iterator found_value = primes.find(mu);
    
    if (found_value != primes.end()) {
      return (*found_value).second;
    }
      
    ZZ n;
    n = (ZZ(INIT_VAL, 1) << mu) - 1;
    while (n > (ZZ(INIT_VAL, 1) << (mu - 1))) {
      if (ProbPrime(n) && n % r == 1) {
	primes.insert(std::pair<int, ZZ>(mu, n));
	return n;
      }
      n--;
    }
    //    std::cout << "mu = " << mu << std::endl;
    //    std::cout << "r = " << r << std::endl;
    // did not find modul that is equal to 1 mod r
    //    std::cout << "Exception, between " << (ZZ(INIT_VAL, 1) << (mu - 1)) << " and " << (ZZ(INIT_VAL, 1) << mu) - 1 << " did not prime find number that = 1 mod " << r << std::endl;
    throw false; // if modul was not found
    return ZZ(INIT_VAL, -1);
  }
 public:
  static int Choose_d(int lambda, int mu, GLWE_Type b) {
    if (b == LWE_Based) {
      return 1;
    }
    // TODO: to be implemented, seems like a reasonable value for noise >= 2^10
    return 8;
  }

  static int Choose_n(int lambda, int mu, GLWE_Type b) {
    if (b == RLWE_Based) {
      return 1;
    }
    // TODO: to be implemented, seems like a reasonable value for noise >= 2^10
    return 8;
  }
 private:
  static int Choose_N(int n, ZZ q) {
    return ceil((2 * n + 1) * NumBits(q));
  }

  // just return the uniform noise that is bounded by B (although should be gaussian - not uniform for security)
  static R_Ring_Number Noise(ZZ q, int d, ZZ B, ZZ p) {
    assert(B >= ZZ(INIT_VAL, 2));
    //    R_Ring_Number res = R_Ring_Number::Uniform_Rand(B, d);
    R_Ring_Number res = NormDistr::sample(q, d, 8);
    //    R_Ring_Number res = R_Ring_Number::Uniform_Rand(ZZ(INIT_VAL, 2), d);
    res.Increase_Modul(q);
    //    R_Ring_Number res(q, d);
    return res;
  }

  static ZZ Choose_B(ZZ q, int d, int N, ZZ p) {
    return ZZ(INIT_VAL, 2); // this value is not used for FHE, the noise is being chosen manually afterwards
  }
public:	
  GLWE_Params Setup(int lambda, int mu, GLWE_Type b, ZZ p = ZZ(INIT_VAL, 2)) const {
    ZZ q;
    q = Choose_q(mu, p);
    int n = Choose_n(lambda, mu, b);
    int N = Choose_N(n, q);
    int d = Choose_d(lambda, mu, b);
    ZZ B = Choose_B(q, d, N, p);

    return GLWE_Params(q, d, n, N, B, p, &Noise);
  }
	
  R_Ring_Vector Secret_Key_Gen(GLWE_Params &params) const {
    R_Ring_Vector sk(params.q, params.d, params.n + 1);
    sk[0] = 1;
    for (int i = 1; i < params.n + 1; i++) {
      sk[i] = params.ksi();
    }
    return sk; 
  }
  
  /*** Matrix enumeration goes as follows (remainder):
       ------------------>
       |
       |		A
       |
       |
       v
  ***/
  R_Ring_Matrix Public_Key_Gen(const GLWE_Params &params, const R_Ring_Vector &sk, R_Ring_Vector *ksi_noise_for_debug = NULL) const {
    /*
    assert(params.q == sk.Get_q());
    // A_prime
    R_Ring_Matrix A_prime(params.q, params.d, params.N, params.n);
    for (int i = 0; i < A_prime.Get_Noof_Rows(); i++) {
      for (int j = 0; j < A_prime.Get_Noof_Columns(); j++) {
	A_prime(i, j) = R_Ring_Number::Uniform_Rand(params.q, params.d);
      }
    }

    R_Ring_Vector ksi_noise(params.q, params.d, params.N);
    for (int i = 0; i < params.N; i++) {
      ksi_noise[i] = params.ksi();
    }

    R_Ring_Matrix pk(params.q, params.d, params.N, params.n + 1);
    pk.Set_Column(0, A_prime * sk.Get_Sub_Vector(1, params.n) + ksi_noise * params.p);
    pk.Set_Block(0, 1, -A_prime);

    if (ksi_noise_for_debug != NULL) {
      (*ksi_noise_for_debug) = ksi_noise;
    }
    */
    //    assert(sk.Get_Dimension() == params.n + 1);
    //    assert(s_prime.Get_Dimension() == params.n);
    //    assert(As.Get_Dimension() == params.N);
    // Old version
    assert(params.q == sk.Get_q());
    // A_prime
    R_Ring_Matrix A_prime(params.q, params.d, params.N, params.n);
    R_Ring_Number sample(params.q, params.d);
    for (int i = 0; i < A_prime.Get_Noof_Rows(); i++) {
      for (int j = 0; j < A_prime.Get_Noof_Columns(); j++) {
        A_prime(i, j) = sample.Uniform_Rand(params.q, params.d);
      }
    }

    R_Ring_Matrix minus_A_prime = -A_prime;
    // s_prime
    assert(sk.Get_Dimension() == params.n + 1);
    R_Ring_Vector s_prime = sk.Get_Sub_Vector(1, params.n);
    assert(s_prime.Get_Dimension() == params.n);

    // ksi_noise
    // R_Ring_Vector ksi_noise(params.q, params.d, params.N);
    R_Ring_Vector ksi_noise(params.q, params.d, params.N);
    for (int i = 0; i < params.N; i++) {
      ksi_noise[i] = params.ksi();
    }

    // pk
    R_Ring_Matrix pk(params.q, params.d, params.N, params.n + 1);
    R_Ring_Vector As = A_prime * s_prime;
    assert(As.Get_Dimension() == params.N);
    pk.Set_Column(0, As + ksi_noise * params.p);
    pk.Set_Block(0, 1, -A_prime);

    if (ksi_noise_for_debug != NULL) {
      (*ksi_noise_for_debug) = ksi_noise;
    }

    return pk;
  }
  
  R_Ring_Vector Encrypt(GLWE_Params &params, R_Ring_Matrix &pk, R_Ring_Number &m, R_Ring_Vector *r_for_debug = NULL) const {
    if (m.Get_q() != params.p) {
      std::cout << m.Get_q() << " " << params.p << std::endl;
    }
    assert(m.Get_q() == params.p);
    
    // m_prime
    R_Ring_Vector m_prime(params.q, params.d, params.n + 1);
    m_prime[0] = R_Ring_Number(params.q, params.d, m.vec);
    for (int i = 1; i < m_prime.Get_Dimension(); i++) {
      m_prime[i] = 0;
    }

    // r
    R_Ring_Vector r(params.q, params.d, params.N);
    for (int i = 0; i < r.Get_Dimension(); i++) {
      r[i] = R_Ring_Number(params.q, params.d, R_Ring_Number::Uniform_Rand(ZZ(INIT_VAL, 2), params.d).vec);
    }

    
    R_Ring_Vector c(params.q, params.d, params.n + 1);
    R_Ring_Matrix pk_transpose = pk.Get_Transpose();
    
    c = m_prime + (pk.Get_Transpose() * r);

    if (r_for_debug != NULL) {
      (*r_for_debug) = r;
    }
    return c;
  }

  // Debug purpose function
  static R_Ring_Number Get_Noise(GLWE_Params &params, R_Ring_Vector &sk, const R_Ring_Vector &c) {
    assert(c.Get_q() == params.q && sk.Get_q() == params.q);
    return c.Dot_Product(sk);
  }
  
  static R_Ring_Number Decrypt(const GLWE_Params &params, const R_Ring_Vector &sk, const R_Ring_Vector &c) {
    assert(c.Get_q() == params.q && sk.Get_q() == params.q);
    R_Ring_Number dot_product = c.Dot_Product(sk);
    R_Ring_Number clamped1 = dot_product.Clamp(params.q);
    R_Ring_Number clamped2 = clamped1.Clamp(params.p);
    return clamped2;
  }
};

#endif /* _GLWE_H_ */

