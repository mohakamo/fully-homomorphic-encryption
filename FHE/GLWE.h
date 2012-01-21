/*
 *  GLWE.h
 *  Created by valerini on 1/5/12.
 */
#ifndef _GLWE_H_
#define _GLWE_H_
#include "R_Ring_Matrix.h"
#include <math.h>
#include "LOG.h"

/*** Params for GLWE_Encrption_Scheme ***/
class GLWE_Params {
public:
  GLWE_Params() {
    q = d = n = N = B = 0;
    noise = NULL;
  }
  GLWE_Params(int q, int d, int n, int N, int B,  R_Ring_Number (*noise)(int q, int d, int B)) {
    this->q = q;
    this->d = d;
    this->n = n;
    this->N = N;
    this->B = B;
    this->noise = noise;
  }
  int q, d, n, N, B;
  /*** Noise distribution ***/
  R_Ring_Number (*noise)(int q, int d, int B);
  
  R_Ring_Number ksi() {
    return noise(q, d, B);
  }
  
  GLWE_Params(const GLWE_Params &r) {
    q = r.q;
    d = r.d;
    n = r.n;
    N = r.N;
    B = r.B;
    noise = r.noise;
  }

  void print() {
    std::cout << "q = " << q << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "n = " << n << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "B = " << B << std::endl;
  }
};

enum GLWE_Type {LWE_Based = 1, RLWE_Based};
typedef R_Ring_Vector GLWE_Secret_Key_Type;
typedef R_Ring_Matrix GLWE_Public_Key_Type;
typedef R_Ring_Vector GLWE_Ciphertext_Type;

/*** Basic GLWE_Encryption_Scheme ***/
class GLWE {
  static int Choose_q(int mu) {
    assert(mu < sizeof(int) * 8);
    return (1 << mu) - 1;
    // return 4095;
  }
 public:
  static int Choose_d(int lambda, int mu, GLWE_Type b) {
    if (b == LWE_Based) {
      return 1;
    }
    // TODO: to be implemented
    return 2;
  }

  static int Choose_n(int lambda, int mu, GLWE_Type b) {
    if (b == RLWE_Based) {
      return 1;
    }
    // TODO: to be implemented
    return 2;
  }
 private:
  static int Choose_N(int n, int q) {
    return ceil((2 * n + 1) * log((double)q));
  }

  static R_Ring_Number Noise(int q, int d, int B) {
    /*
    // TODO: to be implemented
    // int bound = floor(sqrt(q * 0.2 / Choose_N(n, q)));
    // assert(bound > 1);
    int bound = floor((q - 2) / 4.0 / d / (double)N);
    // LOG
    std::cout << "Noise bound = " << bound << std::endl;
    assert(bound >= 2);
    R_Ring_Number res = R_Ring_Number::Uniform_Rand(bound, d);
    res.Increase_Modul(q);
    return res;
    */
    assert(B >= 2);
    R_Ring_Number res = R_Ring_Number::Uniform_Rand(3, d);
    res.Increase_Modul(q);
    //    R_Ring_Number res(q, d);
    return res;
  }

  static int Choose_B(int q, int d, int N) {
    int B = floor((q - 2) / 4.0 / d / (double)N);
    if (B <= 1) {
      std::cout << "suggested q : q / log2(q) >= " << 8 * d * N / ceil(log2(1.0 * q)) + 2 << std::endl;
      std::cout << "q = " << q << std::endl;
      std::cout << "d = " << d << std::endl;
      std::cout << "N = " << N << std::endl;
      std::cout << "B = " << B << std::endl;
    }
    assert(B > 1);

    return B;
  }
public:	
  GLWE_Params Setup(int lambda, int mu, GLWE_Type b) const {
    int q = Choose_q(mu);
    int n = Choose_n(lambda, mu, b);
    int N = Choose_N(n, q);
    int d = Choose_d(lambda, mu, b);
    int B = Choose_B(q, d, N);

    return GLWE_Params(q, d, n, N, B, &Noise);
  }
	
  R_Ring_Vector Secret_Key_Gen(GLWE_Params &params) const {
    R_Ring_Vector sk(params.q, params.d, params.n + 1);
    sk[0] = 1;
    for (int i = 1; i < params.n + 1; i++) {
      sk[i] = params.ksi();
    }
    return sk; 
  }
  
  /*** Matrix enumeration goes as follows:
       ------------------>
       |
       |		A
       |
       |
       v
       
  ***/
  R_Ring_Matrix Public_Key_Gen(GLWE_Params &params, R_Ring_Vector &sk, R_Ring_Vector *ksi_noise_for_debug = NULL) const {
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
    /*std::cout << "minus_A_prime =";
    minus_A_prime.print();
    std::cout << std::endl;

    std::cout << "A_prime = ";
    A_prime.print();
    std::cout << std::endl;*/
    
    // s_prime
    assert(sk.Get_Dimension() == params.n + 1);
    R_Ring_Vector s_prime = sk.Get_Sub_Vector(1, params.n);
    assert(s_prime.Get_Dimension() == params.n);
    /*std::cout << "s_prime = ";
    s_prime.print();
    std::cout << std::endl;*/

    // *** debugging part ***

    // r
    R_Ring_Vector r(params.q, params.d, params.N);
    for (int i = 0; i < r.Get_Dimension(); i++) {
      r[i] = R_Ring_Number(A_prime.Get_q(), A_prime.Get_d(), R_Ring_Number::Uniform_Rand(A_prime.Get_q(), A_prime.Get_d()).Get_vec());
    }

    R_Ring_Matrix check1 = (A_prime * R_Ring_Matrix(s_prime)).Get_Transpose() * R_Ring_Matrix(r);
    R_Ring_Matrix check2 = R_Ring_Matrix(r).Get_Transpose() * (A_prime * R_Ring_Matrix(s_prime));
    R_Ring_Matrix check = check1 + (-check2);

    /*std::cout << "check = ";
    check.print();
    std::cout << std::endl;*/

    // end of debugging part
    
    // ksi_noise
    // R_Ring_Vector ksi_noise(params.q, params.d, params.N);
    R_Ring_Vector ksi_noise(params.q, params.d, params.N);
    for (int i = 0; i < params.N; i++) {
      ksi_noise[i] = params.ksi();
    }

    // pk
    R_Ring_Matrix pk(params.q, params.d, params.N, params.n + 1);
    R_Ring_Vector As = A_prime * s_prime;
    /*std::cout << "As = ";
    As.print();
    std::cout << std::endl;*/
    assert(As.Get_Dimension() == params.N);
    pk.Set_Column(0, As + ksi_noise * 2);
    pk.Set_Block(0, 1, -A_prime);

    if (ksi_noise_for_debug != NULL) {
      (*ksi_noise_for_debug) = ksi_noise;
    }
    return pk;
  }
  
  R_Ring_Vector Encrypt(GLWE_Params &params, R_Ring_Matrix &pk, R_Ring_Number &m, R_Ring_Vector *r_for_debug = NULL) const {
    // TODO: to implement for arbitrary m.q;
    assert(m.Get_q() == 2);
    
    // m_prime
    R_Ring_Vector m_prime(params.q, params.d, params.n + 1);
    m_prime[0] = R_Ring_Number(params.q, params.d, m.Get_vec());
    for (int i = 1; i < m_prime.Get_Dimension(); i++) {
      m_prime[i] = 0;
    }

    // r
    R_Ring_Vector r(params.q, params.d, params.N);
    for (int i = 0; i < r.Get_Dimension(); i++) {
      r[i] = R_Ring_Number(params.q, params.d, R_Ring_Number::Uniform_Rand(2, params.d).Get_vec());
    }

    /*std::cout << "r2 = ";
    r.print();
    std::cout << std::endl;*/
    
    R_Ring_Vector c(params.q, params.d, params.n + 1);
    R_Ring_Matrix pk_transpose = pk.Get_Transpose();
    /*std::cout << "pk_transpose = ";
    pk_transpose.print();
    std::cout << std::endl;*/
    
    c = m_prime + (pk.Get_Transpose() * r);
    /*std::cout << "pk.Get_Transpose() * r = ";
    (pk.Get_Transpose() * r).print();
    std::cout << std::endl;*/
    if (r_for_debug != NULL) {
      (*r_for_debug) = r;
    }
    return c;
  }
  
  static R_Ring_Number Decrypt(GLWE_Params &params, R_Ring_Vector &sk, R_Ring_Vector &c) {
    assert(c.Get_q() == params.q && sk.Get_q() == params.q);
    R_Ring_Number dot_product = c.Dot_Product(sk);
    /*std::cout << "dot_product = ";
    dot_product.print();
    std::cout << std::endl;*/
    //    R_Ring_Number clamped1 = dot_product.Clamp(params.q);
    /*std::cout << "clamped1 = ";
    clamped1.print();
    std::cout << std::endl;*/
    R_Ring_Number clamped2 = dot_product.Clamp(2);
    return clamped2;
  }
};

#endif /* _GLWE_H_ */

