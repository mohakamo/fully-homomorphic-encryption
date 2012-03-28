/**
 *  Regev.h
 *  Created by valerini on 1/5/12.
 */
#ifndef _REGEV_H_
#define _REGEV_H_
#include "NormDistr.h"
#include "R_Ring_Matrix.h"
#include <math.h>
#include <map>
// Getting GLWE_Type enum from GLWE.h
#include "GLWE.h"

/*** Params for Regev Encrption Scheme ***/
class Regev_Params {
public:
  Regev_Params() {
    deviation = d = n = N = 0;
    p = q = ZZ::zero();
    noise = NULL;
  }
  Regev_Params(ZZ q, int d, int n, int N, int deviation, ZZ p,
	       R_Ring_Number (*noise)(ZZ q, int d, int deviation)) {
    this->q = q;
    this->d = d;
    this->n = n;
    this->N = N;
    this->deviation = deviation;
    this->p = p;
    this->noise = noise;
  }
  ZZ q, p;
  int d, n, N, deviation;
  /*** Noise distribution ***/
  R_Ring_Number (*noise)(ZZ q, int d, int deviation);
  
  R_Ring_Number ksi() const {
    return noise(q, d, deviation);
  }
  
  Regev_Params(const Regev_Params &r) {
    q = r.q;
    d = r.d;
    n = r.n;
    N = r.N;
    deviation = r.deviation;
    p = r.p;
    noise = r.noise;
  }

  void print() {
    std::cout << "(q=" << q << ", ";
    std::cout << "d=" << d << ", ";
    std::cout << "n=" << n << ", ";
    std::cout << "N=" << N << ", ";
    std::cout << "deviation=" << deviation << ",";
    std::cout << "p=" << p << ")";
  }
};

typedef R_Ring_Vector Regev_Secret_Key_Type;
typedef R_Ring_Matrix Regev_Public_Key_Type;
typedef R_Ring_Vector Regev_Ciphertext_Type;

/*** Basic Regev_Encryption_Scheme ***/
class Regev {
 public:
  // q is chosen to be prime of length mu such that q = 1 (mod r), where r is the message module
  static ZZ Choose_q(int n) {
    // TODO:
    return to_ZZ(1) << 20;
  }
 public:
  static int Choose_d(GLWE_Type b) {
    if (b == LWE_Based) {
      return 1;
    }
    // TODO: to be implemented, seems like a reasonable value for noise >= 2^10
    return 8;
  }

  static int Choose_n(GLWE_Type b) {
    if (b == RLWE_Based) {
      return 1;
    }
    return 8;
  }
 private:
  static int Choose_N(int n, ZZ q) {
    // Should be N := (n + 1) * (NumBits(q) + O(1));
    return ceil((n + 1) * NumBits(q));
  }

  // just return the uniform noise that is bounded by B (although should be gaussian - not uniform for security)
  static R_Ring_Number Noise(ZZ q, int d, int deviation = 8) {
    // R_Ring_Number res(q, d);
    //    R_Ring_Number res = NormDistr::sample(q, d, deviation);
    R_Ring_Number res = R_Ring_Number::Uniform_Rand(to_ZZ(2), d);
    res.Increase_Modul(q);
    return res;
  }

public:	
  Regev_Params Setup(GLWE_Type b, ZZ p = ZZ(INIT_VAL, 2), int n_opt = -1) const {
    int n;
    assert(p == to_ZZ(2));
    if (n_opt != -1) {
      n = n_opt;
    } else {
      n = Choose_n(b);
    }
    ZZ q = Choose_q(n);
    int N = Choose_N(n, q);
    int d = Choose_d(b);
    // NB! Setting global field
    Ring_Number_d = d;
    int deviation = 8;

    if ((n + 1) * NumBits(q) * 2 >= q / 4) {
      std::cout << "(n + 1) * NumBits(q) * 2 = " << (n + 1) * NumBits(q) * 2 << std::endl;
      std::cout << "q / 4 = " << q / 4 << std::endl;
      //      assert((n + 1) * NumBits(q) * 2 < q / 4);
    }

    return Regev_Params(q, d, n, N, deviation, p, &Noise);
  }
	
  R_Ring_Vector Secret_Key_Gen(Regev_Params &params) const {
    R_Ring_Vector sk(params.q, params.d, params.n);
    for (int i = 0; i < params.n; i++) {
      sk[i] = R_Ring_Number::Uniform_Rand(params.q, params.d);//params.ksi();
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
  R_Ring_Matrix Public_Key_Gen(const Regev_Params &params,
			       const R_Ring_Vector &sk,
			       R_Ring_Vector *ksi_noise_for_debug = NULL) const {
    R_Ring_Matrix A_prime = R_Ring_Matrix::Uniform_Rand(params.q, params.d, params.N, params.n);
    R_Ring_Vector ksi_noise(params.q, params.d, params.N);
    for (int i = 0; i < params.N; i++) {
      ksi_noise[i] = params.ksi();
    }

    R_Ring_Matrix pk(params.q, params.d, params.N, params.n + 1);
    pk.Set_Column(0, A_prime * sk + ksi_noise);
    pk.Set_Block(0, 1, -A_prime);

    if (ksi_noise_for_debug != NULL) {
      (*ksi_noise_for_debug) = ksi_noise;
    }
    return pk;
  }
  
  R_Ring_Vector Encrypt(Regev_Params &params, R_Ring_Matrix &pk, R_Ring_Number &m,
			R_Ring_Vector *r_for_debug = NULL) const {
    // TODO: to change back
    R_Ring_Vector r = R_Ring_Vector::Uniform_Rand(params.q, params.d, params.N, ZZ(INIT_VAL, 2));
    // R_Ring_Vector r(params.q, params.d, params.N);

    R_Ring_Vector m_prime(params.q, params.d, params.n + 1);
    m_prime[0] = R_Ring_Number(params.q, params.d, m.vec) * (params.q / 2);
    //    std::cout << "m_prime = " << m_prime << std::endl;

    R_Ring_Vector c(params.q, params.d, params.n + 1);
    c = m_prime + (pk.Get_Transpose() * r);

    if (r_for_debug != NULL) {
      (*r_for_debug) = r;
    }
    return c;
  }

  static R_Ring_Number Decrypt(const Regev_Params &params, const R_Ring_Vector &sk, const R_Ring_Vector &c) {
    R_Ring_Vector sk_prime(sk.Get_q(), sk.Get_d(), sk.Get_Dimension() + 1);
    sk_prime[0] = 1;
    for (int i = 1; i <= sk.Get_Dimension(); i++) {
      sk_prime[i] = sk[i - 1];
    }
    // Assume that division in NTL rounds down
    R_Ring_Number dot_product = (c.Dot_Product(sk_prime));
    R_Ring_Number rounding_add_constant(dot_product.Get_q(), dot_product.Get_d());
    ZZ q_quater = params.q / 4;
    ZZ q_half = params.q / 2;
    for (int i = 0; i < dot_product.Get_d(); i++) {
      rounding_add_constant[i] = q_quater;
    }
    //        std::cout << "dot_product = " << dot_product << std::endl;
    dot_product = (dot_product + rounding_add_constant) / q_half;
    //        std::cout << "dot_product / q_half = " << dot_product << std::endl;

    return dot_product.Clamp(params.p);
  }

  static R_Ring_Number Get_Noise(Regev_Params &params, R_Ring_Vector &sk, const R_Ring_Vector &c) {
    return c.Dot_Product(sk);
  }
};

#endif /* _REGEV_H_ */

