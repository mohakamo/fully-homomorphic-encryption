#include "SI_HE.h"
#include "LSS.h"
#include <iostream>
#include <time.h>

enum SI_HE_Operation {SI_HE_Addition = 1, SI_HE_Multiplication};

/***
 * Tests section
 ***/
class Test {
private:
  void FAIL() {
    std::cout << "FAIL" << std::endl;
  }

  void PASS() {
    std::cout << "PASS" << std::endl;
  }

  bool test_Zero_Number() {
    std::cout << "testNumber_Zero ";
    R_Ring_Number zero(ZZ(INIT_VAL, 8), Ring_Number_d);
    for (int i = 0; i < Ring_Number_d; i++) {
      if (zero[i] != 0) {
	FAIL();
	return false;
      }
    }
    PASS();
    return true;
  }

  bool test_Nonzero_Number() {
    std::cout << "test_Nonzero_Number ";
    ZZ array_a[] = {ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 0)};
    R_Ring_Number a(ZZ(INIT_VAL, 2), Ring_Number_d, array_a);
    
    for (int i = 0; i < 4; i++) {
      if (a[i] != array_a[i]) {
	FAIL();
	return false;
      }
    }
    PASS();
    return true;
  }

  bool test_Number_Addition() {
    std::cout << "test_Number_Addition ";
    ZZ array_a[] = {ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 0)};
    ZZ array_b[] = {ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 1)};
    R_Ring_Number a(ZZ(INIT_VAL, 2), Ring_Number_d, array_a);
    R_Ring_Number b(ZZ(INIT_VAL, 2), Ring_Number_d, array_b);

    R_Ring_Number sum = a + b;
    sum.Clamp(ZZ(INIT_VAL, 2));
    if (sum[0] != 0 || sum[1] != 1 || sum[2] != 0 || sum[3] != 1) {
      FAIL();
      std::cout << "Expected: {0, 1, 0, 1}" << std::endl << "Got:      {" << sum[0] << ", " <<
	sum[1] << ", " << sum[2] << ", " << sum[3] << "}" << std::endl;
      return false;
    } else {
      PASS();
      return true;
    }
  }

  bool test_Number_Multiplication() {
    std::cout << "test_Number_Multiplication ";
    ZZ array_a[] = {ZZ(INIT_VAL, 2), ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 2), ZZ(INIT_VAL, 1)};
    ZZ array_b[] = {ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 1)};

    ZZ array_res[] = {ZZ(INIT_VAL, -1), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 0)};

    R_Ring_Number a(ZZ(INIT_VAL, 3), Ring_Number_d, array_a);
    R_Ring_Number b(ZZ(INIT_VAL, 3), Ring_Number_d, array_b);
    R_Ring_Number res(ZZ(INIT_VAL, 3), Ring_Number_d, array_res);

    R_Ring_Number mul = a * b;
    mul.Clamp(ZZ(INIT_VAL, 3));

    if (mul != res) {
      FAIL();
      std::cout << "Expected: {-1, 1, 1, 0}" << std::endl << "Got:      {" << mul[0] << ", " <<
	mul[1] << ", " << mul[2] << ", " << mul[3] << "}" << std::endl;
      std::cout << std::endl;
      return false;
    } else {
      PASS();
      return true;
    }
  }

  bool test_Number_Reduction() {
    std::cout << "test_Number_Reduction ";
    int tests[17][3] = {{2, 4, 2}, // 1
		      {3, 4, -1},
		      {6, 4, 2},
		      {-2, 4, 2},
		      {-6, 4, 2},
		      {2, 2, 0},
		      {1, 2, 1},
		      {0, 2, 0},
		      {-1, 2, 1},
		      {-63, 2, 1}, // 10
		      // odd moduls
		      {1, 3, 1},
		      {-1, 3, -1},
		      {0, 3, 0},
		      {-2, 3, 1},
		      {-4, 3, -1},
		      {2, 3, -1},
		      {3, 3, 0}}; // 17
    for (int i = 0; i < 17; i++) {
      if (R_Ring_Number::Clamp(ZZ(INIT_VAL, tests[i][0]), ZZ(INIT_VAL, tests[i][1])) != ZZ(INIT_VAL, tests[i][2])) {
	FAIL();
	std::cout << "R_Ring_Number::Reduce(" << tests[i][0] << ", " << tests[i][1] << ") != " << tests[i][2] << std::endl;
	std::cout << "R_Ring_Number::Reduce(...) = " << R_Ring_Number::Clamp(ZZ(INIT_VAL, tests[i][0]), ZZ(INIT_VAL, tests[i][1])) << std::endl;
	ZZ modul = ZZ(INIT_VAL, tests[i][1]);
	ZZ n = ZZ(INIT_VAL, tests[i][0]);

	// for odd module reduce to range [-(q - 1) / 2, (q - 1) / 2]
	// for even module reduce to range [(q - 2) / 2, q / 2]
	ZZ q_half_b = modul / 2;
	ZZ q_half_a = -(modul - 1) / 2;
	if (n > q_half_b) {
	  n -= ((n - q_half_b - 1) / modul + 1) * modul;
	} else if (n < q_half_a) {
	  n += (-(n - q_half_a + 1) / modul + 1) * modul;
	}
	std::cout << "q_half_b = " << q_half_b << ", q_half_a = " << q_half_a << std::endl;
	return false;
      }
    }
    PASS();
    return true;
  }

  bool test_Number_Equality() {
    std::cout << "test_Number_Equality ";

    ZZ array_a[] = {ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 0)};
    ZZ array_b[] = {ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 0)};
    ZZ array_c[] = {ZZ(INIT_VAL, 1), ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 0), ZZ(INIT_VAL, 0)};
    
    R_Ring_Number a(ZZ(INIT_VAL, 2), Ring_Number_d, array_a),
      b(ZZ(INIT_VAL, 2), Ring_Number_d, array_b),
      c(ZZ(INIT_VAL, 2), 4, array_c);
    if (!(a == b) || a != b || !(a != c) || a == c) {
      FAIL();
      return false;
    }
    PASS();
    return true;
  }

  bool test_Number_Rounding() {
    std::cout << "test_Number_Rounding ";
    ZZ q = to_ZZ(16);
    R_Ring_Number a(q + 1, 4);
    ZZ q_half = q / 2;
    a = q_half;

    R_Ring_Number one(to_ZZ(2), 4);
    one = to_ZZ(1);
    a = ((a * to_ZZ(2))) / q;

    if (a.Clamp(to_ZZ(2)) != one) {
      R_Ring_Number b(q + 1, 4);
      std::cout << "Final a = " << a << std::endl;
      b = q_half;
      std::cout << "a = q / 2 = " << b << std::endl;
      b = b + q_half;
      std::cout << "a = q / 2 + q / 2 = " << b << std::endl;
      b = (b * to_ZZ(2));
      std::cout << "a = 2 * q = " << b << std::endl;
      b = b / q;
      std::cout << "a = a / q = " << b << std::endl;
      b = b.Clamp(to_ZZ(2));
      std::cout << "a = Clamp(a, 2) = " << b << std::endl;
      
      FAIL();
      return false;
    }
    PASS();
    return true;
  }

  bool test_Regev_for_LWE(void) {
    std::cout << "test_Regev_for_LWE ";
    return test_Regev(LWE_Based);
  }

  bool test_Regev(GLWE_Type type) {
    Regev glwe;
    int moduls[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};

    for (int ii = 0; ii < 1; ii++) {
    for (int s = 0; s < 30; s++) {
      int modul = moduls[ii];
      Regev_Params params = glwe.Setup(type, ZZ(INIT_VAL, modul));
      assert(modul < params.q);
      Regev_Secret_Key_Type sk = glwe.Secret_Key_Gen(params);
      R_Ring_Vector ksi;
      Regev_Public_Key_Type pk = glwe.Public_Key_Gen(params, sk, &ksi);
      R_Ring_Number message = R_Ring_Number::Uniform_Rand(ZZ(INIT_VAL, modul), params.d);
      R_Ring_Vector r;
      Regev_Ciphertext_Type c = glwe.Encrypt(params, pk, message, &r);
      R_Ring_Number decoded_message = glwe.Decrypt(params, sk, c);

      if (message != decoded_message) {
	FAIL();
	std::cout << "Message modul = " << modul << std::endl;
	std::cout << "Cipher modul = " << params.q << std::endl;

	std::cout << "Initial message = " << message << "\n";

	ZZ q = ksi.Get_q();

	std::cout << "ksi = " << ksi << "\n";
	
	std::cout << "r = " << r << "\n";

	ksi.Increase_Modul(ksi.Get_q() * 4);
	r.Increase_Modul(r.Get_q() * 4);
	R_Ring_Vector noise_by_elements = ksi * r;
	std::cout << "noise_by_elements " << noise_by_elements << "\n";

	R_Ring_Number noise = ksi.Dot_Product(r);
	std::cout << "noise = " << noise << "\n";

	message.Increase_Modul(noise.Get_q());
	R_Ring_Number test_m = message + noise * ZZ(INIT_VAL, modul);
	std::cout << "test_m = " << test_m << "\n";

	R_Ring_Number th_result = test_m.Clamp(q);
	std::cout << "th_result = " << th_result << "\n";

	th_result = th_result.Clamp(ZZ(INIT_VAL, modul));
	std::cout << "th_result = " << th_result << "\n";

	std::cout << "On attempt #" << s << std::endl;

	std::cout << "Secret key: " << sk << "\n";

	std::cout << "Ciphertext: " << c << "\n";

	std::cout << "Result of coding/deconging message " << message << " is " << decoded_message << "\n";
	return false;
      }
    }
    }
    PASS();
    return true;
  }
 
  bool test_Regev_for_RLWE(void) {
    std::cout << "test_Regev_for_RLWE ";
    return test_Regev(RLWE_Based);
  }

  bool test_SI_HE(GLWE_Type type) {
    SI_HE fhe;
    
    int modules[] = {2};//97, 101, 103, 107, 109};
    Regev_Params params = fhe.Setup(2, type, ZZ(INIT_VAL, modules[0]));
    for (int mi = 0; mi < 1; mi++) {
      ZZ modul = ZZ(INIT_VAL, modules[mi]);
    for (int s = 0; s < 10; s++) {
      //      std::cout << s << std::endl;
      SI_HE_Secret_Key_Type sk = fhe.Secret_Key_Gen(params);
      SI_HE_Public_Key_Type pk = fhe.Public_Key_Gen(params, sk);
      SI_HE_Evaluation_Key_Type evalk = fhe.Evaluation_Key_Gen(params, sk, pk);
      
      R_Ring_Number message = R_Ring_Number::Uniform_Rand(modul, params.d);
      
      SI_HE_Cipher_Text c = fhe.Encrypt(params, &pk, message, &evalk);
      R_Ring_Number decoded_message = c.Decrypt(params, sk);

      
      if (message != decoded_message) {
	FAIL();
	return false;
      }    
    }
    }
    PASS();
    return true;
  }

  bool test_SI_HE_for_RLWE(void) {
    std::cout << "test_SI_HE_RLWE ";
    return test_SI_HE(RLWE_Based);
  }
  
  bool test_SI_HE_for_LWE(void) {
    std::cout << "test_SI_HE_LWE ";
    return test_SI_HE(LWE_Based);
  }

  bool test_Powersof2_BitDecomposition(void) {
    std::cout << "test_Powersof2_BitDecomposition ";
    ZZ q = ZZ(INIT_VAL, 7);
    int d = 3;
    for (int i = 0; i < 30; i++) {
      R_Ring_Vector c = R_Ring_Vector::Uniform_Rand(q, d, 10);
      R_Ring_Vector s = R_Ring_Vector::Uniform_Rand(q, d, 10);
    R_Ring_Number expected = c.Dot_Product(s).Clamp(q);
    R_Ring_Vector bits_c;
    Bit_Decomposition(c, q, bits_c);
    R_Ring_Vector powers2_s;
    Powersof2(s, q, powers2_s);
    R_Ring_Number actual = bits_c.Dot_Product(powers2_s).Clamp(q);

    if (expected != actual) {
      std::cout << "c = " << c << "\n";
      std::cout << "s = " << s << "\n";
      std::cout << "BitDecomp(c, q) = " << bits_c << "\n";
      std::cout << "Powersof2(s, q) = " << powers2_s << "\n";
      std::cout << "Expected result: <c, s> = " << expected << "\n";
      std::cout << "Actual result: <BitDecomp(c, q), Powerspf2(s, q)> = " << actual << "\n";
	
      FAIL();
      std::cout << std::endl;
      return false;
    }
    }
    PASS();
    return true;
  }

  bool test_Number_Scale(void) {
    std::cout << "test_Number_Scale ";
    int d = 1;
    int dimension = 10;
    ZZ q[] = {ZZ(INIT_VAL, 577), ZZ(INIT_VAL, 613), ZZ(INIT_VAL, 983)};
    ZZ p[] = {ZZ(INIT_VAL, 571), ZZ(INIT_VAL, 607), ZZ(INIT_VAL, 977)};
    for (int i = 0; i < 3; i++) {
      for (int s = 0; s < 30; s++) {
	ZZ bound = q[i] / 2 / dimension - q[i] / p[i] * d;
	R_Ring_Vector message = R_Ring_Vector::Uniform_Rand(ZZ(INIT_VAL, 2), d, dimension);
	R_Ring_Vector qv = R_Ring_Vector::Uniform_Rand(q[i], d, dimension, bound);
	R_Ring_Vector pv = R_Ring_Vector(p[i], d, dimension);
	for (int j = 0; j < dimension; j++) {
	  pv[j] = qv[j].Scale(q[i], p[i], ZZ(INIT_VAL, 2));
	}
	message.Increase_Modul(p[i]);
	R_Ring_Number result_p = message.Dot_Product(pv);
	message.Increase_Modul(q[i]);
	R_Ring_Number result_q = message.Dot_Product(qv);
	if (result_q.Get_Clamped(ZZ(INIT_VAL, 2)) != result_p.Get_Clamped(ZZ(INIT_VAL, 2))) {
	  std::cout << "attempt #" << s << std::endl;
	  std::cout << "q = " << q[i] << " p = " << p[i] << std::endl;
	  std::cout << "bound = " << bound << std::endl;
	  std::cout << "qv = " << qv << "\n";

	  qv.Clamp(ZZ(INIT_VAL, 2));
	  std::cout << "qv mod 2 = " << qv << "\n";
	  std::cout << "pv = " << pv << "\n";

	  pv.Clamp(ZZ(INIT_VAL, 2));
	  std::cout << "pv mod 2 = " << pv << "\n";
	  std::cout << "message = " << message << "\n";
	  std::cout << "<qv, message> = " << result_q << "\n";
	  std::cout << "<pv, message> = " << result_p << "\n";

	  FAIL();
	  return false;
	}
      }
    }
    /*
    R_Ring_Number aa_n(ZZ(INIT_VAL, 4194301), 1), bb_n(ZZ(INIT_VAL, 262135), 1);
    aa_n[0] = ZZ(INIT_VAL, -1386462);
    bb_n[0] = ZZ(INIT_VAL, -86652);
    R_Ring_Number res = aa_n.Scale(ZZ(INIT_VAL, 4194301), ZZ(INIT_VAL, 262135), ZZ(INIT_VAL, 3));
    assert(res == bb_n); */
    PASS();
    return true;
  }

  bool test_SI_HE_Switch_Keys() {
    test_SI_HE_Switch_Keys(LWE_Based);
    test_SI_HE_Switch_Keys(RLWE_Based);
  }

  bool test_SI_HE_Switch_Keys(GLWE_Type type) {
    std::cout << "test_SI_HE_Switch_Keys_" << ((type == LWE_Based) ? "LWE " : "RLWE ");
    SI_HE fhe;
    int L = 10;
    int modul = 2;

    for (int s = 0; s < 10; s++) {
      Regev_Params params = fhe.Setup(L, type, ZZ(INIT_VAL, modul));
      SI_HE_Secret_Key_Type sk = fhe.Secret_Key_Gen(params);
      SI_HE_Public_Key_Type pk = fhe.Public_Key_Gen(params, sk);
      SI_HE_Evaluation_Key_Type evalk = fhe.Evaluation_Key_Gen(params, sk, pk);
      
      R_Ring_Number message;
      SI_HE_Cipher_Text c;
      message = R_Ring_Number::Uniform_Rand(to_ZZ(2), params.d);
      c = fhe.Encrypt(params, &pk, message, &evalk);
      c.Add_Secret_Key_Info(&sk);
      assert(c.Decrypt(params, sk) == message);
      
      for (int i = 0; i < L; i++) {
	c.Raise_Level();
	if (c.Decrypt(params, sk) != message) {
	  std::cout << "attempt #" << s << "\n";
	  FAIL();
	  return false;
	}
      }
    }

    PASS();
    return true;
  }

  bool test_Multiple_Refresh(GLWE_Type type) {
    return true;
  }

  bool test_SI_HE_Operations(GLWE_Type type, SI_HE_Operation operation_type, int noof_mults) {
    assert(noof_mults == 1 || noof_mults == 2);
    std::cout << "test_" << noof_mults << "_" << ((operation_type == SI_HE_Addition) ? "Add_" : "Mult_") << ((type == LWE_Based) ? "LWE " : "RLWE ");
    SI_HE fhe;
    int modules[] = {2, 97, 7, 11, 3, 5};
    int L = noof_mults;

    for (int t = 0; t < 1; t++) {
      int modul = modules[t];
      for (int s = 0; s < 10; s++) {
	Regev_Params params = fhe.Setup(L, type, ZZ(INIT_VAL, modul));
	SI_HE_Secret_Key_Type sk = fhe.Secret_Key_Gen(params);
	SI_HE_Public_Key_Type pk = fhe.Public_Key_Gen(params, sk);
	SI_HE_Evaluation_Key_Type evalk = fhe.Evaluation_Key_Gen(params, sk, pk);

	R_Ring_Number message[3];
	SI_HE_Cipher_Text c[3];
	for (int j = 0; j < 3; j++) {
	  message[j] = R_Ring_Number::Uniform_Rand(to_ZZ(2), params.d);
	  c[j] = fhe.Encrypt(params, &pk, message[j], &evalk);
	  c[j].Add_Secret_Key_Info(&sk);
	  assert(c[j].Decrypt(params, sk) == message[j]);
	}
	
	// first level of addition
	SI_HE_Cipher_Text res_c = (operation_type == SI_HE_Addition ? c[0] + c[1] : c[0] * c[1]);
	R_Ring_Number res_m = (operation_type == SI_HE_Addition ? message[0] + message[1] : message[0] * message[1]);
	R_Ring_Number res_m_decr = res_c.Decrypt(params, sk);

	if (res_m_decr != res_m) {
	  std::cout << "One " << (operation_type == SI_HE_Addition ? "addition" : "multiplication") << " failed to decrypt" << std::endl;
	  std::cout << "attempt #" << s << "\n";
	  std::cout << "message[0] = " << message[0] << "\n";
	  std::cout << "message[1] = " << message[1] << "\n";
	  std::cout << "(message[0] o message[1]).Decrypt() = " << res_m_decr << "\n";
	  std::cout << "(message[0] o message[1]) = " << res_m << "\n";
	  std::cout << "c[0] = " << c[0] << "\n";
	  std::cout << "c[1] = " << c[1] << "\n";
	  std::cout << "res_c = " << res_c << "\n";
	  FAIL();
	  return false;
	}
	if (noof_mults == 1) continue;

	SI_HE_Cipher_Text cipher2_copy = c[2];

	SI_HE_Cipher_Text res_c2 = (operation_type == SI_HE_Addition ? res_c + cipher2_copy : res_c * cipher2_copy);
	R_Ring_Number res_m2 = (operation_type == SI_HE_Addition ? res_m + message[2] : res_m * message[2]);
	R_Ring_Number res_m_decr2 = res_c2.Decrypt(params, sk);
	if (res_m_decr2 != res_m2) {
	  std::cout << "Fail on second " << (operation_type == SI_HE_Addition ? "add" : "mult") << std::endl;
	  std::cout << "Attempt #" << s << std::endl;

	  std::cout << "(m1 o m2) o m3 = (" << message[0] << " * " << message[1] << ") * " << message[2] << std::endl;
	  std::cout << "(c1 o c2) o c3 = (" << c[0] << " * " << c[1] << ") * " << c[2] << std::endl;
	  std::cout << "(c1 o c2) o c3 = " << res_c << " * " << cipher2_copy << std::endl;
	  std::cout << "(c1 o c2) o c3 = " << res_c2 << std::endl;
	  std::cout << "((c1 o c2) o c3).Decrypt() = " << res_m_decr2 << std::endl;
	  std::cout << "(c1 o c2) = " << res_c << std::endl;
	  std::cout << "(c1 o c2).Decrypt() = " << res_m << std::endl;
	  std::cout << "c3 = " << c[2] << " => " << cipher2_copy << std::endl;
	  
	  FAIL();
	  return false;
	}
    }}
    PASS();
    return true;
  }    

  bool test_SI_HE_Operations(void) {
      return
	test_SI_HE_Operations(LWE_Based, SI_HE_Addition, 1) &&
	test_SI_HE_Operations(LWE_Based, SI_HE_Addition, 2) &&
	test_SI_HE_Operations(LWE_Based, SI_HE_Multiplication, 1) &&
	test_SI_HE_Operations(LWE_Based, SI_HE_Multiplication, 2) &&

	test_SI_HE_Operations(RLWE_Based, SI_HE_Addition, 1) &&
	test_SI_HE_Operations(RLWE_Based, SI_HE_Addition, 2) &&
	test_SI_HE_Operations(RLWE_Based, SI_HE_Multiplication, 1) &&
	test_SI_HE_Operations(RLWE_Based, SI_HE_Multiplication, 2);
  }

  bool test_SI_HE_LSS() {
    std::cout << "test_SI_HE_LSS ";
    
    GLWE_Type type = LWE_Based; // RLWE_Based
    int message_modul = 2;
    ZZ modul;
    modul = message_modul;
    SI_HE fhe;
    int L = 2;
    Regev_Params params = fhe.Setup(L, type, modul);
    std::cout << "params = ";
    params.print();
    std::cout << std::endl;
    
    int noof_vectors = 5;
    SI_HE_Secret_Key_Type sk = fhe.Secret_Key_Gen(params);
    SI_HE_Public_Key_Type pk = fhe.Public_Key_Gen(params, sk);
    SI_HE_Evaluation_Key_Type evalk = fhe.Evaluation_Key_Gen(params, sk, pk);

    std::vector<R_Ring_Number> messages_x, messages_y;
    std::vector<SI_HE_Cipher_Text> c_x, c_y;
    for (int j = 0; j < noof_vectors; j++) {
      R_Ring_Number mx = R_Ring_Number::Uniform_Rand(modul, params.d);
      R_Ring_Number my = R_Ring_Number::Uniform_Rand(modul, params.d);

      SI_HE_Cipher_Text cx = fhe.Encrypt(params, &pk, mx, &evalk);
      SI_HE_Cipher_Text cy = fhe.Encrypt(params, &pk, my, &evalk);
      messages_x.push_back(mx);
      messages_y.push_back(my);
      cx.Add_Secret_Key_Info(&sk);
      cy.Add_Secret_Key_Info(&sk);

      c_x.push_back(cx);
      c_y.push_back(cy);
    }
      
    Pair<SI_HE_Cipher_Text, SI_HE_Cipher_Text> res_c = Compute_LSS<SI_HE_Cipher_Text>(c_x, c_y);
    Pair<R_Ring_Number, R_Ring_Number> res = Compute_LSS<R_Ring_Number>(messages_x, messages_y);
    Pair<R_Ring_Number, R_Ring_Number> res_d(res_c.first.Decrypt(params, sk),
					     res_c.second.Decrypt(params, sk));
    if (res.first != res_d.first || res.second != res_d.second) {
      std::cout << "real result = (" << res.first << ", " << res.second << ")" << std::endl;
      std::cout << "fhe result = (" << res_d.first << ", " << res_d.second << ")" << std::endl;

      std::cout << "Initial vectors:" << std::endl;
      for (int i = 0; i < messages_x.size(); i++) {
	std::cout << "m" << i << " = (" << messages_x[i] << ", " << messages_y[i] << ")\n";
      }
      FAIL();
      return false;
    }
    PASS();
    return true;
  }

  bool test_LSS(void) {
    std::cout << "test_LSS ";

    for (int d = 1; d < 4; d++) {
      for (int s = 0; s < 30; s++) {
    
    int noof_vectors = 2;
    int message_modul_bound = 10;
    ZZ modul;
    modul = 2 * noof_vectors * noof_vectors * message_modul_bound * message_modul_bound * message_modul_bound  + 1;
    while (ProbPrime(modul)) {
      modul++;
    }

    std::vector<ZZ *> array_m_x(noof_vectors), array_m_y(noof_vectors);
    for (int j = 0; j < noof_vectors; j++) {
      array_m_x[j] = new ZZ [d];
      array_m_x[j][0] = R_Ring_Number::Clamp(ZZ(INIT_VAL, rand()), ZZ(INIT_VAL, message_modul_bound));
      array_m_y[j] = new ZZ [d];
      array_m_y[j][0] = R_Ring_Number::Clamp(ZZ(INIT_VAL, rand()), ZZ(INIT_VAL, message_modul_bound));
    }
    std::vector<R_Ring_Number> messages_x, messages_y;
    for (int j = 0; j < noof_vectors; j++) {
      messages_x.push_back(R_Ring_Number(modul, d, array_m_x[j]));
      messages_y.push_back(R_Ring_Number(modul, d, array_m_y[j]));
    }

    for (int j = 0; j < noof_vectors; j++) {
      delete [] array_m_x[j];
      delete [] array_m_y[j];
    }

    R_Ring_Number den;
      
    Pair<R_Ring_Number, R_Ring_Number> res = Compute_LSS<R_Ring_Number>(messages_x, messages_y, &den);
    if (res.first * messages_x[0] + res.second != messages_y[0] * den ||
	res.first * messages_x[1] + res.second != messages_y[1] * den) {
      std::cout << "d = " << d << ", s = " << s << ", message modul = " << modul << std::endl;
      std::cout << "(x0, y0) = (" << messages_x[0] << ", " << messages_y[0] << ")" << std::endl;
      std::cout << "(x1, y1) = (" << messages_x[1] << ", " << messages_y[1] << ")" << std::endl;
      std::cout << "(a, b) = (" << res.first << ", " << res.second << ")" << std::endl;
      FAIL();
      return false;
    }}}
    PASS();
    return true;
  }

public:
  Test() {
    Ring_Number_d = 4;
  }
  void Run_Tests() {
    if (!test_Zero_Number() ||
	!test_Nonzero_Number() ||
	!test_Number_Addition() ||
	!test_Number_Multiplication() ||
	!test_Number_Equality() ||
	!test_Number_Reduction() ||
	!test_Number_Rounding() ||

	!test_Regev_for_LWE() ||
	!test_Regev_for_RLWE() ||

	!test_SI_HE_for_LWE() ||
	!test_SI_HE_for_RLWE() ||
	!test_Powersof2_BitDecomposition() ||
	!test_Number_Scale() ||
	!test_SI_HE_Switch_Keys() ||

	!test_Multiple_Refresh(LWE_Based) ||
	!test_Multiple_Refresh(RLWE_Based) ||
	!test_SI_HE_Operations() ||
	!test_SI_HE_LSS() ||
	!test_LSS()) {
      std::cout << "Overall tests FAILED" << std::endl;
    } else {
      std::cout << "Overall tests PASSED" << std::endl;
    }
  }
};
