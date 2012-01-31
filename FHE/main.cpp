#include "FHE.h"
#include <iostream>

enum FHE_Operation {FHE_Addition = 1, FHE_Multiplication};

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
    R_Ring_Number zero(8, 8);
    for (int i = 0; i < 8; i++) {
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
    long long array_a[] = {0, 1, 1, 0};
    R_Ring_Number a(2, 4, array_a);
    
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
    long long array_a[] = {0, 1, 1, 0};
    long long array_b[] = {0, 0, 1, 1};
    R_Ring_Number a(2, 4, array_a);
    R_Ring_Number b(2, 4, array_b);

    R_Ring_Number sum = a + b;
    sum.Clamp(2);
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
    long long array_a[] = {2, 0, 2, 1};
    long long array_b[] = {1, 0, 0, 1};

    R_Ring_Number a(3, 4, array_a);
    R_Ring_Number b(3, 4, array_b);

    R_Ring_Number mul = a * b;
    mul.Clamp(3);

    if (mul[0] != -1 || mul[1] != 1 || mul[2] != 1 || mul[3] != 0) {
      FAIL();
      std::cout << "Expected: {-1, 1, 1, 0}" << std::endl << "Got:      {" << mul[0] << ", " <<
	mul[1] << ", " << mul[2] << ", " << mul[3] << "}" << std::endl;
      std::cout << std::endl;
      return false;
    } else {
      PASS();
      std::cout << std::endl;
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
      if (R_Ring_Number::Clamp(tests[i][0], tests[i][1]) != tests[i][2]) {
	FAIL();
	std::cout << "R_Ring_Number::Reduce(" << tests[i][0] << ", " << tests[i][1] << ") != " << tests[i][2] << std::endl;
	return false;
      }
    }
    PASS();
    return true;
  }

  bool test_Number_Equality() {
    std::cout << "test_Number_Equality ";

    long long array_a[] = {1, 1, 0, 0};
    long long array_b[] = {1, 1, 0, 0};
    long long array_c[] = {1, 0, 0, 0};
    
    R_Ring_Number a(2, 4, array_a),
      b(2, 4, array_b),
      c(2, 4, array_c);
    if (!(a == b) || a != b || !(a != c) || a == c) {
      FAIL();
      return false;
    }
    PASS();
    return true;
  }

  bool test_GLWE_for_LWE(void) {
    std::cout << "test_GLWE_for_LWE ";
    return test_GLWE(LWE_Based);
  }

  bool test_GLWE(GLWE_Type type) {
    GLWE glwe;
    int moduls[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};

    for (int ii = 0; ii < 9; ii++) {
    for (int s = 0; s < 30; s++) {
      int modul = moduls[ii];
      GLWE_Params params = glwe.Setup(2, 13, type, modul);
      assert(modul < params.q);
      GLWE_Secret_Key_Type sk = glwe.Secret_Key_Gen(params);
      R_Ring_Vector ksi;
      GLWE_Public_Key_Type pk = glwe.Public_Key_Gen(params, sk, &ksi);
      
      long long *array_m = new long long [params.d];
      for (int j = 0; j < params.d; j++) {
	array_m[j] = R_Ring_Number::Clamp(rand(), modul);
	// array_m[j] = rand() % modul;
      }
      R_Ring_Number message(modul, params.d, array_m);
      R_Ring_Vector r;
      GLWE_Ciphertext_Type c = glwe.Encrypt(params, pk, message, &r);
      R_Ring_Number decoded_message = glwe.Decrypt(params, sk, c);

      delete [] array_m;

      if (message != decoded_message) {
	FAIL();
	std::cout << "Message modul = " << modul << std::endl;
	std::cout << "Cipher modul = " << params.q << std::endl;

	std::cout << "Initial message = ";
	message.print();
	std::cout << std::endl;

	long long q = ksi.Get_q();

	std::cout << "ksi = ";
	ksi.print();
	std::cout << std::endl;
	
	std::cout << "r = ";
	r.print();
	std::cout << std::endl;

	ksi.Increase_Modul(ksi.Get_q() * 4);
	r.Increase_Modul(r.Get_q() * 4);
	R_Ring_Vector noise_by_elements = ksi * r;
	std::cout << "noise_by_elements";
	noise_by_elements.print();
	std::cout << std::endl;
	R_Ring_Number noise = ksi.Dot_Product(r);
	std::cout << "noise = ";
	noise.print();
	std::cout << std::endl;
	message.Increase_Modul(noise.Get_q());
	R_Ring_Number test_m = message + noise * modul;
	std::cout << "test_m = ";
	test_m.print();
	std::cout << std::endl;

	R_Ring_Number th_result = test_m.Clamp(q);
	std::cout << "th_result = ";
	th_result.print();
	std::cout << std::endl;

	th_result = th_result.Clamp(modul);
	std::cout << "th_result = ";
	th_result.print();
	std::cout << std::endl;

	std::cout << "On attempt #" << s << std::endl;

	std::cout << "Secret key: ";
	sk.print();
	std::cout << std::endl;

	std::cout << "Public key: ";
	pk.print();
	std::cout << std::endl;

	std::cout << "Ciphertext: ";
	c.print();
	std::cout << std::endl;

	std::cout << "Result of coding/deconging message ";
	message.print();
	std::cout << " is ";
	decoded_message.print();
	std::cout << "\n";
	return false;
      }
    }
    }
    PASS();
    return true;
  }
 
  bool test_GLWE_for_RLWE(void) {
    std::cout << "test_GLWE_for_RLWE ";
    return test_GLWE(RLWE_Based);
  }

  bool test_FHE(GLWE_Type type) {
    FHE fhe;
    
    for (int s = 0; s < 30; s++) {
      FHE_Params params = fhe.Setup(2, 2, type);
      Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> sk_pk = fhe.Key_Gen(params);
      
      long long *array_m = new long long [params[0].d];
      for (int i = 0; i < params[0].d; i++) {
	array_m[i] = rand() % 2;
      }
      R_Ring_Number message(2, params[0].d, array_m);
      
      FHE_Cipher_Text c = fhe.Encrypt(params, &sk_pk.second, message);
      R_Ring_Number decoded_message = c.Decrypt(params, sk_pk.first);

      delete [] array_m;
      
      if (message != decoded_message) {
	std::cout << "Secret key: ";
	for (std::vector<R_Ring_Vector>::iterator i = sk_pk.first.begin(); i != sk_pk.first.end(); i++) {
	  i->print();
	}
	std::cout << std::endl;
	
	std::cout << "Public key: ";
	for (std::vector<R_Ring_Matrix>::iterator i = sk_pk.second.begin(); i != sk_pk.second.end(); i++) {
	  i->print();
	}
	std::cout << std::endl;
	
	std::cout << "Ciphertext: ";
	c.print();
	std::cout << std::endl;

	std::cout << "Result of coding/deconging message ";
	message.print();
	std::cout << " is ";
	decoded_message.print();
	std::cout << "\n";
	FAIL();
	return false;
      }    
    }
    PASS();
    return true;
  }

  bool test_FHE_for_RLWE(void) {
    std::cout << "test_FHE_RLWE ";
    return test_FHE(RLWE_Based);
  }
  
  bool test_FHE_for_LWE(void) {
    std::cout << "test_FHE_LWE ";
    return test_FHE(LWE_Based);
  }

  bool test_Powersof2_BitDecomposition(void) {
    std::cout << "test_Powersof2_BitDecomposition ";
    int q = 7;
    int d = 3;
    for (int i = 0; i < 30; i++) {
    R_Ring_Vector c = R_Ring_Vector::Uniform_Rand(q, d, 10);
    R_Ring_Vector s = R_Ring_Vector::Uniform_Rand(q, d, 10);
    R_Ring_Number expected = c.Dot_Product(s).Clamp(q);
    R_Ring_Vector bits_c = FHE::Bit_Decomposition(c, q);
    R_Ring_Vector powers2_s = FHE::Powersof2(s, q);
    R_Ring_Number actual = bits_c.Dot_Product(powers2_s).Clamp(q);

    if (expected != actual) {
      std::cout << "c = ";
      c.print();
      std::cout << std::endl;

      std::cout << "s = ";
      s.print();
      std::cout << std::endl;

      std::cout << "BitDecomp(c, q) = ";
      bits_c.print();
      std::cout << std::endl;

      std::cout << "Powersof2(s, q) = ";
      powers2_s.print();
      std::cout << std::endl;

      std::cout << "Expected result: <c, s> = ";
      expected.print();
      std::cout << std::endl;

      std::cout << "Actual result: <BitDecomp(c, q), Powerspf2(s, q)> = ";
      actual.print();
      std::cout << std::endl;
	
      FAIL();
      std::cout << std::endl;
      return false;
    }
    }
    PASS();
    std::cout << std::endl;
    return true;
  }

  bool test_Number_Scale(void) {
    std::cout << "test_Number_Scale ";
    int d = 1;
    int dimension = 10;
    int q[] = {577, 613, 983};
    int p[] = {571, 607, 977};
    for (int i = 0; i < 3; i++) {
      for (int s = 0; s < 30; s++) {
	int bound = q[i] / 2 / dimension - q[i] / p[i] * d;
	R_Ring_Vector message = R_Ring_Vector::Uniform_Rand(2, d, dimension);
	R_Ring_Vector qv = R_Ring_Vector::Uniform_Rand(q[i], d, dimension, bound);
	R_Ring_Vector pv = R_Ring_Vector(p[i], d, dimension);
	for (int j = 0; j < dimension; j++) {
	  pv[j] = qv[j].Scale(q[i], p[i], 2);
	}
	message.Increase_Modul(p[i]);
	R_Ring_Number result_p = message.Dot_Product(pv);
	message.Increase_Modul(q[i]);
	R_Ring_Number result_q = message.Dot_Product(qv);
	if (result_q.Get_Clamped(2) != result_p.Get_Clamped(2)) {
	  std::cout << "attempt #" << s << std::endl;
	  std::cout << "q = " << q[i] << " p = " << p[i] << std::endl;
	  std::cout << "bound = " << bound << std::endl;
	  
	  std::cout << "qv = ";
	  qv.print();
	  std::cout << std::endl;

	  qv.Clamp(2);
	  std::cout << "qv mod 2 = ";
	  qv.print();
	  std::cout << std::endl;

	  std::cout << "pv = ";
	  pv.print();
	  std::cout << std::endl;

	  pv.Clamp(2);
	  std::cout << "pv mod 2 = ";
	  pv.print();
	  std::cout << std::endl;

	  std::cout << "message = ";
	  message.print();
	  std::cout << std::endl;

	  std::cout << "<qv, message> = ";
	  result_q.print();
	  std::cout << std::endl;

	  std::cout << "<pv, message> = ";
	  result_p.print();
	  std::cout << std::endl;
	  FAIL();
	  return false;
	}
      }
    }
    R_Ring_Number aa_n(4194301, 1), bb_n(262139, 1);
    aa_n[0] = -1386462;
    bb_n[0] = -86652;
    R_Ring_Number res = aa_n.Scale(4194301, 262139, 3);
    assert(res == bb_n);
    PASS();
    return true;
  }

  bool test_FHE_Switch_Keys(void) {
    int modules[] = {2, 3, 5, 7, 11};
    GLWE_Type types[2] = {LWE_Based, RLWE_Based};
    
    FHE fhe;

    for (int j = 0; j < 2; j++) {
      std::cout << "test_FHE_Switch_Keys for " << (j == 0 ? "LWE" : "RLWE") << " ";
    for (int modul_index = 0; modul_index < 5; modul_index++) {
      int modul = modules[modul_index];
    for (int s = 0; s < 30; s++) {
      FHE_Params params = fhe.Setup(2, 2, types[j], modul);
      Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> sk_pk = fhe.Key_Gen(params);
      
      long long *array_m = new long long [params[0].d];
      for (int i = 0; i < params[0].d; i++) {
	array_m[i] = R_Ring_Number::Clamp(rand(), modul);
      }
      R_Ring_Number message(modul, params[0].d, array_m);
      
      FHE_Cipher_Text c = fhe.Encrypt(params, &sk_pk.second, message);
      // check initial correction
      R_Ring_Number decoded_initial_message = c.Decrypt(params, sk_pk.first);
      if (decoded_initial_message != message.Get_Clamped(modul)) {
	std::cout << "Initial encryption is incorrect" << std::endl;
	std::cout << "Current message modul = " << modul << std::endl;
	std::cout << "attempt #" << s << std::endl;
	std::cout << "message =";
	message.print();
	std::cout << std::endl;
	
	std::cout << "decoded_initial_message =";
	decoded_initial_message.print();
	std::cout << std::endl;

	FAIL();
	return false;
      }
	
      // switch keys several times, see if we preserve correctness
      for (int t = 0; t < 1; t++) {
	Pair<R_Ring_Vector, int> c1 = c.Copy_Cipher();
	Pair<R_Ring_Vector, int> tensored_c1(R_Ring_Vector(c1.first.Get_q(), c1.first.Get_d(), (c1.first.Get_Dimension() * (c1.first.Get_Dimension() + 1)) / 2), c1.second);
	for (int i = 0; i < c1.first.Get_Dimension(); i++) {
	  tensored_c1.first[i] = c1.first[i];
	}
	R_Ring_Vector secret_tensored = sk_pk.first[c.Get_Cipher().second];
	
	secret_tensored = secret_tensored.Tensor_Product(secret_tensored);
	R_Ring_Vector secret_tensored_prime = FHE::Bit_Decomposition(secret_tensored, secret_tensored.Get_q());

	R_Ring_Vector powersof2_tensored_c1 = FHE::Powersof2(tensored_c1.first, tensored_c1.first.Get_q());
	R_Ring_Number dot1 = tensored_c1.first.Dot_Product(secret_tensored), dot2 = powersof2_tensored_c1.Dot_Product(secret_tensored_prime);

	if (dot1.Get_Clamped(dot1.Get_q()) != dot2.Get_Clamped(dot2.Get_q())) {
	  std::cout << "attempt #" << s << " level #" << t << std::endl;
	  assert(tensored_c1.first.Get_q() == secret_tensored.Get_q());

	  std::cout << "<cipher, sk> = ";
	  sk_pk.first[c.Get_Cipher().second].Dot_Product(c.Get_Cipher().first).print();
	  std::cout << std::endl;

	  std::cout << "s = ";
	  sk_pk.first[c.Get_Cipher().second].print();
	  std::cout << std::endl;

	  std::cout << "s' = ";
	  secret_tensored.print();
	  std::cout << std::endl;

	  std::cout << "<c, s'> = <";
	  tensored_c1.first.print();
	  std::cout << ", ";
	  secret_tensored.print();
	  std::cout << "> = ";
	  dot1.print();
	  std::cout << std::endl;

	  std::cout << "<c1, s''> = <";
	  powersof2_tensored_c1.print();
	  std::cout << ", ";
	  secret_tensored_prime.print();
	  std::cout << "> = ";
	  dot2.print();
	  std::cout << std::endl;

	  std::cout << 
	  std::cout << "!" << std::endl;
	  FAIL();
	  return false;
	}
	R_Ring_Vector tensored_c1_old = tensored_c1.first;
	
	c.Refresh(tensored_c1, &sk_pk.second);
	FHE_Cipher_Text tc(tensored_c1, &sk_pk.second);
	R_Ring_Number decoded_message = tc.Decrypt(params, sk_pk.first);
	if (message.Get_Clamped(modul) != decoded_message) {
	  std::cout << "modul = " << modul << std::endl;
	  std::cout << "attempt #" << s << " level #" << t << std::endl;
	  std::cout << "message =";
	  message.print();
	  std::cout << std::endl;

	  std::cout << "decoded_message =";
	  decoded_message.print();
	  std::cout << std::endl;

	  std::cout << "Modules ladder = (";
	  for (int i = 0; i < params.size(); i++) {
	    if (i != 0) {
	      std::cout << ", ";
	    }
	    std::cout << params[i].q;
	  }
	  std::cout << ")" << std::endl;

	  /*	  std::cout << "secret key = ";
	  for (int i = 0; i < sk_pk.first.size(); i++) {
	    sk_pk.first[i].print();
	    std::cout << std::endl;
	  }

	  std::cout << "pulic key = ";
	  for (int i = 0; i < sk_pk.second.size(); i++) {
	    sk_pk.second[i].print();
	    std::cout << std::endl;
	    }*/
	  std::cout << "The amount of previous noise = ";
	  c.Get_Noise(params, sk_pk.first).print();
	  std::cout << std::endl;
	  std::cout << "Cipher level = " << c.Get_Cipher().second << std::endl;
	  std::cout << "Cipher module = " << c.Get_Cipher().first.Get_q() << std::endl;

	  std::cout << "The amount of current noise = ";
	  tc.Get_Noise(params, sk_pk.first).print();
	  std::cout << std::endl;
	  std::cout << "Cipher level = " << tc.Get_Cipher().second << std::endl;
	  std::cout << "Cipher module = " << tc.Get_Cipher().first.Get_q() << std::endl;

	  std::cout << "dot1 = "; dot1.print(); std::cout << std::endl;
	  std::cout << "dot2 = "; dot2.print(); std::cout << std::endl;

	  std::cout << "powersof2_tensored_c1 = "; powersof2_tensored_c1.print(); std::cout << std::endl;

	  long long q = c.Get_Cipher().first.Get_q(), p = tc.Get_Cipher().first.Get_q();
	  
	  std::cout << "<c, s> before scale = ";
	  powersof2_tensored_c1.Dot_Product(secret_tensored_prime).print();
	  std::cout << std::endl;

	  R_Ring_Vector cc = FHE_Cipher_Text::Scale(powersof2_tensored_c1, c.Get_Cipher().first.Get_q(), tc.Get_Cipher().first.Get_q(), modul);
	  std::cout << "<c, s> after scale = ";
	  cc.Increase_Modul(c.Get_Cipher().first.Get_q());
	  cc.Dot_Product(secret_tensored_prime).Clamp(tc.Get_Cipher().first.Get_q()).print();
	  std::cout << std::endl;

	  int s_norm = 0;

	  for (int si = 0; si < secret_tensored_prime.Get_Dimension(); si++) {
	    s_norm += secret_tensored_prime[si][0];
	  }
	  
	  assert(p < q);
	  int bound = q / 2 - q / p * modul / 2 * c.Get_Cipher().first.Get_d() * s_norm;
	  std::cout << "bound = " << bound << std::endl;

	  std::cout << "s = "; secret_tensored_prime.print(); std::cout << std::endl;
	  // std::cout << "cc = "; cc.print(); std::cout << std::endl;

	  long long maxx = q * q * c.Get_Cipher().first.Get_d();
	  std::cout << "q = " << q << ", maxx = " << maxx << std::endl;
	  powersof2_tensored_c1.Increase_Modul(maxx);
	  secret_tensored_prime.Increase_Modul(maxx);

	  std::cout << "<c, s> before scale = ";
	  powersof2_tensored_c1.Dot_Product(secret_tensored_prime).print();
	  std::cout << "  ";
	  powersof2_tensored_c1.Dot_Product(secret_tensored_prime).Clamp(modul).print();
	  std::cout << std::endl;

	  std::cout << "<c, s> after scale = ";
	  cc.Increase_Modul(maxx);
	  cc.Dot_Product(secret_tensored_prime).print();
	  std::cout << "   ";
	  cc.Dot_Product(secret_tensored_prime).Clamp(modul).print();
	  std::cout << std::endl;

	  std::cout << "cscale = ";
	  cc.print();
	  std::cout << std::endl;

	  std::cout << "c1 = "; c1.first.print(); std::cout << std::endl;
	  std::cout << "tensored_c1 = "; tensored_c1_old.print(); std::cout << std::endl;

	  FAIL();
	  return false;
	}
	c = tc;
      }

      delete [] array_m;
    }
    }
    PASS();
    std::cout << std::endl;}
    return true;
  }

  bool test_FHE_One_Add(GLWE_Type type, FHE_Operation operation_type, int modul = 2) {
    FHE fhe;

    for (int ii = 0; ii < 30; ii++) {
      FHE_Params params = fhe.Setup(2, 2, type, modul);
    Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> sk_pk = fhe.Key_Gen(params);
    
    long long *array_m[2];
    for (int j = 0; j < 2; j++) {
      array_m[j] = new long long [params[0].d];
      for (int i = 0; i < params[0].d; i++) {
	array_m[j][i] = R_Ring_Number::Clamp(rand(), modul);
      }
    }
        
    R_Ring_Number message1(modul, params[0].d, array_m[0]), message2(modul, params[0].d, array_m[1]), message3(modul, params[0].d);
    if (operation_type == FHE_Addition) {
      message3 = (message1 + message2);
      message3.Clamp(modul);
    } else if (operation_type == FHE_Multiplication) {
      message3 = (message1 * message2);
      message3.Clamp(modul);
    }
     
    FHE_Cipher_Text c1 = fhe.Encrypt(params, &sk_pk.second, message1), c2 = fhe.Encrypt(params, &sk_pk.second, message2);
    R_Ring_Number m1 = c1.Decrypt(params, sk_pk.first), m2 = c2.Decrypt(params, sk_pk.first);
    if (m1 != message1 || m2 != message2) {
      std::cout << "m1 = "; m1.print(); std::cout << std::endl;
      std::cout << "message1 = "; message1.print(); std::cout << std::endl;
      std::cout << "m2 = "; m2.print(); std::cout << std::endl;
      std::cout << "message2 = "; message2.print(); std::cout << std::endl;
      std::cout << "Initial decryption failed" << std::endl;
      FAIL();
      return false;
    }
    FHE_Cipher_Text res_c;
    if (operation_type == FHE_Addition) {
      res_c = c1 + c2;
    } else if (operation_type == FHE_Multiplication) {
      res_c = c1 * c2;
    }
    
    R_Ring_Number decoded_message = res_c.Decrypt(params, sk_pk.first);
    if (message3 != decoded_message) {
      std::cout << "attempt #" << ii << std::endl;
      std::cout << "modul " << modul << std::endl;

      std::cout << "message1 " << ((operation_type == FHE_Addition) ? "+" : "*") << " message2 = ";
      message1.print();
      std::cout << " " << ((operation_type == FHE_Addition) ? "+" : "*") << " ";
      message2.print();
      std::cout << std::endl;
      
      std::cout << "message = ";
      message3.print();
      std::cout << std::endl;

      std::cout << "decoded_message = ";
      decoded_message.print();
      std::cout << std::endl;

      std::cout << "c1 = ";
      c1.print();
      std::cout << std::endl;

      std::cout << "c2 = ";
      c2.print();
      std::cout << std::endl;

      std::cout << "final amount of noise = mod " << res_c.Get_Cipher().first.Get_q() << " = ";
      res_c.Get_Cipher().first.Dot_Product(sk_pk.first[res_c.Get_Cipher().second]).print();
      std::cout << std::endl;

      std::cout << "original amount of noise = <c1, s1> mod " << c1.Get_Cipher().first.Get_q() << " = ";
      c1.Get_Cipher().first.Dot_Product(sk_pk.first[c1.Get_Cipher().second]).print();
      std::cout << std::endl;

      std::cout << "original amount of noise = <c2, s2> mod " << c2.Get_Cipher().first.Get_q() << " = ";
      c2.Get_Cipher().first.Dot_Product(sk_pk.first[c2.Get_Cipher().second]).print();
      std::cout << std::endl;

      std::cout << "Modules ladder = (";
      for (int j = 0; j < params.size(); j++) {
	if (j != 0) { std::cout << ", "; }
	params[j].print();
      }
      std::cout << ")" << std::endl;

      std::cout << "Current module = " << params[res_c.Get_Cipher().second].q << std::endl;

      std::cout << "sk = ";
      sk_pk.first[c1.Get_Cipher().second].print();
      std::cout << std::endl;

      std::cout << "c1.q = " << c1.Get_Cipher().first.Get_q() << std::endl;
      std::cout << "c2.q = " << c2.Get_Cipher().first.Get_q() << std::endl;
      int c1_s_dot = 0;
      for (int i = 0; i < c1.Get_Cipher().first.Get_Dimension(); i++) {
	c1_s_dot += c1.Get_Cipher().first[i][0] * sk_pk.first[c1.Get_Cipher().second][i][0];
      }
      std::cout << "<c1, s> = " << c1_s_dot << " = ";
      c1.Get_Cipher().first.Dot_Product(sk_pk.first[c1.Get_Cipher().second]).print();
      std::cout << std::endl;
      int c2_s_dot = 0;
      for (int i = 0; i < c2.Get_Cipher().first.Get_Dimension(); i++) {
	c2_s_dot += c2.Get_Cipher().first[i][0] * sk_pk.first[c1.Get_Cipher().second][i][0];
      }
      std::cout << "<c2, s> = " << c2_s_dot << " = ";
      c2.Get_Cipher().first.Dot_Product(sk_pk.first[c2.Get_Cipher().second]).print();
      std::cout << std::endl;
      R_Ring_Vector sk_tensored = sk_pk.first[c1.Get_Cipher().second].Tensor_Product(sk_pk.first[c1.Get_Cipher().second]);
      int dimension = c1.Get_Cipher().first.Get_Dimension();
      int new_dimension = (dimension * (dimension + 1)) / 2;
      R_Ring_Vector c3(c1.Get_Cipher().first.Get_q(), c1.Get_Cipher().first.Get_d(), new_dimension);
      if (operation_type == FHE_Multiplication) {
	int index = 0;
	for (int i = 0; i < dimension; i++) {
	  c3[index++] = c1.Get_Cipher().first[i] * c2.Get_Cipher().first[i];
	  for (int j = i + 1; j < dimension; j++) {
	    c3[index++] = c1.Get_Cipher().first[i] * c2.Get_Cipher().first[j] + c1.Get_Cipher().first[j] * c2.Get_Cipher().first[i];
	  }
	}
      } else if (operation_type == FHE_Addition) {
	for (int i = 0; i < dimension; i++) {
	  c3[i] = c1.Get_Cipher().first[i] + c2.Get_Cipher().first[i];
	}
      }

      std::cout << "<sk_tensored, c3> mod " << c1.Get_Cipher().first.Get_q() << " = <";
      sk_tensored.print();
      std::cout << ", ";
      c3.print();
      std::cout << "> = ";
      R_Ring_Number before_reduction = sk_tensored.Dot_Product(c3);
      before_reduction.print();
      std::cout << std::endl;

      R_Ring_Number my_number(modul, 1);
      my_number[0] = LLONG_MAX;
      // 9223372036854775807
      // 2078548109380307949
      // 32477229033819617
      std::cout << "LLONG_MAX = "; my_number.print(); std::cout << std::endl;
      

      int p = res_c.Get_Cipher().first.Get_q(), q = c3.Get_q();

      R_Ring_Vector c_powers2 = FHE::Powersof2(c3, c3.Get_q());
      assert(c3.Get_q() == sk_tensored.Get_q());
      R_Ring_Vector sk_tensored_bitdec = FHE::Bit_Decomposition(sk_tensored, c3.Get_q());
      std::cout << "<c_powers2, sk_tensored_bitdec> = ";
      c_powers2.Dot_Product(sk_tensored_bitdec).print();
      std::cout << " = ";
      c_powers2.Dot_Product(sk_tensored_bitdec).Clamp(modul).print();
      std::cout << " mod " << modul;
      std::cout << std::endl;

      R_Ring_Vector c_scale = res_c.Scale(c_powers2, q, p, modul);
      for (int i = 0; i < c_scale.Get_Dimension(); i++) {
	assert(R_Ring_Number::Clamp(c_scale[i][0], modul) == R_Ring_Number::Clamp(c_powers2[i][0], modul));
      }
      std::cout << "Scale(" << c3.Get_q() << ", " << res_c.Get_Cipher().first.Get_q() << ", " << modul << ")" << std::endl;
      std::cout << "<c_scale, sk_tensored_bit_dec> = ";
      c_scale.Increase_Modul(c3.Get_q());
      c_scale.Dot_Product(sk_tensored_bitdec).print();
      std::cout << std::endl;

      for (int i = 0; i < c_powers2.Get_Dimension(); i++) {
	assert((c_scale[i] * sk_tensored_bitdec[i]).Get_Clamped(modul) ==
	       (c_powers2[i] * sk_tensored_bitdec[i]).Get_Clamped(modul));
      }
      
      R_Ring_Number r1(c3.Get_q(), c_powers2.Get_d()), r2(c3.Get_q(), c_scale.Get_d());

      std::cout << "q = " << q << " p = " << p << std::endl;
      long long r1_int = 0, r2_int = 0;
      for (int i = 0; i < c_powers2.Get_Dimension(); i++) {
	R_Ring_Number r1_new, r2_new;
	R_Ring_Number a1 = c_powers2[i] * sk_tensored_bitdec[i], a2 = c_scale[i] * sk_tensored_bitdec[i];
	r1_new = r1 + a1;
	r2_new = r2 + a2;
	r1_int += a1[0];
	r2_int += a2[0];
	
	/* if (r1_new.Get_Clamped(modul) != r2_new.Get_Clamped(modul)) {
	  std::cout << "r1_new != r2_new in position #" << i << std::endl;
	  std::cout << "r1_new = ";

	  r1_new.print();
	  std::cout << std::endl << "r2_new = ";
	  r2_new.print();
	  std::cout << std::endl;

	  std::cout << "r1 = ";
	  r1.print();
	  std::cout << std::endl << "r2 = ";
	  r2.print();
	  std::cout << std::endl;

	  std::cout << "a1 = ";
	  a1.print();
	  std::cout << std::endl << "r2 = ";
	  a2.print();
	  std::cout << std::endl;

	  break;
	  } */
	r1 = r1_new;
	r2 = r2_new;
      }
      std::cout << "r1_int = " << r1_int << std::endl;
      std::cout << "r2_int = " << r2_int << std::endl;

      std::cout << "r1 = ";
      r1.print();
      std::cout << std::endl;
      std::cout << "r2 = ";
      r2.print();
      std::cout << std::endl;

      R_Ring_Vector mod_p_zero_vector = (c_scale - c_powers2).Get_Clamped(modul);
      for (int i = 0; i < mod_p_zero_vector.Get_Dimension(); i++) {
	assert(mod_p_zero_vector[i][0] == 0);
      }

      std::cout << "mod_p_zero_vector = ";
      mod_p_zero_vector.print();
      std::cout << std::endl;

      std::cout << "c_sclae = ";
      c_scale.print();
      std::cout << std::endl << "c_powers2 = ";
      c_powers2.print();
      std::cout << std::endl << "sk_tensored_bitdec = ";
      sk_tensored_bitdec.print();
      std::cout << std::endl;

      R_Ring_Vector c_powers2_clamped = c_powers2.Get_Clamped(modul), c_scale_clamped = c_scale.Get_Clamped(modul);
      std::cout << "c_powers2_clamped = ";
      c_powers2_clamped.print();
      std::cout << std::endl << "c_scale_clamped   = ";
      c_scale_clamped.print();
      std::cout << std::endl;
      if (c_powers2_clamped != c_scale_clamped) {
	std::cout << "ALERT c_powers2_clamped != c_scale_clamped" << std::endl;
      }

      std::cout << "Result = ";
      before_reduction.Clamp(modul);
      before_reduction.print();
      std::cout << std::endl;

      R_Ring_Number result = res_c.Get_Cipher().first.Dot_Product(sk_pk.first[res_c.Get_Cipher().second]);

      std::cout << "<s, c> = ";
      result.print();
      std::cout << " = ";
      result.Clamp(result.Get_q());
      result.print();
      std::cout << " mod " << result.Get_q();
      std::cout << " = ";
      result.Clamp(modul);
      result.print();
      std::cout << " mod " << modul;
      std::cout << std::endl;
      

      //      FHE::Refresh(c3, sk_pk.second);

      FAIL();
      return false;
    }
    }
    PASS();
    return true;
  }

  bool test_FHE_Operations(void) {
    int modules[] = {2, 3, 5, 7, 31};
    char tests_names[][27] = {"test_FHE_One_Add_for_LWE  ",
			    "test_FHE_One_Add_for_RLWE ",
			    "test_FHE_One_Mult_for_LWE ",
			    "test_FHE_One_Mult_for_RLWE"};
    GLWE_Type types[] = {LWE_Based, RLWE_Based, LWE_Based, RLWE_Based};
    FHE_Operation operation[] = {FHE_Addition, FHE_Addition, FHE_Multiplication, FHE_Multiplication};

    for (int i = 0; i < 4; i++) {
      std::cout << tests_names[i] << " ";
    for (int j = 0; j < 5; j++) {
      std::cout << std::endl << "Module " << modules[j] << std::endl;
      if (!test_FHE_One_Add(types[i], operation[i], modules[j])) {
	return false;
      }
    }
    }
    return true;
  }

  bool test_FHE_LSS() {
    std::cout << "test_FHE_LSS ";/*
    FHE fhe;
    FHE_Params params = fhe.Setup(3, 2, type);
    int noof_vectors = 5;
    Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> sk_pk = fhe.Key_Gen(params);
    
    long long *array_m[5];
    for (int j = 0; j < 5; j++) {
      array_m[j] = new long long [params[0].d];
      for (int i = 0; i < params[0].d; i++) {
	array_m[j][i] = rand() % 2;
      }
    }
    
    R_Ring_Number message1(2, params[0].d, array_m[0]), message2(2, params[0].d, array_m[1]), message3(2, params[0].d);
    if (operation_type == FHE_Addition) {
      message3 = message1 + message2;
    } else if (operation_type == FHE_Multiplication) {
      message3 = message1 * message2;
    }
     
    FHE_Cipher_Text c1 = fhe.Encrypt(params, &sk_pk.second, message1), c2 = fhe.Encrypt(params, &sk_pk.second, message2);*/
    PASS();
    return true;
  }

public:
  void Run_Tests() {
    if (!test_Zero_Number() ||
	!test_Nonzero_Number() ||
	!test_Number_Addition() ||
	!test_Number_Multiplication() ||
	!test_Number_Equality() ||
	!test_Number_Reduction() ||

	!test_GLWE_for_LWE() ||
	!test_GLWE_for_RLWE() ||

	!test_FHE_for_LWE() ||
	!test_FHE_for_RLWE() ||
	!test_Powersof2_BitDecomposition() ||
	!test_Number_Scale() ||
	!test_FHE_Switch_Keys() ||
	!test_FHE_Operations() ||
	!test_FHE_LSS()) {
      std::cout << "Overall tests FAILED" << std::endl;
    } else {
      std::cout << "Overall tests PASSED" << std::endl;
    }
  }
};

int main (int argc, char * const argv[]) {
  std::cout << std::endl;

  Test tests;
  tests.Run_Tests();

  /*  std::cout << "-1 % 2 == " << (-1) % 2 << std::endl;
  std::cout << "-2 % 2 == " << (-2) % 2 << std::endl;
  std::cout << "-3 % 2 == " << (-3) % 2 << std::endl;

  std::cout << "(int)(-3.5) == " << (int)(-3.5) << std::endl;
  std::cout << "(int)(3.5) == " << (int)(3.5) << std::endl; */
  return 0;	
}
