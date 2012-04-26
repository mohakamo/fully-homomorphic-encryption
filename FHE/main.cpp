#include "FHE.h"
#include "LSS.h"
#include "Tests_SI_HE.h"
#include "Timing.h"
#include "Tests_ZZX_Ring_Number.h"
#include <iostream>
#include <time.h>

using namespace std;

class Timing_Mult {
  SI_HE fhe;
  Regev_Params params;
  SI_HE_Secret_Key_Type sk;
  SI_HE_Public_Key_Type pk;
  SI_HE_Evaluation_Key_Type evalk;
  ZZ modul;
  int L;
  GLWE_Type type;

  int my_bound;

public:
  Timing_Mult() : L(1) {}
  /**
   * Setup SI_HE scheme
   * @param dimension - vector dimension
   **/
  bool Setup(bool prompt = false) {
    clock_t start, total_start = clock();

    type = RLWE_Based;
    modul = to_ZZ(1) << 16;
    Regev_Params temp_params = fhe.Setup(L, type, modul);
    params = temp_params;

    long long size_q = NumBits(params.q);
    long long multiplicative_constant = params.d * params.n * size_q / 8; 
		long long exp_evalk_size = (params.n + 1) * (params.n + 1) * params.n * size_q * size_q * size_q * multiplicative_constant;
		cout << "Estimated |evalk| = " << exp_evalk_size / 1024 / 1024 << " MB\n";
		char c = 'y';

		start = clock();
    sk = fhe.Secret_Key_Gen(params);
		if (prompt) {
			cout << "Secret key generates in time: " << (clock() - start) / (double)CLOCKS_PER_SEC << "\nWant to continue? (y/n)";
			cin >> c;
			if (c == 'n' || c == 'N') return 0;
		}

		start = clock();
    pk = fhe.Public_Key_Gen(params, sk);
		if (prompt) {
			cout << "Public key generates in time: " << (clock() - start) / (double)CLOCKS_PER_SEC << "\nWant to continue? (y/n)";
			cin >> c;
			if (c == 'n' || c == 'N') return 0;
		}

		start = clock();
    evalk = fhe.Evaluation_Key_Gen(params, sk, pk);
		if (prompt) {
			cout << "Evaluation key generates in time: " << (clock() - start) / (double)CLOCKS_PER_SEC << "\nWant to continue? (y/n)";
			cin >> c;
			if (c == 'n' || c == 'N') return 0;
		}

    long long evalk_size = evalk.size() * evalk[0].Get_Noof_Columns() * evalk[0].Get_Noof_Rows() * multiplicative_constant;
    cout << "|evalk| = " << evalk_size / 1024 / 1024 << " MB\n";

    cout << "Total time for key generation = " <<  (clock() - total_start) / (double)CLOCKS_PER_SEC << std::endl;
		return 1;
  }

  /**
   * Run_LSS function runs LSS on randomly generated data with parameters of FHE scheme that are adjusted to 
   *  get better performance
   * @param dimension - the dimension of the vectors being randomly generated
   * @param bound - bound on absolute value of vector's elements
   * @return time spent on computing LSS, not including scheme setup
   **/
  void Run_Mult(void) {
    clock_t start, total_start = clock();
    double time = 0;

    std::vector<R_Ring_Number> messages_x, messages_y;
    std::vector<SI_HE_Cipher_Text> c_x, c_y;
    start = clock();
    ZZ sub_entr = (modul % 2 == 0) ? (modul - 2) / 2 : (modul - 1) / 2;
		R_Ring_Number mx(modul, params.d);
		mx[0] = RandomBnd(modul) - sub_entr;//R_Ring_Number::Uniform_Rand(modul, params.d);
		R_Ring_Number my(modul, params.d);
		my[0] = RandomBnd(modul) - sub_entr;//R_Ring_Number::Uniform_Rand(modul, params.d);

		cout << "Calculating: " << mx[0] << " * " << my[0] << "\n";

		SI_HE_Cipher_Text cx = fhe.Encrypt(params, &pk, mx, &evalk);
		SI_HE_Cipher_Text cy = fhe.Encrypt(params, &pk, my, &evalk);
		std::cout << "Size of cipher = " << cx.getSizeInBytes() << "\n";
		//		std::cout << "Size of evalk = " << evalk.Get_Dimension() * evalk[0].Get_Noof_Columns() * evalk[0].Get_Noof_Rows() * NumBits(
		cx.Add_Secret_Key_Info(&sk);
		cy.Add_Secret_Key_Info(&sk);

    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to encrypt = " << time << std::endl;
  
    start = clock();
    SI_HE_Cipher_Text res_c = cx * cy;

    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to compute on the cipher = " << time << std::endl;

    start = clock();
		R_Ring_Number res = mx * my;
    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to compute in the clear = " << time << std::endl;

    start = clock();
    R_Ring_Number res_d = res_c.Decrypt(params, sk);
    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to decrypt = " << time << std::endl;

    if (res != res_d) {
      std::cout << "FAILUR!" << std::endl;
      std::cout << res << " != " << res_d << "\n";
    } else {
      std::cout << "PASS!" << std::endl;
      std::cout << res << " == " << res_d << "\n";
    }

    std::cout << "Total time = " <<  (clock() - total_start) / (double)CLOCKS_PER_SEC << std::endl;
  }  
};

class Timing_LSS {
  SI_HE fhe;
  Regev_Params params;
  SI_HE_Secret_Key_Type sk;
  SI_HE_Public_Key_Type pk;
  SI_HE_Evaluation_Key_Type evalk;
  ZZ modul;
  int L;
  GLWE_Type type;

  int my_max_dimension;
  int my_bound;

public:
  Timing_LSS() : L(2) {}
  /**
   * Setup SI_HE scheme
   * @param dimension - vector dimension
   **/
  void Setup(int dimension, int modul_i, int n = -1) {
    type = RLWE_Based;
    modul = to_ZZ(modul_i);
    Regev_Params temp_params = fhe.Setup(1, type, modul, n);
    params = temp_params;

    long long size_q = NumBits(params.q);
    long long multiplicative_constant = params.d * params.n * size_q / 8; 
    long long exp_evalk_size = (params.n + 1) * (params.n + 1) * params.n * size_q * size_q * size_q * multiplicative_constant;

    cout << "Estimated |evalk| = " << exp_evalk_size / 1024 / 1024 << " MB\n";

    sk = fhe.Secret_Key_Gen(params);
    long long sk_size = sk.size() * sk[0].Get_Dimension() * multiplicative_constant;
    cout << "|sk| = " << sk_size << " B\n";

    pk = fhe.Public_Key_Gen(params, sk);
    long long pk_size = pk.Get_Noof_Columns() * pk.Get_Noof_Rows() * multiplicative_constant;
    cout << "|pk| = " << pk_size / 1024 << " B\n";

    evalk = fhe.Evaluation_Key_Gen(params, sk, pk);
    long long evalk_size = evalk.size() * evalk[0].Get_Noof_Columns() * evalk[0].Get_Noof_Rows() * multiplicative_constant;
    cout << "|evalk| = " << evalk_size / 1024 / 1024 << " MB\n";
  }

  void Run_Mult(void) {
    clock_t start;
    double time;
    
    cout << "modul = " << modul << "\n";
    R_Ring_Number mx = R_Ring_Number::Uniform_Rand(modul, params.d);
    R_Ring_Number my = R_Ring_Number::Uniform_Rand(modul, params.d);

    start = clock();
    SI_HE_Cipher_Text cx = fhe.Encrypt(params, &pk, mx, &evalk);
    SI_HE_Cipher_Text cy = fhe.Encrypt(params, &pk, my, &evalk);
    cx.Add_Secret_Key_Info(&sk);
    cy.Add_Secret_Key_Info(&sk);

    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    cout << "Time to encrypt = " << time << endl;
  
    start = clock();
    SI_HE_Cipher_Text cres = cx * cy;
    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    cout << "Time to mult on cipher = " << time << endl;
    
    start = clock();
    R_Ring_Number mres = mx * my;
    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    cout << "Time to mult in clear = " << time << endl;

    R_Ring_Number res = cres.Decrypt(params, sk);

    if (res != mres) {
      cout << "FAILUR!" << endl;
      cout << "(" << res << ") != (" << mres << ")\n";
    } else {
      cout << "PASS!" << endl;
      cout << "(" << res << ")\n";
    }
  }

  /**
   * Run_LSS function runs LSS on randomly generated data with parameters of FHE scheme that are adjusted to 
   *  get better performance
   * @param dimension - the dimension of the vectors being randomly generated
   * @param bound - bound on absolute value of vector's elements
   * @return time spent on computing LSS, not including scheme setup
   **/
  void Run_LSS(int dimension) {
    clock_t start, total_start = clock();
    double time = 0;

    std::vector<R_Ring_Number> messages_x, messages_y;
    std::vector<SI_HE_Cipher_Text> c_x, c_y;
    start = clock();
    cout << "modul = " << modul << "\n";
    for (int j = 0; j < dimension; j++) {
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
    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to encrypt = " << time << std::endl;
  
    start = clock();
    Pair<SI_HE_Cipher_Text, SI_HE_Cipher_Text> res_c = Compute_LSS<SI_HE_Cipher_Text>(c_x, c_y);

    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to compute on the cipher = " << time << std::endl;

    start = clock();
    Pair<R_Ring_Number, R_Ring_Number> res = Compute_LSS<R_Ring_Number>(messages_x, messages_y);
    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to compute in the clear = " << time << std::endl;

    start = clock();
    Pair<R_Ring_Number, R_Ring_Number> res_d(res_c.first.Decrypt(params, sk), res_c.second.Decrypt(params, sk));
    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to decrypt = " << time << std::endl;

    if (res.first != res_d.first || res.second != res_d.second) {
      std::cout << "FAILUR!" << std::endl;
      std::cout << "(" << res.first << ", " << res.second << ") != (" << res_d.first << ", " << res_d.second <<	")\n";
    } else {
      std::cout << "PASS!" << std::endl;
      std::cout << "(" << res.first << ", " << res.second << ") == (" << res_d.first << ", " << res_d.second << ")\n";
    }

    std::cout << "Total time = " <<  (clock() - total_start) / (double)CLOCKS_PER_SEC << std::endl;
  }  
};

int main (int argc, char * const argv[]) {
  std::cout << std::endl;
  /*
  std::cout << "(int)(-1.999) = " << (int)(-1.999) << std::endl;
  std::cout << "(int)(-0.0002) = " << (int)(-0.0002) << std::endl;

  std::cout << "(ZZ)(1 / 16) = " << to_ZZ(1 + 8) / to_ZZ(16) << std::endl;
  std::cout << "(ZZ)(9 / 16) = " << to_ZZ(9 + 8) / to_ZZ(16) << std::endl;
  std::cout << "(ZZ)(-1 / 16) = " << to_ZZ(-1 + 8) / to_ZZ(16) << std::endl;
  std::cout << "(ZZ)(-9 / 16) = " << to_ZZ(-9 + 8) / to_ZZ(16) << std::endl;
  return 0;
  */

  /*
  std::cout << "Norm distr samples: ";
  for (int i = 0; i < 20; i++) {
    std::cout << NormDistr::sample_standard(0, 8) << " ";
  }
  std::cout << std::endl;
  */

  /*  Test tests;
  tests.Run_Tests();
  return 0;
  */

  /*
  std::cout << "check new ";
  long long nn = 371558400 * sizeof(R_Ring_Number);
  char *a = new char [nn];
  delete [] a;
  std::cout << " success";
  std::cout << "\ncheck new 2 ";
  long long n = 30963200;
  R_Ring_Number *arr = new R_Ring_Number[n];
  delete [] arr;
  std::cout << " success\n";
  return 0;
  */

  /*
  SI_HE fhe;
  clock_t start;
  
  GLWE_Type type = RLWE_Based;
  for (int i = 3; i < 10; i++) {
    int n = 1 << i;
    Regev_Params params = fhe.Setup(2, type, ZZ(INIT_VAL, 2), n);
    start = clock();
    SI_HE_Secret_Key_Type sk = fhe.Secret_Key_Gen(params);
    SI_HE_Public_Key_Type pk = fhe.Public_Key_Gen(params, sk);
    SI_HE_Evaluation_Key_Type evalk = fhe.Evaluation_Key_Gen(params, sk, pk);
    std::cout << "n = " << n << ", time = " << (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "|sk| = " << sk.size() * sk[0].Get_Dimension() << std::endl;
    std::cout << "|pk| = " << pk.Get_Noof_Columns() * pk.Get_Noof_Rows() << std::endl;
    std::cout << "|evalk| = " << evalk.size() * evalk[0].Get_Noof_Columns() * evalk[0].Get_Noof_Rows() << std::endl;
    std::cout << "****************************" << std::endl;
  }

  return 0;
  */

  /*  Timing_Operations timing;
  timing.TimingOperations();
  return 0;
  */

  /*
  Timing_LSS timingLSS;
  const int max_dimension = 6; // till (20, 14) for dimension = 3, modul = 10, noise bound not found
  int elements_modul[] = {2, 4, 8, 16, 32};

  for (int dimension = 2; dimension <= max_dimension; dimension++) {
    for (int modul_i = 0; modul_i < 2; modul_i++) {
      int modul = elements_modul[modul_i];
      timingLSS.Setup(dimension, modul);
      std::cout << "dimension = " << dimension << ", modul = " << modul << std::endl << std::flush;
      timingLSS.Run_LSS(dimension);
      std::cout << std::endl;
    }
  }
  */
  
	/*
  Timing_LSS timingLSS;
  int ns[] = {8, 16, 32, 64, 128, 256};

  for (int i = 0; i < 6; i++) {
    int n = ns[i];
    cout << "###############    " << n << "   ##############\n";
    timingLSS.Setup(2, 16, n);
    timingLSS.Run_Mult();
    std::cout << std::endl;
  }
*/

	/*	long long size_q = NumBits(to_ZZ(1) << 16);
	long long multiplicative_constant = 2048 * size_q / 8; 
	long long exp_evalk_size = 4 * 16 * 16 * 16 * 2048 * 16 / 8;
	cout << "Real paramters |evalk| = " << exp_evalk_size / 1024 / 1024 << " MB\n";
	*/

	Timing_Mult timing;
	if (timing.Setup(true))	timing.Run_Mult();

  return 0;	
}
