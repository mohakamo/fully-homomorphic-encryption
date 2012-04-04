#include "FHE.h"
#include "LSS.h"
#include "Tests_SI_HE.h"
#include "Timing.h"
#include "Tests_ZZX_Ring_Number.h"
#include <iostream>
#include <time.h>

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
  void Setup(int dimension) {
    type = RLWE_Based;
    ZZ modul = to_ZZ(2);
    Regev_Params temp_params = fhe.Setup(L, type, modul);
    params = temp_params;

    sk = fhe.Secret_Key_Gen(params);
    pk = fhe.Public_Key_Gen(params, sk);
    evalk = fhe.Evaluation_Key_Gen(params, sk, pk);
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

  /*
  Test tests;
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

  Timing_LSS timingLSS;
  const int max_dimension = 6; // till (20, 14) for dimension = 3, modul = 10, noise bound not found
  //  int elements_modul[] = {1, 127, 2, 3, 4, 5, 6};

  for (int dimension = 2; dimension <= max_dimension; dimension++) {
    for (int modul_i = 0; modul_i < 1; modul_i++) {
      //      int modul = elements_modul[modul_i];
      timingLSS.Setup(dimension);
      std::cout << "dimension = " << dimension << "\n";// << ", modul = " << modul << std::endl;
      timingLSS.Run_LSS(dimension);
      std::cout << std::endl;
    }
  }
  
  return 0;	
}
