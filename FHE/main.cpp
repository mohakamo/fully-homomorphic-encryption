// Define below comments all assertions (hopefully) :), no debug information
#define SLASH
#define assert SLASHSLASH assert

#include "FHE.h"
#include "LSS.h"
#include "Tests.h"
#include <iostream>
#include <time.h>

class Timing_LSS {
  FHE fhe;
  FHE_Params params;
  FHE_Secret_Key_Type sk;
  FHE_Public_Key_Type pk;
  ZZ modul;
  int L;
  GLWE_Type type;

  int my_max_dimension;
  int my_bound;

  void IncreaseL() {
    L++;
    Setup(my_max_dimension, my_bound);
  }
public:
  Timing_LSS() : L(3) {}
  /**
   * Setup FHE scheme
   * @param max_dimension - maximum desired dimension
   * @param bound - bound on absolute value of vector's elements (required to choose module ladder)
   **/
  void Setup(int max_dimension, int bound, bool justPrintParams = false) {
    type = RLWE_Based;
    if (bound % 2 == 0) bound = bound / 2;
    else bound = (bound - 1) / 2;
    my_max_dimension = max_dimension;
    my_bound = bound;

    // choosing message module
    ZZ zz_dimension = ZZ(INIT_VAL, max_dimension);
    ZZ zz_bound = ZZ(INIT_VAL, bound);
    ZZ lower_bound_on_modul;
    lower_bound_on_modul = 2 * zz_dimension * zz_dimension * zz_bound * zz_bound * zz_bound + 1;
    //    std::cout << "lower_bound_on_modul = " << lower_bound_on_modul << std::endl;
    //    std::cout << "upper_bound_on_modul = " << 2 * lower_bound_on_modul << std::endl;
    ZZ upper_bound_on_modul;
    upper_bound_on_modul = 2 * lower_bound_on_modul; // by Bertrand's postulate we should find a prime in this range, otherwise something is going wrong....
    modul = lower_bound_on_modul;
    while (!ProbPrime(modul) && modul < upper_bound_on_modul) {
      modul++;
    }
    if (!ProbPrime(modul)) {
      std::cout << "Modul not found" << std::endl;
      exit(1);
    }
    //std::cout << "Message modul chosen = " << modul << std::endl;
    if (justPrintParams) {
      fhe.Print_Possible_Parameters(L, type, modul);
      return;
    }

    FHE_Params temp_params = fhe.Setup(3, L, type, modul);
    //    std::cout << "temp_params.size() = " << temp_params.size() << std::endl;
    params = temp_params;
    
    
    Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> sk_pk = fhe.Key_Gen(params);
    sk = sk_pk.first;
    pk = sk_pk.second;

    //    std::cout << "params.size() = " << params.size() << std::endl;
  }

  /**
   * Run_LSS function runs LSS on randomly generated data with parameters of FHE scheme that are adjusted to 
   *  get better performance
   * @param dimension - the dimension of the vectors being randomly generated
   * @param bound - bound on absolute value of vector's elements
   * @return time spent on computing LSS, not including scheme setup
   **/
  void Run_LSS(int dimension, int bound) {
    clock_t start, total_start = clock();
    double time = 0;

    std::cout << "params.size() = " << params.size() << std::endl;

    std::vector<R_Ring_Number> messages_x, messages_y;
    std::vector<FHE_Cipher_Text> c_x, c_y;
    start = clock();
    R_Ring_Number mx(modul, params[0].d), my(modul, params[0].d);
    for (int j = 0; j < dimension; j++) {
      mx = RandomBnd(ZZ(INIT_VAL, bound));
      my = RandomBnd(ZZ(INIT_VAL, bound));
      messages_x.push_back(mx);
      messages_y.push_back(my);
      c_x.push_back(fhe.Encrypt(params, &pk, messages_x[j]));
      assert(c_x[j].ThNoise != 0);
      c_x[c_x.size() - 1].Add_Secret_Key_Info(&sk);
      c_y.push_back(fhe.Encrypt(params, &pk, messages_y[j]));
      assert(c_y[j].ThNoise != 0);
      c_y[c_y.size() - 1].Add_Secret_Key_Info(&sk);
    }
    time = (clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Time to encrypt = " << time << std::endl;
  
    start = clock();
    Pair<FHE_Cipher_Text, FHE_Cipher_Text> res_c = Compute_LSS<FHE_Cipher_Text>(c_x, c_y);

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
      std::cout << "params = ";
      for (int i = 0; i < params.size(); i++) {
	params[i].print();
	if (i + 1 != params.size()) std::cout << ", ";
      }
      std::cout << std::endl;
      
      std::cout << "sk dimensions = (";
      for (int i = 0; i < sk.size(); i++) {
	std::cout << sk[i].Get_Dimension();
	if (i + 1 != sk.size()) {
	  std::cout << ", ";
	}
      }
      std::cout << ")" << std::endl;
      
      std::cout << "pk dimensions = (";
      for (int i = 0; i < pk.size(); i++) {
	std::cout << "(" << pk[i].Get_Noof_Rows() << ", " << pk[i].Get_Noof_Columns() << ")";
	if (i + 1 != pk.size()) {
	  std::cout << ", ";
	}
      }
      std::cout << ")" << std::endl;
    }

    std::cout << "Total time = " <<  (clock() - total_start) / (double)CLOCKS_PER_SEC << std::endl;
  }  
};

int main (int argc, char * const argv[]) {
  std::cout << std::endl;

  /*
  Test tests;
  tests.Run_Tests();
  return 0;
  */

  Timing_LSS timingLSS;
  const int max_dimension = 10; // till (20, 14) for dimension = 3, modul = 10, noise bound not found
  int elements_modul[] = {3, 5, 7, 9, 11, 13};

  for (int dimension = 5; dimension <= max_dimension; dimension++) {
    for (int modul_i = 0; modul_i < 6; modul_i++) {
      int modul = elements_modul[modul_i];
      timingLSS.Setup(dimension, modul, false);
      std::cout << "dimension = " << dimension << ", modul = " << modul << std::endl;
      timingLSS.Run_LSS(modul, dimension);
      std::cout << std::endl;
    }
  }
  
  return 0;	
}
