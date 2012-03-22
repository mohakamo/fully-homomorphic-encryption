#ifndef _TIMING_H_
#define _TIMING_H_

#include "FHE.h"
#include <iostream>
#include <fstream>
#include <time.h>

class Timing_Operations {
  FHE fhe;
  FHE_Params params;
  FHE_Secret_Key_Type sk;
  FHE_Public_Key_Type pk;
  ZZ modul;
  int L;
  GLWE_Type type;

public:
  /**
   * Printing out parameters of the scheme: q0, mu
   * @param t - message modul, should be prime
   * @param L - number of levels of multiplication
   * @param n
   */
  void PrintParameters(int t, int L, int n) {
    type = RLWE_Based;

    params = fhe.Setup(0, L, type, to_ZZ(t), n);
    
    /*    std::cout << "params = ";
    for (int i = 0; i < params.size(); i++) {
      params[i].print();
      if (i + 1 != params.size()) std::cout << ", ";
    }
    std::cout << std::endl;*/
  }

  void FillTable() {
    int prms[][2] = {
      {1, 512},
      {2, 1024},
      {3, 2048},
      {4, 2048},
      {4, 4096},
      {5, 4096},
      {10, 8192},
      {15, 16384}};
    int t = 2;
    type = RLWE_Based;
    ofstream myfile;
    myfile.open("res_timing.txt");

    for (int i = 0; i < 8; i++) {
      int D = prms[i][0];
      int n = prms[i][1];
      std::cout << "D = " << D << " n = " << n << std::endl << "  ";
      params = fhe.Setup(0, D, type, to_ZZ(t), n);
      myfile << D << ", " << n << ", " << NumBits(params[0].q) << ", " << NumBits(params[D].q) << std::endl;
      std::cout << "####################" << std::endl;
    }
    myfile.close();
  }

  void TimingOperations() {
    int prms[][2] = {
      {1, 8},
      {1, 512},
      {2, 1024},
      {3, 2048},
      {4, 2048},
      {4, 4096},
      {5, 4096},
      {10, 8192},
      {15, 16384}};
    int t = 2;
    type = RLWE_Based;
    ofstream myfile;
    myfile.open("res_timing.txt");
    clock_t start;

    for (int i = 0; i < 1; i++) {
      int D = prms[i][0];
      int n = prms[i][1];
      std::cout << "D = " << D << " n = " << n << std::endl << "  ";
      params = fhe.Setup(0, D, type, to_ZZ(t), n);
      myfile << D << ", " << n << ", " << NumBits(params[0].q) << ", " << NumBits(params[D].q) << std::endl;

      std::cout << "Key gen .... " << std::endl;

      Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> sk_pk = fhe.Key_Gen(params);
      sk = sk_pk.first;
      pk = sk_pk.second;

      std::cout << "Encryption ... " << std::endl;
      
      std::vector<R_Ring_Number> messages;
      std::vector<FHE_Cipher_Text> ca, cm; // ciphers for add and ciphers for mul

      for (int j = 0; j < D + 1; j++) {
	R_Ring_Number m(to_ZZ(t), params[0].d);
	m = RandomBnd(ZZ(INIT_VAL, t));
	ca.push_back(fhe.Encrypt_for_Level(params, &pk, m, j == 0 ? D : D - j + 1));
	assert(ca[j].ThNoise != 0);
	ca[ca.size() - 1].Add_Secret_Key_Info(&sk);

	m = RandomBnd(ZZ(INIT_VAL, t));
	cm.push_back(fhe.Encrypt_for_Level(params, &pk, m, j == 0 ? D : D - j + 1));
	assert(cm[j].ThNoise != 0);
	cm[cm.size() - 1].Add_Secret_Key_Info(&sk);
      }

      for (int j = 1; j < D + 1; j++) {
	std::cout << "#level " << j << " ";
	start = clock();
	ca[0] = ca[0] + ca[j];
	std::cout << "Time for addition: " << (clock() - start) / (double)CLOCKS_PER_SEC << "s" << std::endl;
	start = clock();
	cm[0] = cm[0] * cm[j];
	std::cout << "Time for multiplication: " << (clock() - start) / (double)CLOCKS_PER_SEC << "s" << std::endl;
      }
      
      std::cout << "####################" << std::endl;
    }
    myfile.close();
  }
};

#endif /* _TIMING_H_ */
