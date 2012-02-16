#include "FHE.h"
#include <iostream>

int main (int argc, char * const argv[]) {
  int L = 2; // number of multiplication levels in the circuit + 1
  ZZ modul = ZZ(INIT_VAL, 31); // boundary on numbers involved
  FHE fhe;
  FHE_Params params = fhe.Setup(3, L, RLWE_Based, modul);
  Pair<FHE_Secret_Key_Type, FHE_Public_Key_Type> sk_pk = fhe.Key_Gen(params);

  ZZ m1 = ZZ(INIT_VAL, 5), m2 = ZZ(INIT_VAL, -3);
  std::cout << m1 << " * " << m2 << " = ";
  
  FHE_Cipher_Text c1 = fhe.Encrypt(params, &sk_pk.second, m1);
  FHE_Cipher_Text c2 = fhe.Encrypt(params, &sk_pk.second, m2);

  FHE_Cipher_Text sum = c1 * c2;
  ZZ result = (ZZ)sum.Decrypt(params, sk_pk.first);

  std::cout << result << std::endl;
  
  return 0;
}
