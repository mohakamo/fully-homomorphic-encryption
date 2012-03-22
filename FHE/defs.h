#ifndef _DEFS_H_
#define _DEFS_H_
#include "NTL/ZZ.h"

NTL_CLIENT

extern int Ring_Number_d;

class R_Ring_Vector;

void Bit_Decomposition(const R_Ring_Vector &x, ZZ q, R_Ring_Vector &result);
void Powersof2(const R_Ring_Vector &x, ZZ q, R_Ring_Vector &result);

#endif /* _DEFS_H_ */
