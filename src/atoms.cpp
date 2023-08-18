
#include "atoms.h"

void 
Atoms::resize(int N_){
  N = N_;
  x.resize(N_);
  v.resize(N_);
  F.resize(N_);
  C.resize(N_);
  mole.resize(N_);
  noise.resize(N_);
}
 