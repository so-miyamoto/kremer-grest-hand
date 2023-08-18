/**
 * @file rand.h
 * @author Souta Miyamoto (miyamoon98@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-10-27
 * 
 * @copyright Copyright (c) 2021
 * 
 * 
 * usage:
 * 
 * int main(){
 *   int seed = 10007;
 *   mrand::init( seed );
 *   for(int i=0; i<10; i++){
 *     std::cout<<mrand::rand()<<"\n"; // 0-1
 *     std::cout<<mrand::nrand()<<"\n"; // normal distribution N(0,1)
 *   }
 *   mrand::finalize();
 *   return 0;
 * }
 * 
 * 
 */
#include<cmath>
#include<random>

#ifdef _OPENMP
# include <omp.h>
#endif

#include"dSFMT.h"


namespace mrand{

  void init(const int seed);
  void finalize();
  void gen_nrand(double* array,const int size);
  double nrand();
  double rand();

}// namespace
