/**
 * @file rand.cpp
 * @author Souta Miyamoto (miyamoon98@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-10-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include"rand.h"

namespace mrand{

  namespace{
    const int arysize = (2*DSFMT_N);
#ifdef _OPENMP
    dsfmt_t* dsfmt;    
    double* rand_buffer;
    double* nrand_buffer;
    int* rest_nran_buf;
    int* rest_ran_buf;
#else
    dsfmt_t dsfmt;    
    double  rand_buffer[arysize];
    double nrand_buffer[arysize];
    int rest_nran_buf;
    int rest_ran_buf;
#endif
  }

  void init(const int seed){
#ifdef _OPENMP
    std::mt19937 mt(seed);
    int num_threads = omp_get_max_threads();
    dsfmt = new dsfmt_t [num_threads];    
    rand_buffer = new double [num_threads*arysize];
    nrand_buffer = new double [num_threads*arysize];
    rest_nran_buf = new int [num_threads];
    rest_ran_buf = new int [num_threads];
    for (int i = 0; i < num_threads; i++)
    {
      dsfmt_init_gen_rand(&dsfmt[i],mt());
      rest_ran_buf[i] = 0;
      rest_nran_buf[i] = 0;
    }    
#else
    dsfmt_init_gen_rand(&dsfmt,seed);
    rest_ran_buf = 0;
    rest_nran_buf = 0;
#endif
  }

  void finalize()
  {
#ifdef _OPENMP
    delete[] dsfmt;
    delete[] rand_buffer;
    delete[] nrand_buffer;
    delete[] rest_nran_buf;
    delete[] rest_ran_buf;
#endif
  }

  void gen_rand(double* array,const int size)
  {
    assert(size%2==0);
#ifdef _OPENMP
    dsfmt_fill_array_close_open(&dsfmt[omp_get_thread_num()],array,size);
#else
    dsfmt_fill_array_close_open(&dsfmt,array,size);
#endif
  }

  void gen_nrand(double* array,const int size)
  {
    assert(size%2==0);
#ifdef _OPENMP
    dsfmt_fill_array_close_open(&dsfmt[omp_get_thread_num()],array,size);
#else
    dsfmt_fill_array_close_open(&dsfmt,array,size);
#endif
    constexpr double TWO_PI = 2.0*M_PI;
    for( int i=0; i<size; i+=2 ){
      const double sqrt_2log_r0 =  std::sqrt( std::fabs( 2.0 * std::log(array[i]) ));
      const double TWO_PI_r1    =  TWO_PI*array[i+1];
      array[i  ] = sqrt_2log_r0 * std::cos( TWO_PI_r1 );
      array[i+1] = sqrt_2log_r0 * std::sin( TWO_PI_r1 );
    }
  }

  double rand(){
#ifdef _OPENMP
    int omp_i = omp_get_thread_num();
    if( rest_ran_buf[omp_i]==0 ){
      gen_rand(&rand_buffer[omp_i*arysize],arysize);
      rest_ran_buf[omp_i] = arysize;
    }
    return rand_buffer[omp_i*arysize+(--rest_ran_buf[omp_i])];
#else
    if( rest_ran_buf==0 ){
      gen_rand(rand_buffer,arysize);
      rest_ran_buf = arysize;
    }
    return rand_buffer[--rest_ran_buf];
#endif
  }

  double nrand(){
#ifdef _OPENMP
    int omp_i = omp_get_thread_num();
    if( rest_nran_buf[omp_i]==0 ){
      gen_nrand(&nrand_buffer[omp_i*arysize],arysize);
      rest_nran_buf[omp_i] = arysize;
    }
    return nrand_buffer[omp_i*arysize+(--rest_nran_buf[omp_i])];
#else
    if( rest_nran_buf==0 ){
      gen_nrand(nrand_buffer,arysize);
      rest_nran_buf = arysize;
    }
    return nrand_buffer[--rest_nran_buf];
#endif
  }

}// namespace rand