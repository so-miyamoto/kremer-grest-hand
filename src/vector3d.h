/**
 * @file vector3d.h
 * @author Souta Miyamoto
 * @brief  simple 3d Vector type
 * @version 0.1
 * @copyright Copyright (c) 2020
 * 
 */

#pragma once

#include <array>
#include <iostream>
#include <cmath>

template <typename T>
using Vec = std::array<T,3>;

using Vecd = Vec<double>;
using Veci = Vec<int>;

// alias

// usage / test
/*
int main(){
  Vecd a;  a = Vecd{0.0,1.0,2.0};
  Vecd b(Vecd{1.0,2.0,3.0});

  a += b;
  std::cout<<"a    = "<<  a    <<std::endl;
  std::cout<<"b    = "<<  b    <<std::endl;
  std::cout<<"a+b  = "<< (a+b) <<std::endl;
  std::cout<<"a-b  = "<< (a-b) <<std::endl;
  std::cout<<"a*b  = "<< (a*b) <<std::endl;
  std::cout<<"a/b  = "<< (a/b) <<std::endl;
  std::cout<<"a<b  = "<< (a<b) <<std::endl;
  std::cout<<"a>b  = "<< (a>b) <<std::endl;
  std::cout<<"a==b = "<< (a==b)<<std::endl;
  std::cout<<"a!=b = "<< (a!=b)<<std::endl;

  return 0;
}
*/


// define operators +=, -=, *=, /=, +, -, *, /

#define VECTOR_OPS_ROLLING(op,op2)                                      \
  template <typename T>                                                 \
  inline Vec<T>& operator op (Vec<T>& left, const Vec<T> right)         \
  {                                                                     \
    left[0] op right[0];                                                \
    left[1] op right[1];                                                \
    left[2] op right[2];                                                \
    return left;                                                        \
  }                                                                     \
  template <typename T>                                                 \
  inline Vec<T> operator op2 (Vec<T> left, const Vec<T> right){         \
    return left op right;                                               \
  }                                                                     \
  template <typename T>                                                 \
  inline Vec<T>& operator op (Vec<T>& left, const T right)              \
  {                                                                     \
    left[0] op right;                                                   \
    left[1] op right;                                                   \
    left[2] op right;                                                   \
    return left;                                                        \
  }                                                                     \
  template <typename T>                                                 \
  inline Vec<T> operator op2 (Vec<T> left, const T right){              \
    return left op right;                                               \
  }                                                                     \
  template <typename T>                                                 \
  inline Vec<T> operator op2 (T left, Vec<T> right){                    \
    return right op left;                                               \
  }

VECTOR_OPS_ROLLING(+=,+); 
VECTOR_OPS_ROLLING(-=,-); 
VECTOR_OPS_ROLLING(*=,*); 
VECTOR_OPS_ROLLING(/=,/); 

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const Vec<T>& v){
  os << v[0] << " " << v[1] << " " << v[2];
  return os;
}

template<typename T>
inline std::istream& operator<<(std::istream& is, const Vec<T>& v){
  is >> v[0] >> v[1] >> v[2];
  return is;
}

template<typename T>
inline T 
norm2(const Vec<T>& v){return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; }

template<typename T>
inline T 
norm(const Vec<T>& v){return std::sqrt(norm2(v)); }

// }; // namespace Vector

