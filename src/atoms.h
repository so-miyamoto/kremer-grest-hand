/**
 * @file atoms.h
 * @author Souta Miyamoto
 * @brief 
 * @version 0.1
 * @date 2022-09-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include "vector3d.h"
#include <vector>  

struct Atoms{
  
  static constexpr double a = 1.0; // diameter
  static constexpr double m = 1.0; // mass
  std::vector<Vecd> x; // position 
  std::vector<Vecd> v; // velocity
  std::vector<Vecd> F; // force
  std::vector<int>  C; // component
  std::vector<int>  mole; // molecular
  std::vector<Vecd>  noise; // noise
  int N; // num of atoms

  void 
  resize(int N);
};