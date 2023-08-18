/**
 * @file pair_potential.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-09-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include<cmath>
#include "vector3d.h"

struct PairPotential {

  inline double 
  potential_LJ(const int Ci, const int Cj, const double r, const int n=6) const
  {
    double inv_rn = 1.0/std::pow(r,n);
    return 4.0*_LJeps[Ci][Cj]*( inv_rn*inv_rn - inv_rn );
  }

  inline double 
  potential_LJ_grad(const int Ci, const int Cj, const double r, const int n=6) const
  {
    double inv_r = 1.0/r;
    double inv_rn = 1.0/std::pow(r,n);
    return 4.0*_LJeps[Ci][Cj]*( - (2*n)*inv_rn* inv_rn*inv_r + n*inv_rn*inv_r );
  }


  inline double 
  potential_LJ_12_6(const int Ci, const int Cj, const double r) const
  {
    double r2 = r*r;
    double inv_r6 = 1.0/(r2*r2*r2);
    return 4.0*_LJeps[Ci][Cj]*( inv_r6*inv_r6 - inv_r6 );
  }

  inline double 
  potential_LJ_grad_12_6(const int Ci, const int Cj, const double r) const
  {
    double inv_r = 1.0/r;
    double inv_r2 = inv_r*inv_r;
    double inv_r6 = inv_r2*inv_r2*inv_r2;
    return _LJeps[Ci][Cj]*( - 48.0*inv_r6* inv_r6*inv_r + 24.0*inv_r6*inv_r );
  }

  void 
  set_LJparams(std::vector<std::vector<double>>& cutoff, std::vector<std::vector<double>>& LJeps)
  {
    _num_component = cutoff.size();
    _cutoff.resize(_num_component);
    _cutoff2.resize(_num_component);
    _LJeps.resize(_num_component);
    _LJtranc.resize(_num_component);
    for (int i = 0; i < _num_component; i++)
    {
      _cutoff[i].resize(_num_component);
      _cutoff2[i].resize(_num_component);
      _LJeps[i].resize(_num_component);
      _LJtranc[i].resize(_num_component);
    }    
    for (int i = 0; i < _num_component; i++){
      for (int j = 0; j < _num_component; j++){
        _cutoff[i][j]  = cutoff[i][j];
        _cutoff2[i][j] = cutoff[i][j]*cutoff[i][j];
        _LJeps[i][j] = LJeps[i][j];
        _LJtranc[i][j] = potential_LJ_12_6(i,j,_cutoff[i][j]);
    }}
    
  }

  std::vector<std::vector<double>>*
  get_cutoff()  
  {
    return &(_cutoff);
  }

  double 
  get_cutoff(const int Ci, const int Cj) const 
  {
    return _cutoff[Ci][Cj];
  }

  double
  get_LJtranc(const int Ci, const int Cj) const  
  {
    return _LJtranc[Ci][Cj];
  }

  bool 
  is_cutted_squared(const int Ci, const int Cj, double& r2) const 
  {
    return (r2 > _cutoff2[Ci][Cj]);
  }

  bool 
  is_cutted(const int Ci, const int Cj, const double& r) const 
  {
    return (r > _cutoff[Ci][Cj]);
  }

private:
  int _num_component;
  std::vector<std::vector<double>> _cutoff;
  std::vector<std::vector<double>> _cutoff2;
  std::vector<std::vector<double>> _LJeps;
  std::vector<std::vector<double>> _LJtranc;


};