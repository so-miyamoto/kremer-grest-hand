/**
 * @file bond.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-09-16
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#define LENGTH_RATE_TO_MAXIMUM 0.99

struct BondPotential{

  BondPotential(){
    set_bond_params(0,0.0,0.0,100000.0);
  }

  void
  set_bond_params(const int type, const double spring_const, const double equi_length, const double max_length){
    _type = type;
    _spring_const = spring_const;
    _equi_length = equi_length;
    _max_length = max_length;
  };


  inline double
  bond_potential(const double r) const {
    if(_type==0){
      return hookean_bond_potential(r);
    }else if(_type==1){
      return fene_bond_potential(r);
    }else{
      throw std::runtime_error("unknown bond!");
    }
  };

  inline double
  bond_potential_grad(const double r) const {
    if(_type==0){
      return hookean_bond_potential_grad(r);
    }else if(_type==1){
      return fene_bond_potential_grad(r);
    }else{
      throw std::runtime_error("unknown bond!");
    }
  };

  double 
  hookean_bond_potential(const double r) const {
    const double dr = r-_equi_length;
    return 0.5*_spring_const*dr*dr;
  };

  double 
  hookean_bond_potential_grad(const double r) const {
    const double dr = r-_equi_length;
    return _spring_const*dr;
  };

  double 
  fene_bond_potential(const double r) const {
    const double dr = std::min(r/_max_length, LENGTH_RATE_TO_MAXIMUM);
    return -0.5*_spring_const*_max_length*_max_length*std::log(1.0-dr*dr);
  };

  double 
  fene_bond_potential_grad(const double r) const {
    const double dr = std::min(r/_max_length, LENGTH_RATE_TO_MAXIMUM);
    return 0.5*_spring_const*r/(1-dr*dr);
  };

private:
  int _type; // 0:hookean 1:fene
  double _spring_const;
  double _equi_length;
  double _max_length;

};