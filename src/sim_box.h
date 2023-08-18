/**
 * @file sim_box.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-09-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include"vector3d.h"

struct SimBox{
  
  inline void
  init_box(double Lx, double Ly, double Lz)
  {
    _Lx = Lx;
    _Ly = Ly;
    _Lz = Lz;
    _hLx = 0.5*Lx;
    _hLy = 0.5*Ly;
    _hLz = 0.5*Lz;
    _volume = Lx*Ly*Lz;
  }

  inline double get_Lx() const {return _Lx;}
  inline double get_Ly() const {return _Ly;}
  inline double get_Lz() const {return _Lz;}
  inline double get_volume() const {return _volume;}

  inline void
  adjust_relative_vector(Vecd& q) const 
  {
    constexpr int x=0,y=1,z=2;
    if (q[x] <  -_hLx) q[x] += _Lx; // left x
    if (q[y] <  -_hLy) q[y] += _Ly; // left y
    if (q[z] <  -_hLz) q[z] += _Lz; // left z
    if (q[x] >= +_hLx) q[x] -= _Lx; // right x
    if (q[y] >= +_hLy) q[y] -= _Ly; // right y
    if (q[z] >= +_hLz) q[z] -= _Lz; // right z  
  }

  inline void
  adjust_position(Vecd& q) const 
  {
    constexpr int x=0,y=1,z=2;
    if (q[x] <  0.0) q[x] += _Lx; // left x
    if (q[y] <  0.0) q[y] += _Ly; // left y
    if (q[z] <  0.0) q[z] += _Lz; // left z
    if (q[x] >= _Lx) q[x] -= _Lx; // right x
    if (q[y] >= _Ly) q[y] -= _Ly; // right y
    if (q[z] >= _Lz) q[z] -= _Lz; // right z  
  }



private:

  double _Lx,_Ly,_Lz;
  double _hLx,_hLy,_hLz;
  double _volume;


};