/**
 * @file cell_list.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-09-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include "vector3d.h"

struct CellList{

  CellList();

  void 
  set_boxsize(double Lx,double Ly, double Lz);

  void 
  cut_box(double minimum_cell_length);

  void
  register_addresses(std::vector<Vecd>& pos);

  std::vector<int>*
  get_neighbors(const Vecd p);


private:
  bool _initialized_boxsize;
  bool _cutted_box;
  double _Lx, _Ly, _Lz;
  double _cellLx, _cellLy, _cellLz;
  int _num_cells, 
      _num_cells_x, 
      _num_cells_y, 
      _num_cells_z;
  std::vector<std::vector<int>> _cells_list;

  int 
  _idx(int ix, int iy, int iz) const ;

};
