/**
 * @file md_system.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-09-16
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include "atoms.h"
#include "cell_list.h"
#include "pair_potential.h"
#include "sim_box.h"
#include "systemIO.h"
#include "vector3d.h"
#include "observer.h"
#include "bond.h"


struct MDsystem{

  Atoms atoms;
  SimBox simbox;
  CellList celllist;
  PairPotential pp;
  BondPotential bp;


  MDsystem();


  void
  adjust_positions(const double max_dr);


  void
  load_configure(SystemIO& sysio, const char* fname);

  void
  load_pair_potential(std::vector<std::vector<double>>& cutoff, std::vector<std::vector<double>>& LJeps, const double margine);


  inline bool
  finished_configuration() const {return _finish_configuration;}


private:
  int _configure;
  double _margine;
  double _margine_life;
  bool _finish_configuration;

};

