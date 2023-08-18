#pragma once

#include "atoms.h"
#include "sim_box.h"
#include "cell_list.h"
#include "pair_potential.h"
#include "bond.h"
#include <cstdio>

struct Observer
{

private:
  FILE *_fp;

public:
  Observer();
  ~Observer();

  void 
  open_logfile(const char* log_fname);

  void 
  watch(const int step, Atoms &atoms, SimBox &simbox, CellList& celllist ,PairPotential &pp, BondPotential& bp);

  Vecd 
  average_momentum(Atoms &atoms);

  std::vector<double> 
  potential_energy(Atoms &atoms, SimBox &simbox, CellList& celllist ,PairPotential &pp, BondPotential& bp);


  double 
  kinetic_energy(Atoms &atoms);

  inline double 
  temperature(Atoms &atoms)
  {
    return kinetic_energy(atoms) / 1.5;
  }
};

