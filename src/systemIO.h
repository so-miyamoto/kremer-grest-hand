#pragma once

#include "atoms.h"
#include "sim_box.h"

#include "pair_potential.h"

#include <cstdio>
#include <string>

struct SystemIO
{

  SystemIO();
  ~SystemIO();

  void 
  open_dumpfile(const char* fname);

  void 
  close_dumpfile();

  void 
  dump(int step, Atoms &atoms, SimBox &simbox) const;

  void 
  load_state(const char* fname, Atoms &atoms, SimBox &simbox) const;

  // void 
  // load_pair_potential(const char* fname, PairPotential& pp);

  private:
    FILE* _fp;

};
