#include "md_system.h"


MDsystem::MDsystem(){
  _configure = 0;
  _finish_configuration = false;
}


void
MDsystem::adjust_positions(const double max_dr){
  for (int i = 0; i < atoms.N; i++)
  {
    simbox.adjust_position(atoms.x[i]);
  }
  _margine_life -= max_dr;
  if( _margine_life < 0.0 ){
    celllist.register_addresses(atoms.x);
    _margine_life = _margine;
  }
}


void
MDsystem::load_configure(SystemIO& sysio, const char* fname){
  sysio.load_state(fname,atoms,simbox);
  celllist.set_boxsize(simbox.get_Lx(),simbox.get_Ly(),simbox.get_Lz());  
  if( _configure != 0 ){
    throw std::logic_error("MDsystem::invalid order of load_configure");
  }
  _configure++;
}

void
MDsystem::load_pair_potential(std::vector<std::vector<double>>& cutoff, std::vector<std::vector<double>>& LJeps, const double margine){
  pp.set_LJparams(cutoff,LJeps);
  double maximum_cutoff = 0.0;
  const int num_component = cutoff.size();
  for (int i = 0; i < num_component; i++)
  for (int j = 0; j < num_component; j++)
  {
    maximum_cutoff = std::max(maximum_cutoff,cutoff[i][j]);
  }
  
  celllist.cut_box(maximum_cutoff+margine);
  _margine = margine;
  _margine_life = margine;
  adjust_positions(margine+1.0);
  if( _configure != 1 ){
    throw std::logic_error("MDsystem::invalid order of set_pair_potential");
  }
  _configure++;
  _finish_configuration = true;
}

