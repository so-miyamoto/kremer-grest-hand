#include "observer.h"



Observer::Observer()
{
  _fp = nullptr;
}

Observer::~Observer()
{
  if(_fp!=nullptr){
    std::fclose(_fp);
    std::fprintf(stdout, "# observer close a file: stat.log\n");
  }
}
void 
Observer::open_logfile(const char* log_fname){
  std::fprintf(stdout, "# observer initialized\n");
  _fp = std::fopen(log_fname, "w");
  std::fprintf(_fp, "# step V_all V_bond V_pair  K  Eall  T  P  px  py  pz\n");
  std::fprintf(stdout, "# observer open a new file: %s\n",log_fname);    
}

Vecd 
Observer::average_momentum(Atoms &atoms)
{
  Vecd pm = {0.0, 0.0, 0.0};
  for (auto &a : atoms.v)
  {
    pm[0] += a[0];
    pm[1] += a[1];
    pm[2] += a[2];
  }
  pm[0] /= (atoms.m*atoms.N);
  pm[1] /= (atoms.m*atoms.N);
  pm[2] /= (atoms.m*atoms.N);
  return pm;
}

double 
Observer::kinetic_energy(Atoms &atoms)
{
  double K = 0.0;
  for (auto &a : atoms.v)
  {
    K += norm2(a);
  }
  return 0.5 * atoms.m * K / atoms.N;
}

void 
Observer::watch(const int step, Atoms &atoms, SimBox &simbox, CellList& celllist ,PairPotential &pp, BondPotential& bp)
{
  std::vector<double> obs;
  std::vector<double> V_P_T = potential_energy(atoms, simbox, celllist, pp, bp);
  const double V_all = V_P_T[0];
  const double V_bond = V_P_T[1];
  const double V_pair = V_P_T[2];
  const double K = V_P_T[4] * 1.5;
  const double P = V_P_T[3];
  const double T = V_P_T[4];
  Vecd mom = average_momentum(atoms);
  std::fprintf(_fp, "%d %f %f %f %f %f %f %f %f %f %f\n", step, V_all, V_bond, V_pair, K, V_all + K, T, P, mom[0], mom[1], mom[2]);
  return;
}

std::vector<double> 
Observer::potential_energy(Atoms &atoms, SimBox &simbox, CellList& celllist ,PairPotential &pp, BondPotential& bp)
{
  
  const int N = atoms.N;
  const double Vol = simbox.get_volume();
  const double density = N / Vol;
  const double T = temperature(atoms);

  double V_all = 0.0;
  double V_pair = 0.0;
  double V_bond = 0.0;
  double virial = 0.0;

  // bond interaction
  for (int ai = 0; ai < N-1; ai++)
  {
    int i = ai;
    int j = ai+1;
    if( atoms.mole[i] != atoms.mole[j] ) continue;
    Vecd xij = atoms.x[i] - atoms.x[j]; 
    simbox.adjust_relative_vector(xij);
    double r = norm(xij);
    V_bond += bp.bond_potential(r);
    virial += - r * bp.bond_potential_grad(r);
  }


  // pair interaction
  for (int i = 0; i < N; i++)
  {
    auto& xi = atoms.x[i];
    const int Ci = atoms.C[i];
    std::vector<int>* jlist = celllist.get_neighbors(xi);
    const int Nj = jlist->size();
    for (int nj=0; nj<Nj; nj++)
    {
      const int j = (*jlist)[nj];
      if (!(i < j)) continue;
      // if( atoms.mole[i] == atoms.mole[j] ) continue;
      auto& xj = atoms.x[j];
      const int Cj = atoms.C[j];
      Vecd xij = xi - xj;
      simbox.adjust_relative_vector(xij);
      double r2 = norm2(xij);
      if ( pp.is_cutted_squared(Ci,Cj,r2) ) continue;
      double r = std::sqrt(r2);
      V_pair += pp.potential_LJ_12_6(Ci,Cj,r) - pp.get_LJtranc(Ci,Cj);
      virial += - r * pp.potential_LJ_grad_12_6(Ci,Cj,r);
    }
  }
  virial /= (3.0 * Vol);
  V_all = V_bond + V_pair;
  
  return std::vector<double>{V_all/N, V_bond/N, V_pair/N, density * T + virial, T}; // return V, P, T
}
