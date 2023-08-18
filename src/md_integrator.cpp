/**
 * @file compute_force.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-09-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include "md_integrator.h"


void 
Simulator::_normal_run(MDsystem& mdsys, Integrator& inte, Observer& obs, SystemIO& sysio)
{
  if( !mdsys.finished_configuration() )
    throw std::logic_error("MDsystem before initialized");

  std::fprintf(stdout,"# activate standard integrator\n");
  
  for (int step = 0; step < inte.max_step; step++)
  {
    if(inte.velocity_scaling){
      bool do_scaling = (step % inte.scaling_steps == 0) and (step < inte.final_scaling_step);
      if( do_scaling ){
        double alpha = std::sqrt(inte.temperature/obs.temperature(mdsys.atoms));
        _scaling_velocity(mdsys.atoms,alpha);
        std::printf("# step=%d scaling velocity with %f\n",step,alpha);
      }
    }
    if( step % inte.obs_step == 0 ){
      std::printf("# step=%d\n",step);
      obs.watch(step,mdsys.atoms,mdsys.simbox,mdsys.celllist,mdsys.pp,mdsys.bp);
      sysio.dump(step,mdsys.atoms,mdsys.simbox);
    }

    // update position t -> t + dt/2
    mdsys.adjust_positions(2.0*_update_positions_half(mdsys.atoms,0.5*inte.dt));

    _compute_force(mdsys.atoms,mdsys.simbox,mdsys.celllist,mdsys.pp, mdsys.bp);

    // update velocity t -> t + dt
    if( inte.use_nose_hoover ){
      _update_velocities_nose_hoover(mdsys.atoms,inte.dt,inte.temperature,obs.temperature(mdsys.atoms),inte.zeta_nose,inte.tau_nose);
    }//  _update_velocities(mdsys.atoms,inte.dt);
    _update_velocities(mdsys.atoms,inte.dt);

    // update position t + dt/2 -> t + dt
    mdsys.adjust_positions(2.0*_update_positions_half(mdsys.atoms,0.5*inte.dt));

  }
  
}

void 
Simulator::_langevin_run(MDsystem& mdsys, Integrator& inte, Observer& obs, SystemIO& sysio)
{
  if( !mdsys.finished_configuration() )
    throw std::logic_error("MDsystem before initialized");

  std::fprintf(stdout,"# activate langevin integrator\n");

  _compute_force(mdsys.atoms,mdsys.simbox,mdsys.celllist,mdsys.pp, mdsys.bp);

  for (int step = 0; step < inte.max_step; step++)
  {
    if(inte.velocity_scaling){
      bool do_scaling = (step % inte.scaling_steps == 0) and (step < inte.final_scaling_step);
      if( do_scaling ){
        double alpha = std::sqrt(inte.temperature/obs.temperature(mdsys.atoms));
        _scaling_velocity(mdsys.atoms,alpha);
        std::printf("# step=%d scaling velocity with %f\n",step,alpha);
      }
    }
    if( step % inte.obs_step == 0 ){
      std::printf("# step=%d\n",step);
      obs.watch(step,mdsys.atoms,mdsys.simbox,mdsys.celllist,mdsys.pp,mdsys.bp);
      sysio.dump(step,mdsys.atoms,mdsys.simbox);
    }


    // compute noise beta^n+1
    _gen_thermal_noise_langevin(mdsys.atoms,inte.dt,inte.temperature,inte.langevin_friction);

    // update positions t -> t + dt
    mdsys.adjust_positions(2.0*_update_positions_langevin(mdsys.atoms,inte.dt,inte.langevin_b));

    // kick velocity by F^n
    _update_velocities_langevin_1(mdsys.atoms,inte.dt,inte.langevin_a,inte.langevin_b);

    // compute force F^n+1
    _compute_force(mdsys.atoms,mdsys.simbox,mdsys.celllist,mdsys.pp, mdsys.bp);

    // kick velocity by F^n+1
    _update_velocities_langevin_2(mdsys.atoms,inte.dt);

  }
  
}

// void
// Simulator::_update_velocities_langevin(Atoms& atoms, const double dt, double temperature, double gamma){
//   static std::mt19937 mt(42);
//   double T = temperature;
//   double D = std::sqrt(2.0*gamma*T/dt);
//   std::normal_distribution<double> nd(0.0,D);
//   for (int i = 0; i < atoms.N; i++)
//   {
//     atoms.v[i][0] += (atoms.F[i][0] - gamma*atoms.m*atoms.v[i][0] + nd(mt) )*(dt/atoms.m);
//     atoms.v[i][1] += (atoms.F[i][1] - gamma*atoms.m*atoms.v[i][1] + nd(mt) )*(dt/atoms.m);
//     atoms.v[i][2] += (atoms.F[i][2] - gamma*atoms.m*atoms.v[i][2] + nd(mt) )*(dt/atoms.m);
//   }
//   return ;
// }

void
Simulator::_update_velocities_nose_hoover(Atoms& atoms, const double dt, double temperature, double now_temperature, double& zeta_nose, const double tau_nose){
  double dT = now_temperature - temperature;
  zeta_nose += dT / (tau_nose*tau_nose) * dt;
  for (int i = 0; i < atoms.N; i++)
  {
    atoms.v[i] -= atoms.v[i] * zeta_nose*dt;
  }
  return ;
}

void
Simulator::_update_velocities(Atoms& atoms, const double dt){
  for (int i = 0; i < atoms.N; i++)
  {
    atoms.v[i] += atoms.F[i]*(dt/atoms.m);
  }
  return ;
}

void
Simulator::_update_velocities_langevin_1(Atoms& atoms, const double dt, const double langevin_a, const double langevin_b){
  for (int i = 0; i < atoms.N; i++)
  {
    atoms.v[i] = ((atoms.v[i] + atoms.F[i]*(0.5*dt/atoms.m)) * langevin_a) + (atoms.noise[i]*(langevin_b/atoms.m));
  }
  return ;
}

void
Simulator::_update_velocities_langevin_2(Atoms& atoms, const double dt){
  for (int i = 0; i < atoms.N; i++)
  {
    atoms.v[i] += (atoms.F[i]*(0.5*dt/atoms.m));
  }
  return ;
}

double
Simulator::_update_positions_half(Atoms& atoms, const double hdt){
  double max_move2 = 0.0;
  for (int i = 0; i < atoms.N; i++)
  {
    auto dx = atoms.v[i]*hdt;
    max_move2 = std::max(max_move2,norm2(dx));
    atoms.x[i] += dx;
  }
  return std::sqrt(max_move2);
}

double
Simulator::_update_positions_langevin(Atoms& atoms, const double dt,const double langevin_b){
  double max_move2 = 0.0;
  for (int i = 0; i < atoms.N; i++)
  {
    auto dx = ( (atoms.v[i]*dt) + (atoms.F[i]*(dt*dt*0.5/atoms.m)) + (atoms.noise[i]*(0.5*dt/atoms.m)));
    dx *= langevin_b;
    max_move2 = std::max(max_move2,norm2(dx));
    atoms.x[i] += dx;
  }
  return std::sqrt(max_move2);
}

void
Simulator::_gen_thermal_noise_langevin(Atoms& atoms, const double dt,const double temperature, const double friction){
  const double coeff = std::sqrt(2.0*temperature*friction*dt);

  mrand::gen_nrand(&atoms.noise[0][0],3*atoms.N);

  for(int i = 0; i < atoms.N; i++)
  {
    atoms.noise[i] *= coeff;
  }
  return ;
}


void
Simulator::_compute_force(Atoms& atoms, SimBox& simbox, CellList& celllist,PairPotential& pp, BondPotential& bp)
{
  const int N = atoms.N;
  for (int i = 0; i < N; i++)
    atoms.F[i] = Vecd{0.0,0.0,0.0};

  // bond interaction 
  
  for (int ai = 0; ai < N-1; ai++)
  {
    int i = ai;
    int j = ai+1;
    if( atoms.mole[i] != atoms.mole[j] ) continue;
    Vecd xij = atoms.x[i] - atoms.x[j]; 
    simbox.adjust_relative_vector(xij);
    double r = norm(xij);
    const double mdUdr_r = - bp.bond_potential_grad(r)/r;
    atoms.F[i] += (xij*mdUdr_r);      
    atoms.F[j] -= (xij*mdUdr_r);      
  }
  
  
  // pair interaction 
  for (int i = 0; i < N; i++)
  {
    auto& xi = atoms.x[i];
    const int Ci = atoms.C[i];
    std::vector<int>* jlist = celllist.get_neighbors(xi);
    const int Nj = jlist->size();

    for (int nj = 0; nj < Nj; nj++)
    {
      const int j = (*jlist)[nj];
      if( !(i < j) ) continue;
      // if( atoms.mole[i] == atoms.mole[j] ) continue;
      auto& xj = atoms.x[j];
      const int Cj = atoms.C[j];
      Vecd xij = xi - xj;
      simbox.adjust_relative_vector(xij);
      double r = norm(xij);
      if( pp.is_cutted(Ci,Cj,r) ) continue;     
      double mdUdr_r = - pp.potential_LJ_grad_12_6(Ci,Cj,r)/r;
      atoms.F[i] += (xij*mdUdr_r);      
      atoms.F[j] -= (xij*mdUdr_r);      
    }
  }
  // for (int i = 0; i < N; i++){
  //   std::cout<<atoms.F[i]<<std::endl;
  // }

  return;
}

void 
Simulator::_scaling_velocity(Atoms& atoms, const double alpha){
  for (int i = 0; i < atoms.N; i++)
    atoms.v[i] *= alpha;
  return;
};
