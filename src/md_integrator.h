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

#pragma once
#include "md_system.h"
#include "atoms.h"
#include "sim_box.h"
#include "cell_list.h"
#include "pair_potential.h"
#include "bond.h"
#include "rand.h"


struct Integrator{

  inline void
  set_params(const int max_step_, const int obs_step_, const double dt_)
  {
    max_step = max_step_;
    obs_step = obs_step_;
    dt = dt_;
  }

  inline void
  set_params_for_temperature(
    const int scaling_steps_, 
    const double temperature_,
    const int final_scaling_step_,
    const bool use_langevin_,
    const double langevin_friction_,
    const bool use_nose_hoover_,
    const double tau_nose_           )
  {
    velocity_scaling = ( scaling_steps_ > 0 );
    scaling_steps = scaling_steps_;
    temperature = temperature_;
    final_scaling_step = final_scaling_step_;

    use_langevin = use_langevin_;
    langevin_friction = langevin_friction_;
    langevin_b = 1.0 / ( 1.0 + langevin_friction*0.5*dt );
    langevin_a = ( 1.0 - langevin_friction*0.5*dt ) * langevin_b;

    zeta_nose = 0.0;
    tau_nose = tau_nose_;
    use_nose_hoover = use_nose_hoover_;
  };

  int max_step;
  int obs_step;
  double dt;

  bool velocity_scaling;
  int scaling_steps;
  int final_scaling_step;
  double temperature;

  bool use_langevin;
  double langevin_friction;
  double langevin_a;
  double langevin_b;

  bool use_nose_hoover;
  double zeta_nose;
  double tau_nose;
};

struct Simulator{

  inline void
  run(MDsystem& mdsys, Integrator& inte, Observer& obs, SystemIO& sysio){
    if( inte.use_langevin ){
      _langevin_run(mdsys,inte,obs,sysio);
    }else{
      _normal_run(mdsys,inte,obs,sysio);
    }
  };


private:

  void 
  _normal_run(MDsystem& mdsys, Integrator& inte, Observer& obs, SystemIO& sysio);

  void 
  _langevin_run(MDsystem& mdsys, Integrator& inte, Observer& obs, SystemIO& sysio);

  void
  _update_velocities(Atoms& atoms, const double dt);

  // void
  // _update_velocities_langevin(Atoms& atoms, const double dt, double temperature, double gamma);

  void
  _update_velocities_nose_hoover(Atoms& atoms, const double dt, double temperature, double now_temperature, double& zeta_nose, const double tau_nose);

  void
  _update_velocities_langevin_1(Atoms& atoms, const double dt, const double langevin_a, const double langevin_b);

  void
  _update_velocities_langevin_2(Atoms& atoms, const double dt);

  double
  _update_positions_half(Atoms& atoms, const double hdt);

  double
  _update_positions_langevin(Atoms& atoms, const double dt ,const double langevin_b);

  void
  _compute_force(Atoms& atoms, SimBox& simbox, CellList& celllist,PairPotential& pp, BondPotential& bp);

  void
  _gen_thermal_noise_langevin(Atoms& atoms, const double dt,const double temperature, const double friction);

  void 
  _scaling_velocity(Atoms& atoms, const double alpha);

};