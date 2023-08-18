#include "md_system.h"
#include "observer.h"
#include "systemIO.h"
#include "md_integrator.h"
#include "timer.h"
#include "rand.h"
#include "yaml-cpp/yaml.h"

#include <iostream>
#include <fstream>
#include <random>

int main(int argc, char** argv)
{
  try{
  
  if( argc != (1+5) ) throw std::runtime_error("input 5 arguments this.exe [params.yml] [init.dump] [run.dump] [final.dump] [logfile.log]");
  auto param_fname = argv[1];
  auto init_fname  = argv[2];
  auto run_fname   = argv[3];
  auto final_fname = argv[4];
  auto log_fname   = argv[5];

  auto yamlf = YAML::LoadFile(param_fname);

  std::fprintf(stdout,"# read potentials configure from yaml\n");

  YAML::Node&& potentialf = yamlf["potential"];
  int    bond_type     = potentialf["bond"]["type"].as<int>();
  double k_hookean     = potentialf["bond"]["spring_constant"].as<double>(); 
  double equi_length   = potentialf["bond"]["equi_length"].as<double>();
  double max_length    = potentialf["bond"]["max_length"].as<double>();

  double margine       = potentialf["pair"]["margine"].as<double>(); 
  int    types         = potentialf["pair"]["num_types"].as<int>();
  std::vector<std::vector<double>> cutoff, LJeps;
  for (int i = 1; i <= types; i++)
  {
    std::string num = std::to_string(i);
    cutoff.emplace_back( potentialf["pair"]["cutoff"][num].as<std::vector<double>>() ); 
    LJeps .emplace_back( potentialf["pair"]["LJeps"][num] .as<std::vector<double>>() );
  }

  std::fprintf(stdout,"# read integrate configure from yaml\n");

  int max_step             = yamlf["integrate"]["max_step"].as<int>();  
  int obs_step             = yamlf["integrate"]["obs_step"].as<int>(); 
  double dt                = yamlf["integrate"]["dt"]      .as<double>();

  std::fprintf(stdout,"# read thermostat configure from yaml\n");

  double temperature       = yamlf["thermostat"]["temperature"].as<double>(); 
  int scaleT               = yamlf["thermostat"]["scaling_interval"].as<int>(); 
  int max_scaling_step     = yamlf["thermostat"]["scaling_final_step"].as<int>();
  bool use_langevin        = yamlf["thermostat"]["use_langevin"].as<bool>();  
  double langevin_friction = yamlf["thermostat"]["langevin_friction"].as<double>();  
  int use_nose_hoover      = yamlf["thermostat"]["use_nose_hoover"].as<bool>();
  int nose_tau             = yamlf["thermostat"]["nose_tau"].as<double>();

  int seed             = yamlf["option"]["seed"].as<int>();
  bool use_random_seed = yamlf["option"]["use_random_seed"].as<bool>();
  if( use_random_seed ){
    std::random_device rg;
    seed = rg();
  }


  std::fprintf(stdout,"# set main objects for simulation\n");

  ChronoTimer timer;
  timer.start();

  mrand::init(seed);

  SystemIO sysio;
  sysio.open_dumpfile(run_fname);

  MDsystem mdsys;
  mdsys.load_configure(sysio,init_fname);
  mdsys.load_pair_potential(cutoff,LJeps,margine);
  mdsys.bp.set_bond_params(bond_type,k_hookean,equi_length,max_length);

  Observer obs;
  obs.open_logfile(log_fname);

  Integrator inte;
  inte.set_params(max_step,obs_step,dt);
  inte.set_params_for_temperature(
    scaleT,temperature,max_scaling_step,
    use_langevin,langevin_friction,
    use_nose_hoover,nose_tau);

  Simulator simulator;
  simulator.run(mdsys,inte,obs,sysio);

  std::cout<<"# finish simulator run"<<std::endl;

  sysio.close_dumpfile();

  sysio.open_dumpfile(final_fname);
  sysio.dump(max_step,mdsys.atoms,mdsys.simbox);

  timer.end();
  std::cout<<timer.elapsed_sec()<<" [sec]"<<std::endl;
  std::cout<<timer.elapsed_msec()<<" [msec]"<<std::endl;
  timer.stamp_time();

  mrand::finalize();

  std::cout<<"# === end of simulation :) ==="<<std::endl;

  }catch(YAML::Exception& e){
    std::cout<<e.what()<<std::endl;
  }catch(std::exception& e){
    std::cout<<e.what()<<std::endl;
  }
  return 0;
}
