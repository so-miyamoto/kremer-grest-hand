#include "systemIO.h"


SystemIO::SystemIO()
{ 
  _fp = nullptr;
}

SystemIO::~SystemIO()
{
  close_dumpfile();
}

void 
SystemIO::open_dumpfile(const char* fname)
{
  _fp = std::fopen(fname,"w");
  std::fprintf(stdout,"# opened dumpfile: %s\n",fname);
}

void 
SystemIO::close_dumpfile()
{
  if( _fp != nullptr ){
    std::fclose(_fp);
  }
  std::fprintf(stdout,"# closed dumpfile\n");
}

void 
SystemIO::dump(int step, Atoms &atoms, SimBox &simbox) const
{
  std::fprintf(_fp, "ITEM: TIMESTEP\n");
  std::fprintf(_fp, "%d\n", step);
  std::fprintf(_fp, "ITEM: NUMBER OF ATOMS\n");
  std::fprintf(_fp, "%d\n", atoms.N);
  std::fprintf(_fp, "ITEM: BOX BOUNDS xx yy zz\n");
  std::fprintf(_fp, "0.0 %f\n", simbox.get_Lx());
  std::fprintf(_fp, "0.0 %f\n", simbox.get_Ly());
  std::fprintf(_fp, "0.0 %f\n", simbox.get_Lz());
  std::fprintf(_fp, "ITEM: ATOMS id mol type x y z vx vy vz diameter\n");
  for (int i = 0; i < atoms.N; i++)
  {
    std::fprintf(_fp, "%d %d %d %f %f %f %f %f %f %f\n",
     i, atoms.mole[i], atoms.C[i], atoms.x[i][0], atoms.x[i][1], atoms.x[i][2], atoms.v[i][0], atoms.v[i][1], atoms.v[i][2], atoms.a);
  }
  return;
}

void 
SystemIO::load_state(const char* fname, Atoms &atoms, SimBox &simbox) const
{
  FILE *fp = std::fopen(fname, "r");
  char dummy[64];
  std::fprintf(stdout, "# reading %s\n", fname);

  int n, m, k, mole;
  double f;
  std::fgets(dummy, 64, fp);                      // std::fprintf(fp,"ITEM: TIMESTEP\n");
  std::fscanf(fp, "%d\n", &n);                    // std::fprintf(fp,"%d\n",step);
  std::fgets(dummy, 64, fp);                      // std::fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  std::fscanf(fp, "%d\n", &n);                    // std::fprintf(fp,"%lu\n",atoms.size());
  std::fgets(dummy, 64, fp);                      // std::fprintf(fp,"ITEM: BOX BOUNDS xx yy zz\n");
  Vecd dL;
  std::fscanf(fp, "%lf %lf\n", &f, &dL[0]); // std::fprintf(fp,"0.0 %f\n",simbox.L[0]);
  std::fscanf(fp, "%lf %lf\n", &f, &dL[1]); // std::fprintf(fp,"0.0 %f\n",simbox.L[1]);
  std::fscanf(fp, "%lf %lf\n", &f, &dL[2]); // std::fprintf(fp,"0.0 %f\n",simbox.L[2]);
  simbox.init_box(dL[0],dL[1],dL[2]);
  std::fgets(dummy, 64, fp);                      // std::fprintf(fp,"ITEM: ATOMS id type x y z vx vy vz diameter\n");

  atoms.resize(n);
  for (int i = 0; i < n; i++)
  {
    double qx, qy, qz;
    double px, py, pz;
    std::fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf %lf\n", &m, &mole, &k, &qx, &qy, &qz, &px, &py, &pz, &f);
    atoms.mole[i] = mole;
    atoms.C[i] = k;
    atoms.x[i] = Vecd{qx, qy, qz};
    atoms.v[i] = Vecd{px, py, pz};
  }
  std::fclose(fp);

  std::fprintf(stdout, "# atom size = %d\n", n);
  std::fprintf(stdout, "# box size = {%f, %f, %f}\n", dL[0], dL[1], dL[2]);
  std::fprintf(stdout, "# number density = %f\n", n / simbox.get_volume());
  std::fprintf(stdout, "# finish loading\n");
  return;
}


// void 
// SystemIO::load_pair_potential(const char* fname, PairPotential& pp){
//   FILE *fp = std::fopen(fname, "r");
//   int num_component;
//   std::vector<std::vector<double>> cutoff;
//   std::vector<std::vector<double>> LJeps;

//   std::fprintf(stdout, "# reading %s\n", fname);

//   std::fscanf(fp,"%d\n",&num_component);
//   std::fprintf(stdout, "# num_component = %d\n",num_component);

//   cutoff.resize(num_component);
//   LJeps.resize(num_component);
//   for (int i = 0; i < num_component; i++)
//   {
//     cutoff[i].resize(num_component);
//     LJeps[i].resize(num_component);
//   }

//   std::fprintf(stdout, "# cutoff = \n");
//   for (int i = 0; i < num_component; i++){
//     for (int j = 0; j < num_component; j++){
//       std::fscanf(fp,"%lf",&cutoff[i][j]);
//       std::fprintf(stdout, "%f ",cutoff[i][j]);
//     }
//     std::fscanf(fp,"\n");
//     std::fprintf(stdout, "\n");
//   }

//   std::fprintf(stdout, "# LJepsilon = \n");
//   for (int i = 0; i < num_component; i++){
//     for (int j = 0; j < num_component; j++){
//       std::fscanf(fp,"%lf",&LJeps[i][j]);
//       std::fprintf(stdout, "%f ",LJeps[i][j]);
//     }
//     std::fscanf(fp,"\n");
//     std::fprintf(stdout, "\n");
//   }
//   std::fclose(fp);

//   pp.set_LJparams(cutoff,LJeps);
// };

