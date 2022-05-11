#include "read_config.h"
#include "traj_reader.h"
#include "structures.h"
#include "welford.h"
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>

int main(int argc, char *argv[]) {

  // read parameters from json file
  std::string config_file = argv[1];
  std::cout << "Reading configuration file..." << std::endl;
  Read_config config(config_file);
  config.get_config();

  std::string datafile = config.j["datafile"];
  std::cout << "Reading data from " << datafile << std::endl;

  int num_atoms = config.j["atoms"];
  std::cout << num_atoms << " atoms" << std::endl;

  int num_frames = config.j["frames"];
  std::cout << num_frames << " frames" << std::endl;

  int num_skipframes = config.j["skipframes"];
  std::cout << "skipping " << num_skipframes << " frames" << std::endl;


  //lines below should all remove // to uncomment
  // create a resultSet obect to house our results
  resultSet<double> results;
  simFrame<double> frame1;
  std::vector<std::vector<size_t>> nbl;
  //NNeighbors<double>(filename, num_atoms, num_frames, num_nbs);
  // only results.avg is needed from variance00WK at thi point, that functionality
  // will be broken out into a more compact method at a later date
  
  
  clock_t t0, t1, t2;
  double t1sum=0.0;
  double t2sum=0.0;

  // start timing
  t0 = clock();
  /*
  // Initialization timing
  t1 = clock();
  t1sum = ((double)(t1-t0))/CLOCKS_PER_SEC;
  std::cout << "Init took " << t1sum << " seconds.  Begin compute" << std::endl;
  */

  variance00WK<double>(datafile, num_atoms,
  num_frames, num_skipframes, results);
  
  // serial timing
  t2 = clock();
  t2sum = ((double)(t2-t0))/CLOCKS_PER_SEC;
  std::cout << "Done. Compute took " << t2sum << " seconds" << std::endl;

  // values below not needed anymore but still being reported for reference
  std::cout << std::fixed << std::setprecision(16);
  // std::cout << "variance00= " << results.variance << std::endl;
  // std::cout << "std00= " << pow(results.variance,0.5) << std::endl;

  std::cout << "var00= " << results.variance << std::endl;
  std::cout << "varx00= " << results.var_xyz[0] << std::endl;
  std::cout << "vary00= " << results.var_xyz[1] << std::endl;
  std::cout << "varz00= " << results.var_xyz[2] << std::endl;
  std::cout << "std00= " << pow(results.variance,0.5) << std::endl;


  frame1 = results.avg;
  std::cout << "completed var00 calculation" << std::endl;
}
