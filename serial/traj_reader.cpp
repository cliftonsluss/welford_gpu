#include "welford.h"
#include "traj_reader.h"

Trajectory::Trajectory(std::string &filename, const size_t num_atoms,
    const size_t header) {
  Trajectory::header = header;
  Trajectory::num_atoms = num_atoms;
  Trajectory::inputfile.open(filename, std::ifstream::in);
}

Trajectory::Trajectory(std::string &filename) {
  Trajectory::outputfile.open(filename, std::ofstream::out);
  Trajectory::frame_num = 0;
}

void Trajectory::getNextFrame(simFrame<double> &frame) {
  std::string temp;
  double x, y, z, xc, yc, zc;
  size_t n, nt;
  for (size_t i = 0; i < Trajectory::header; i++) {
    getline(inputfile, temp);
    // std::cout << temp << std::endl;
  }

  // read min max in all dimensions, then redetermine min max to center
  // simulation box around the origin. This is to eliminate overrun
  // errors for atoms close to 0 at the edge of a box.
  inputfile >> frame.box.xmin >> frame.box.xmax;
  // std::cout << frame.box.xmin << std::endl;
  // std::cout << "box min/max\n";
  // std::cout << frame.box.xmax << std::endl;
  frame.box.xlen = frame.box.xmax - frame.box.xmin;
  frame.box.xmin = -frame.box.xlen/2.0;
  frame.box.xmax = frame.box.xlen/2.0;
  inputfile >> frame.box.ymin >> frame.box.ymax;
  frame.box.ylen = frame.box.ymax - frame.box.ymin;
  frame.box.ymin = -frame.box.ylen/2.0;
  frame.box.ymax = frame.box.ylen/2.0;
  inputfile >> frame.box.zmin >> frame.box.zmax;
  frame.box.zlen = frame.box.zmax - frame.box.zmin;
  frame.box.zmin = -frame.box.zlen/2.0;
  frame.box.zmax = frame.box.zlen/2.0;
  // std::cout << "xlen = " << frame.box.xlen << std::endl;
  // std::cout << "ylen = " << frame.box.ylen << std::endl;
  // std::cout << "zlen = " << frame.box.zlen << std::endl;
  inputfile >> temp >> temp >> temp >> temp >> temp >> temp >> temp;
  frame.pts.resize(num_atoms);
  frame.points.resize(num_atoms);
  frame.atms.resize(num_atoms);
  frame.num_atoms = num_atoms;

  for (size_t i = 0; i < num_atoms; i++) {
    inputfile >> n >> nt >> xc >> yc >> zc;
    // std::cout << n << std::endl;
    frame.atms[n-1].atom_num = n;
    frame.atms[n-1].atom_type = nt;
    x = xc*frame.box.xlen + frame.box.xmin;
    y = yc*frame.box.ylen + frame.box.ymin;
    z = zc*frame.box.zlen + frame.box.zmin;
    // frame.pts[n-1].x = x;
    // frame.pts[n-1].y = y;
    // frame.pts[n-1].z = z;
    frame.points[n-1] = {x,y,z};
  }
  getline(inputfile, temp);
}

void Trajectory::skipFrames(const size_t sframes) {
  std::string temp;
  size_t slines = sframes*(Trajectory::header + 4 + num_atoms);
  for (size_t i = 0; i < slines; i++) {
    getline(inputfile, temp);
  }
  // statement below lets us check skipFrames functionality
  // std::cout << temp << '\n';
}




