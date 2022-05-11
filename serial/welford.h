#include <cstdlib>
#include <cmath>
#include "structures.h"
#include "traj_reader.h"


class RunningStat
{
  public:
    RunningStat() : m_n(0) {}
    void Clear()
    {
      m_n = 0;
    }
    void Push(double x)
    {
      m_n++;
      // See Knuth TAOCP vol 2, 3rd edition, page 232
      if (m_n == 1)
      {
        m_oldM = m_newM = x;
        m_oldS = 0.0;
      }
      else
      {
        m_newM = m_oldM + (x - m_oldM)/m_n;
        m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
        // set up for next iteration
        m_oldM = m_newM;
        m_oldS = m_newS;
      }
    }
    size_t NumDataValues() const
    {
      return m_n;
    }
    double Mean() const
    {
      return (m_n > 0) ? m_newM : 0.0;
    }
    double Variance() const
    {
      return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
    }
    double StandardDeviation() const
    {
      return sqrt( Variance() );
    }

  private:
    size_t m_n;
    double m_oldM, m_newM, m_oldS, m_newS;
};

//Welford method expanded to 3d point cloud
// https://www.johndcook.com/blog/standard_deviation/
template <typename T>
void variance00WK(std::string &filename, int num_atoms, int num_frames,
    int num_skipframes, resultSet<T> &result) {
  double xa, ya, za, xb, yb, zb, rb, xlen, ylen, zlen, old_avgx, old_avgy,
    old_avgz, varx, vary, varz, diff_sqrd, nsamples, variance, old_avgr,
    new_avgr, varrb, old_rb;
  double invRootThree = 1/pow(3,0.5);
  simFrame<double> frame0;
  simFrame<double> frame;
  Trajectory traj(filename, num_atoms, 5);
  traj.skipFrames(num_skipframes);
  traj.getNextFrame(frame);
  frame0 = frame;
  result.avg.points.resize(num_atoms, std::vector<double>(3));
  result.avg.num_atoms = num_atoms;
  result.avg.box.xmin = frame.box.xmin;
  result.avg.box.xmax = frame.box.xmax;
  result.avg.box.ymin = frame.box.ymin;
  result.avg.box.ymax = frame.box.ymax;
  result.avg.box.zmin = frame.box.zmin;
  result.avg.box.zmax = frame.box.zmax;

  xlen = frame.box.xlen;
  ylen = frame.box.ylen;
  zlen = frame.box.zlen;

  std::vector<double> len = {xlen,ylen,zlen};


  // for all points in first frame place them as initial values of the average
  // frame.
  // Points are already scaled by Trajectory::getNextFrame method
  for (int j = 0; j < num_atoms; j++) {
    //cout << frame0.pts[j].x << std::endl;
    result.avg.points[j][0] = frame0.points[j][0];
    result.avg.points[j][1] = frame0.points[j][1];
    result.avg.points[j][2] = frame0.points[j][2];
  }
  diff_sqrd = 0;
  varx = 0;
  vary = 0;
  varz = 0;

  double div_3 = 1.0/3.0;


  // std::vector<std::vector<double>

  for (int i = 1; i < num_frames; i++) {
    traj.getNextFrame(frame);

    int k = 0;
    for (int j = 0; j < num_atoms; j++) {
      // we want to make sure that if an atom moves through a periodic
      // boundary that we apply the minium image criterium and 'lock'
      // the atom to the side of the boundary it was on in the first
      // frame before we sample its position.
      xa = frame0.points[j][0];
      xb = frame.points[j][0];
      ya = frame0.points[j][1];
      yb = frame.points[j][1];
      za = frame0.points[j][2];
      zb = frame.points[j][2];
      if ((xa - xb) < (-xlen*0.5)){
	      xb = xb - xlen;
      }
      if ((xa - xb) > (xlen*0.5)){
	      xb = xb + xlen;
      }
      if ((ya - yb) < (-ylen*0.5)){
	      yb = yb - ylen;
      }
      if ((ya - yb) > (ylen*0.5)){
	      yb = yb + ylen;
      }
      if ((za - zb) < (-zlen*0.5)){
	      zb = zb - zlen;
      }
      if ((za - zb) > (zlen*0.5)){
	      zb = zb + zlen;
      }

      rb = pow(xb*xb + yb*yb + zb*zb, 0.5)*invRootThree;


      // welford code for reference
      // m_newM = m_oldM + (x - m_oldM)/m_n;
      // m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
      // // set up for next iteration
      // m_oldM = m_newM;
      // m_oldS = m_newS;

      // store old values
      old_avgx = result.avg.points[j][0];
      old_avgy = result.avg.points[j][1];
      old_avgz = result.avg.points[j][2];
      // xdiff = xb -

      old_avgr = pow(old_avgx*old_avgx + old_avgy*old_avgy
      + old_avgz*old_avgz, 0.5)*invRootThree;

      // calculate running average
      result.avg.points[j][0] = old_avgx + (xb - old_avgx)/i;
      result.avg.points[j][1] = old_avgy + (yb - old_avgy)/i;
      result.avg.points[j][2] = old_avgz + (zb - old_avgz)/i;

      new_avgr = old_avgr + (rb - old_avgr)/i;

      //calculate running variance
      varx = varx + (xb - old_avgx)*(xb - result.avg.points[j][0]);
      vary = vary + (yb - old_avgy)*(yb - result.avg.points[j][1]);
      varz = varz + (zb - old_avgz)*(zb - result.avg.points[j][2]);

      varrb = varrb + (rb - old_avgr)*(rb - new_avgr);

      diff_sqrd = diff_sqrd + ((rb - old_avgr)*(rb - new_avgr));


    }
  }
// we want to make sure the points in our average frame are
// contained in the original simulation box limits so we apply
// the minium image criterium

  for (int j = 0; j < num_atoms; j++){
    for (int k = 0; k < 3; k++){
      if (result.avg.points[j][k] < result.avg.box.xmin) {
  result.avg.points[j][k] = result.avg.points[j][k] + len[k];
      }
      if (result.avg.points[j][k] < result.avg.box.xmin) {
  result.avg.points[j][k] = result.avg.points[j][k] + len[k];
      }
    }
  }
  result.avg.box.xlen = xlen;
  result.avg.box.ylen = ylen;
  result.avg.box.zlen = zlen;


  nsamples = num_frames*num_atoms;
  result.variance = varrb/(nsamples-1);
  result.var_xyz.push_back(varx/(nsamples-1));
  result.var_xyz.push_back(vary/(nsamples-1));
  result.var_xyz.push_back(varz/(nsamples-1));

  return;
}



