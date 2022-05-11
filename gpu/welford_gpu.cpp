#include <cstdlib>
#include <cmath>
#include "structures.h"
#include "traj_reader.h"
#include "welford_gpu.h"


// error checking macro lifted from OCLF nVidia training class
#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

      // welford code for reference
      // m_newM = m_oldM + (x - m_oldM)/m_n;
      // m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
      // // set up for next iteration
      // m_oldM = m_newM;
      // m_oldS = m_newS;

__global__ void box_adjust(Point *p,
		           const Point *p0,
			   Box<double> *b,
			   int n) {
  
  int i = threadIdx.x+blockDim.x*blockIdx.x;
  if (i < n){
    if ((p0[i].x - p[i].x) < (-b->xlen*0.5)){
      p[i].x = p[i].x - b->xlen;
    }
    if ((p0[i].x - p[i].x) > (b->xlen*0.5)){
      p[i].x = p[i].x + b->xlen;
    }
    if ((p0[i].y - p[i].y) < (-b->ylen*0.5)){
      p[i].y = p[i].y - b->ylen;
    }
    if ((p0[i].y - p[i].y) > (b->ylen*0.5)){
      p[i].y = p[i].y + b->ylen;
    }
    if ((p0[i].z - p[i].z) < (-b->zlen*0.5)){
      p[i].z = p[i].z - b->zlen;
    }
    if ((p0[i].z - p[i].z) > (b->zlen*0.5)){
      p[i].z = p[i].z + b->zlen;
    }
  }
}

__global__ void welford_kernel(const Point *pts,
			       Point *avg_pts,
			       int m_n,
			       int dn) {
  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  if (idx < dn) {
    avg_pts[idx].x = avg_pts[idx].x + (pts[idx].x - avg_pts[idx].x)/m_n;
    avg_pts[idx].y = avg_pts[idx].y + (pts[idx].y - avg_pts[idx].y)/m_n;
    avg_pts[idx].z = avg_pts[idx].z + (pts[idx].z - avg_pts[idx].z)/m_n;  
  }
}	

//Welford method expanded to 3d point cloud
// https://www.johndcook.com/blog/standard_deviation/
void variance00WK_gpu(std::string &filename,
		      int num_atoms,
		      int num_frames,
		      int num_skipframes,
		      resultSet<double> &result) {
  
  simFrame<double> frame0;
  simFrame<double> frame;
  
  // create pointer and allocate memory for data
  Point *d_pts, *d_pts0, *d_avg_pts;
  int N, M;

  Box<double> *box, *d_box;
  N = num_atoms;

  Trajectory traj(filename, num_atoms, 5);
  traj.skipFrames(num_skipframes);
  traj.getNextFrame(frame);
  frame0 = frame;
  
//  std::cout << "read first frame" << std::endl;

  cudaMalloc((void**)&d_pts, sizeof(Point)*num_atoms);
  cudaMalloc((void**)&d_avg_pts, sizeof(Point)*num_atoms);
  cudaMalloc((void**)&d_pts0, sizeof(Point)*num_atoms);
  cudaMalloc((void**)&d_box, sizeof(Box<double>)*num_atoms);

//  std::cout << "allocated device memory" << std::endl;
  
  cudaMemcpy(d_avg_pts, frame0.pts.data(), sizeof(Point)*num_atoms, cudaMemcpyHostToDevice);
  cudaCheckErrors("copy average points");
  cudaMemcpy(d_pts0, frame0.pts.data(), sizeof(Point)*num_atoms, cudaMemcpyHostToDevice);
  cudaCheckErrors("copy pts0");
  box = &frame0.box;
  cudaMemcpy(d_box, box, sizeof(Box<double>), cudaMemcpyHostToDevice);
  cudaCheckErrors("copy box");
 
//  std::cout << "copied memory to device" << std::endl;


  M = 512;
  int T = (N + M - 1)/M;
  
  for (int i = 2; i < num_frames; i++) {
    traj.getNextFrame(frame);
    
//    std::cout << "read frame " << i << std::endl;

    cudaMemcpy(d_pts, frame.pts.data(), sizeof(Point)*num_atoms, cudaMemcpyHostToDevice);
    
    box_adjust<<<T, M>>>(d_pts, d_pts0, d_box, N);
    cudaCheckErrors("box adjust kernel");
    cudaDeviceSynchronize();
    
//    std::cout << "adjusted box for frame " << i << std::endl;

    welford_kernel<<<T, M>>>(d_pts, d_avg_pts, i, N);
    cudaCheckErrors("welford kernel");
    cudaDeviceSynchronize();

//    std::cout << "calcuted welford for frame " << i << std::endl;

  }
  
  cudaFree(d_pts);
  cudaFree(d_avg_pts);
  cudaFree(d_pts0);
  cudaFree(d_box);
}

