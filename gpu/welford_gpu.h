#include "structures.h"

#ifndef WELFORD_GPU_H
#define WELFORD_GPU_H

void variance00WK_gpu(std::string &filename,
		      int num_atoms,
		      int num_frames,
		      int num_skipframes,
		      resultSet<double> &result);

#endif
