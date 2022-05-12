# welford_gpu

First step to bring the Welford implementation from SFunk to GPU.

For a set of available silcon data consisting of 54872 atoms and 16905 frames the gpu code improves performance by ~15%. A larger iron system consisting of 93312 atoms and 3844 frames shows a similar performance boost with the gpu code.

This code is a first step toward realizing real gains in performance and is not considered to be a complete work at this time. 5-12-2022

