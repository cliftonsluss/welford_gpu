CC = g++
CFLAGS = -std=c++11
#CC	 = clang++

objects = main.o welford_gpu.o traj_reader.o read_config.o

all: $(objects)
	nvcc $(objects) -o gpu

%.o: %.cpp
	nvcc -x cu -I. -dc $< -o $@

clean:
	rm *.o gpu 
