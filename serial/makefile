CC = g++
CFLAGS = -std=c++11
#CC	 = clang++

serial: main.o traj_reader.o read_config.o 
	$(CC) -o $@ $^

main.o: main.cpp welford.h
	$(CC) $(CFLAGS) -c -g main.cpp 

traj_reader.o: traj_reader.cpp traj_reader.h
	$(CC) $(CFLAGS) -c -g traj_reader.cpp 

read_config.o: read_config.cpp read_config.h
	$(CC) $(CFLAGS) -c -g read_config.cpp


clean:
	rm *.o serial  
