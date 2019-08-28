#CXX = g++ -O3 -std=c++11 -mavx -mfma
CXX = icpc  -O3 -std=c++11 
#CXX = g++ -O3 -std=c++11 -fopenmp -mavx -mfma
#CXX = clang++-8 -O3
#CXX = g++ -g

PRO = ising

OBJ =   io.o init_conf.o bond_energy.o sweep.o

$(PRO) : main.o $(OBJ)
	 $(CXX) $(OBJ) main.o -o $(PRO) -L/storagehome/zhangshuai/soft/lib/boost/lib -lboost_program_options

main.o : main.cpp
	$(CXX) -c main.cpp

io.o : io.cpp 
	$(CXX) -c io.cpp

init_conf.o : init_conf.cpp
	$(CXX) -c init_conf.cpp

bond_energy.o : bond_energy.cpp
	$(CXX) -c bond_energy.cpp
sweep.o : sweep.cpp
	$(CXX) -c sweep.cpp

clean :
	rm *.o ising

