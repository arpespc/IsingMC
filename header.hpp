#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include <vector>
#include <map>
#include <string>
#include <Eigen/Eigen>

void read_para(std::map<std::string, double> &params);
// void init_conf(Eigen::MatrixXi &spin_conf, const int Nx, const int Ny, );
void init_conf(Eigen::MatrixXi &spin_conf, const int Nx, const int Ny, boost::random::mt19937 & seed_spin);
double bond_energy(std::array<Eigen::MatrixXi, 2> &spin_lattice, int s, int i, int j, std::map<std::string, double> & params);
void sweep(std::array<Eigen::MatrixXi, 2> & spin_lattice, const int Nx, const int Ny, std::map<std::string, double> & params,  boost::random::mt19937 & rng);
