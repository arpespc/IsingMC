#ifndef _RANDOM_H_
#define _RANDOM_H_
#include <boost/random.hpp>
#endif

#ifndef _STRING_H_
#define _STRING_H_
#include <string>
#endif

#ifndef _EIGEN_H_
#define _EIGEN_H_
#include <Eigen/Eigen>
#endif

void read_para(std::map<std::string, double> &params);
// void init_conf(Eigen::MatrixXi &spin_conf, const int Nx, const int Ny, );
void init_conf(Eigen::MatrixXd &spin_conf, const int Nx, const int Ny, boost::random::mt19937 & seed_spin);
double bond_energy(std::array<Eigen::MatrixXd, 2> &spin_lattice, int s, int i, int j, std::map<std::string, double> & params);
void sweep(std::array<Eigen::MatrixXd, 2> & spin_lattice, const int Nx, const int Ny, std::map<std::string, double> & params,  boost::random::mt19937 & rng);
