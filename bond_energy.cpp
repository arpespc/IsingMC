#include "header.hpp"

double bond_energy (std::array<Eigen::MatrixXi, 2> &spin_lattice, int s, int i, int j, std::map<std::string, double> & params) {
    // for this lattice, we have two sub lattice, label by s, which blongs to {0, 1}
    // and i, j denote the site of sub lattice s.
    // so the spin of (i, j) site of s sub lattice is spin_lattice[s](i, j)

    // the exchange energy of nearest  and second nearest J1 and J2
    double J1 = params["J1"];
    double J2 = params["J2"];
    int s1 = (s + 1) % 2; // the other sub lattice index
    int Nx = params["Nx"];
    int Ny = params["Ny"];
    i = i % Nx;
    j = j % Ny;
    
    // the nearest neigbhor bond energy is
    double E1 = 0.0;
    E1 = -J1 * spin_lattice[s](i % Nx, j % Ny) * \
         ( spin_lattice[s]((i+1) % Nx, j % Ny) + spin_lattice[s]((i - 1 + Nx) % Nx, j));
    
    
    // the secend nearest bond energy is :
    double E2 = 0;
    E2 = - J2 * spin_lattice[s](i % Nx, j % Ny) * \
        ( spin_lattice[s1](i, j) + spin_lattice[s1]((i + 1) % Nx, j % Ny) + \
        spin_lattice[s1](i % Nx, (j+1) % Ny) + spin_lattice[s1]((i + 1) % Nx, (j + 1) % Ny) );
    
    return (E1 + E2);
}
