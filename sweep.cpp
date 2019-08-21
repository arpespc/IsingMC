#include "header.hpp"

void sweep(std::array<Eigen::MatrixXi, 2> & spin_lattice, const int Nx, const int Ny, std::map<std::string, double> & params,  boost::random::mt19937 & rng) {
    double beta = params["beta"];
    for (int s = 0; s <= 1; s++) {
        for (int i = 0; i <= Nx - 1; i++) {
            for (int j = 0; j <= Ny - 1; j++) {
                double E = bond_energy(spin_lattice, s, i, j, params);
                double dE = -E - E; // J (-M1) * M2 - J * M1 * M2 = E2 - E1
                double p_tran = std::exp(-beta * dE);
                // Random number generator
                boost::random::uniform_real_distribution<double> prob(0, 1);
                if ( prob(rng) < p_tran) {
                    spin_lattice[s](i, j) = - spin_lattice[s](i, j);
                }
            }
        }
    }

}
