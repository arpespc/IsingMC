#define _STRING_H_
#include "header.hpp"

void init_conf(Eigen::MatrixXi &spin_conf, const int Nx, const int Ny, boost::random::mt19937 & seed_spin) {
    boost::random::uniform_int_distribution<> rand_spin(0,1);
    // boost::variate_generator<boost::mt19937, boost::uniform_int<>> rand_spin(seed_spin, dist);

    for (int i =0; i <= Nx - 1; i++) {
        for (int j = 0; j <= Ny - 1; j++){
            int s = 2 * rand_spin(seed_spin) - 1;
            spin_conf(i, j) = s;
        }
    }
}
