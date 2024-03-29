#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <map>
#include "header.hpp"

boost::random::mt19937 rng(std::time(0));
// scientific constant
const double kb = 8.617333262145 * 1e-5 ; // Boltzman constant eV / K

int main(int argc, char* argv[]) {
    
    // comand line options
    namespace  po = boost::program_options;
    po::options_description opts("options");
    opts.add_options()
        ("help,h", "help info")
        ("temp,T", po::value<double>() ,"MC temprature, override the value from  the input file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);

    if (vm.count("help")) {
        std::cout << opts << std::endl;
        return -1;
    }

    // imput params
    std::map<std::string, double> params;
    read_para(params);
    const int Nx = (int) params["Nx"]; // number of size
    const int Ny = (int) params["Ny"]; // number of size
    const long Num_thermo = params["Num_thermo"];
    const long Num_sample = params["Num_sample"];
    const double S = params["S"]; // muB
    double t = 0;
    if (vm.count("temp")) {
        t = vm["temp"].as<double>();
    } else {
        t = params["T"];
    }
    const double T = t; // a little ugly...
    const double beta = 1.0 / (kb * T);
    params.insert(std::pair<std::string, double>("beta", beta));
    const bool snapshot = (int) params["snapshot"];
    int snap_interval =  0;
    std::string snap_folder;
    if (snapshot) {
        snap_folder = "./snapshot/";
        std::string conmand = "mkdir -p " + snap_folder;
        int mkdir = std::system(conmand.c_str());
        snap_interval = (int) params["snap_interval"];
    }
    const double num_site =  Nx * Ny;
    const double Anneal_blow_T = params["Anneal_blow_T"];
    const double Anneal_T = params["Anneal_T"];
    const double Anneal_num = params["Anneal_num"];

    // out put the basic input infomation
    std::cout << "############################## Critical Input Infomation #################################" << std::endl;
    std::cout << "MC Temperature      : " << T << std::endl;
    std::cout << "lattice size  Nx    : " << Nx << std::endl;
    std::cout << "lattice size  Ny    : " << Ny << std::endl;
    std::cout << "Number of site      : " << Ny << std::endl;
    std::cout << "snapshot            : " << snapshot << std::endl;
    std::cout << "snap_interval       : " << snap_interval << std::endl;
    std::cout << "Thermo number       : " << Num_thermo << std::endl;
    std::cout << "MC Sample number    : " << Num_sample << std::endl;
    std::cout << "##########################################################################################" << std::endl;
    std::cout << std::endl;

    double Ms_sample = 0; // the magnetization  of one sample
    // initial the the basic output 
    Eigen::MatrixXd M_mat = Eigen::MatrixXd::Zero(Nx, Ny); // The total magetic moment of theh whole lattice = spin1 + spin2
    double Ms = 0; // magnetization
    double Ms2 = 0;
    double chi = 0; // the magnetic susceptbility
    double errorbar = 0;

    // Initial Spin  Configure
    std::array<Eigen::MatrixXd, 2> spin_lattice;
    Eigen::MatrixXd spin_conf = Eigen::MatrixXd::Zero(Nx, Ny);
    Eigen::MatrixXd spin_conf1 = Eigen::MatrixXd::Ones(Nx, Ny); // The spin configuration of of first sub lattice
    Eigen::MatrixXd spin_conf2 = Eigen::MatrixXd::Ones(Nx, Ny); // The spin configuration of the second sub lattice
    Eigen::MatrixXd M_mat_sample = Eigen::MatrixXd::Zero(Nx, Ny); // The total magetic moment of theh whole lattice = spin1 + spin2 of this sample

    init_conf(spin_conf, Nx, Ny, rng);
    spin_lattice[0] =  S * spin_conf;
    init_conf(spin_conf, Nx, Ny, rng);
    spin_lattice[1] = S * spin_conf;

    // Annealing 
    if ( T < Anneal_blow_T ) {
        std::cout << "############################## Annealing #################################################" << std::endl;
        params["T"] = Anneal_T;
        std::cout << "Anneal_blow_T             : " << Anneal_blow_T << std::endl;
        std::cout << "Anneal_T                  : " << Anneal_T << std::endl;
        std::cout << "Anneal_num                : " << Anneal_num << std::endl;
        std::cout << "Now T is                  : " << params["T"] << std::endl;
        for (int anneal_step = 0; anneal_step <= (Anneal_num - 1); anneal_step++) {
            sweep(spin_lattice, Nx, Ny, params, rng);
            if ((anneal_step % 5000) == 0)  {
                std::cout << "Annealing step: " << anneal_step << std::endl;
            }
        }
        params["T"] = T;
        std::cout << "After Anneal, T is back to : " << params["T"] << std::endl;
        std::cout << "################################ Annealing finished #######################################" << std::endl;
        std::cout << std::endl;

    }

    // Thermolazition to equilition
    std::cout << "################################# Thermolization #########################################" << std::endl;
    for (int thermo_step = 0; thermo_step <= (Num_thermo - 1); thermo_step++) {
        sweep(spin_lattice, Nx, Ny, params, rng);
        if ((thermo_step % 5000) == 0)  {
            std::cout << "Thermo step: " << thermo_step << std::endl;
        }

        if (snapshot) {
            if  ((thermo_step % snap_interval) == 0) {
                std::string snap_file = snap_folder +  "M_thermo" + "_" + std::to_string(thermo_step);
                std::ofstream fob(snap_file);
                fob << spin_lattice[0] + spin_lattice[1] << std::endl;
            }
        }
    }
    std::cout << "################################# Thermo finished ########################################" << std::endl;

    // MC process
    std::cout << "################################# MC Sampling ############################################" << std::endl;
    for (int step = 0; step <= (Num_sample - 1); step++) {
        sweep(spin_lattice, Nx, Ny, params, rng);
        // After casting, the original type of Matrix will not change;
        M_mat_sample =  spin_lattice[0].cast<double>() +  spin_lattice[1].cast<double>();
        M_mat =  M_mat_sample + M_mat;
        Ms_sample =   M_mat_sample.sum() / num_site ;
        Ms +=  Ms_sample;
        Ms2 +=  std::pow(Ms_sample, 2);
        if ((step % 5000) == 0)  {
            std::cout << "MC step: " << step << std::endl;
        }

        if (snapshot) {
            if  ((step % snap_interval) == 0) {
                std::string snap_file = snap_folder +  "M_mc" + "_" + std::to_string(step);
                std::ofstream fob(snap_file);
                fob << spin_lattice[0] + spin_lattice[1] << std::endl;
            }
        }
    }    
    std::cout << "################################ Sampling finished #######################################" << std::endl;

    // calculate  the output values
    M_mat = M_mat / Num_sample;
    Ms = Ms / Num_sample;
    Ms = std::abs(Ms); // just for plot, because of the Z2  symmetry of spin
    Ms2 = Ms2 / Num_sample;  // mean value of squre of magnetization  of eche sampling
    errorbar = Ms2 - std::pow(Ms, 2);
    chi = beta * errorbar;
    errorbar  = std::pow(errorbar, 0.5);
    

    std::string M_mat_folder = "./M_mat/";
    std::string M_mat_file = M_mat_folder +  "M_mat_" + std::to_string(T);
    std::string conmand = "mkdir -p " + M_mat_folder;
    int mkdir = std::system(conmand.c_str());
    std::ofstream fob(M_mat_file);
    fob << M_mat << std::endl;

    boost::format  fmt_Msdat("%7.2f       %12.9f     %14.9f     %14.9f ");
    fmt_Msdat%T%std::abs(Ms)%errorbar%chi;
    std::ofstream Ms_file("Ms.dat", std::ios::app);
    Ms_file << fmt_Msdat << std::endl;

    return 0;
}
