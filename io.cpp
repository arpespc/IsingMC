#include "header.hpp"

void read_para( std::map<std::string, double> &params){
    std::ifstream param_file("./params");
    if ( !param_file.is_open() ) {
        std::cerr << "There is not params file !!!" << std::endl;
        exit(-1);
    }
    std::string current_line;
    while (getline(param_file, current_line)) {
        boost::trim(current_line);
        std::vector<std::string> s_list;
        boost::split(s_list, current_line, boost::is_any_of(":=#"), boost::token_compress_on);
        s_list[0] = boost::trim_copy(s_list[0]);
        s_list[1] = boost::trim_copy(s_list[1]);
        params[s_list[0]] = stod(s_list[1]);
    }
}
