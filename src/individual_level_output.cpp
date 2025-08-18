/*
 * mosquito_biology.cpp
 *
 *  Created on: 12 Aug 2025
 *      Author: JasonRWood
 */
#include "individual_level_output.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

//[[Rcpp::export]]
void print_to_csv(
    const std::string filename,
    const int timestep,
    const std::vector<int> personal_indicies,
    const std::string process,
    const std::vector<std::string> categories
){  
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app);
    for (auto i = 0u; i < personal_indicies.size(); ++i){
        outfile << timestep << "," << personal_indicies[i] << ","  << process << "," << categories[i] << "\n";
    }
    outfile.close();
    return;
}