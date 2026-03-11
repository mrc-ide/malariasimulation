/*
 * individual_level_output.cpp
 *
 *  Created on: 12 Aug 2025
 *      Author: JasonRWood
 */
#include "individual_level_output.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unordered_map>

//[[Rcpp::export]]
void print_to_csv(
    const std::string filename,
    const int timestep,
    const std::vector<int> personal_indicies,
    const int process,
    const std::vector<int> categories,
    const int turnon_time
){  
    if (timestep < turnon_time){
        return;
    }
    // std::unordered_map<std::string, int> umap = {
    //     {"S", 0},
    //     {"U", 1},
    //     {"A", 2},
    //     {"D", 3},
    //     {"Tr", 4}
    // };
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app);
    for (auto i = 0u; i < personal_indicies.size(); ++i){
        // outfile << timestep << "," << personal_indicies[i] << ","  << process << "," << umap[categories[i]] << "\n";
        outfile << timestep << "," << personal_indicies[i] << ","  << process << "," << categories[i] << "\n";
    }
    outfile.close();
    return;
}