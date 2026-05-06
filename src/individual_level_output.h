/*
 * individual_level_output.h
 *
 *  Created on: 12 Aug 2025
 *      Author: JasonRWood
 */

#ifndef SRC_INDIVIDUAL_LEVEL_OUTPUT_H_
#define SRC_INDIVIDUAL_LEVEL_OUTPUT_H_
#include <vector>
#include <string>

void print_to_csv(
    const std::string filename,
    const int timestep,
    const std::vector<int> personal_indicies,
    const std::string process,
    const std::vector<int> categories,
    const int turnon_time
);

void print_for_snapshot(
    const std::string filename,
    const int timestep,
    const std::vector<int> personal_indicies,
    const std::vector<int> ages,
    const std::vector<int> categories
);

#endif /* SRC_INDIVIDUAL_LEVEL_OUTPUT_H_ */
