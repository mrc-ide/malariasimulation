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
    const int timestep,
    const std::vector<int> personal_indicies,
    const std::string process,
    const std::vector<std::string> categories
);

#endif /* SRC_INDIVIDUAL_LEVEL_OUTPUT_H_ */
