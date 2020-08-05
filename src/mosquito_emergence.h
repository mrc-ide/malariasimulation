/*
 * mosquito_emergence.h
 *
 *  Created on: 4 Aug 2020
 *      Author: gc1610
 */

#ifndef SRC_MOSQUITO_EMERGENCE_H_
#define SRC_MOSQUITO_EMERGENCE_H_

#include <individual.h>

Rcpp::XPtr<process_t> create_mosquito_emergence_process_cpp(std::string, std::string, std::string);

#endif /* SRC_MOSQUITO_EMERGENCE_H_ */
