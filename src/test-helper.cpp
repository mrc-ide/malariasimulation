/*
 * test-helper.cpp
 *
 *  Created on: 20 Aug 2020
 *      Author: gc1610
 */

#include "test-helper.h"

individual_index_t individual_index(size_t size, std::vector<size_t> i) {
    return individual_index_t(size, i.cbegin(), i.cend());
}
