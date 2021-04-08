/*
 * history.cpp
 *
 *  Created on: 08 Apr 2021
 *      Author: gc1610
 */

#include <Rcpp.h>
#include <algorithm>
#include "history.h"

History::History() : max_size(-1) {}
History::History(size_t max_size) : max_size(max_size) {}

void History::push(double time, double value) {
    values.insert({time, value});
    if (max_size != -1) {
        while(values.size() > max_size) {
            values.erase(values.begin());
        }
    }
}

double History::at(double time) const {
    auto it = values.lower_bound(time);
    if (it == values.end()) {
        Rcpp::stop("`time` is later than the stored history");
    }

    // Check if we've landed on the exact time
    if (it->first == time) {
        return it->second;
    }

    // Find the element before
    auto after_element = *it;
    while(it->first > time) {
        // Check if we're at the start of the history
        if (it == values.begin()) {
            Rcpp::stop("`time` is before our stored history");
        }
        it--;
    }

    // Interpolate
    return it->second + (
        (time - it->first) / (after_element.first - it->first)
    ) * (after_element.second - it->second);
}
