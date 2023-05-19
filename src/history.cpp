/*
 * history.cpp
 *
 *  Created on: 08 Apr 2021
 *      Author: gc1610
 */

#include <Rcpp.h>
#include <algorithm>
#include "history.h"

History::History() : max_size(-1), has_default(false) {}
History::History(size_t max_size) : max_size(max_size), has_default(false) {}
History::History(size_t max_size, double default_value)
    : max_size(max_size), has_default(true), default_value(default_value) {}

void History::push(double value, double time) {
    values.insert({time, value});
    if (max_size != -1) {
        while(values.size() > max_size) {
            values.erase(values.begin());
        }
    }
}

double History::at(double time, bool linear) const {
    auto it = values.lower_bound(time);
    if (it == values.end()) {
        if (!linear) {
          it--;
          return it->second;
        }
        if (has_default) {
            return default_value;
        }
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
            if (has_default) {
                return default_value;
            }
            Rcpp::stop("`time` is before our stored history");
        }
        it--;
    }

    if(linear) {
      // Interpolate
      return it->second + (
          (time - it->first) / (after_element.first - it->first)
      ) * (after_element.second - it->second);
    }
    
    return it->second;
}

//[[Rcpp::export]]
Rcpp::XPtr<History> create_history(size_t size, double default_value) {
    return Rcpp::XPtr<History>(new History(size, default_value), true);
}

//[[Rcpp::export]]
double history_at(Rcpp::XPtr<History> history, double timestep, bool linear) {
    return history->at(timestep, linear);
}

//[[Rcpp::export]]
void history_push(Rcpp::XPtr<History> history, double value, double timestep) {
    return history->push(value, timestep);
}
