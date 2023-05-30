/*
 * timeseries.cpp
 *
 *  Created on: 08 Apr 2021
 *      Author: gc1610
 */

#include <Rcpp.h>
#include <algorithm>
#include "timeseries.h"

Timeseries::Timeseries() : max_size(-1), has_default(false) {}
Timeseries::Timeseries(size_t max_size) : max_size(max_size), has_default(false) {}
Timeseries::Timeseries(size_t max_size, double default_value)
    : max_size(max_size), has_default(true), default_value(default_value) {}

void Timeseries::push(double value, double time) {
    values.insert({time, value});
    if (max_size != -1) {
        while(values.size() > max_size) {
            values.erase(values.begin());
        }
    }
}

double Timeseries::at(double time, bool linear) const {
    auto it = values.lower_bound(time);
    if (it == values.end()) {
        if (values.size() > 0 && !linear) {
          it--;
          return it->second;
        }
        if (has_default) {
            return default_value;
        }
        Rcpp::stop("`time` is later than the stored timeseries");
    }

    // Check if we've landed on the exact time
    if (it->first == time) {
        return it->second;
    }

    // Find the element before
    auto after_element = *it;
    while(it->first > time) {
        // Check if we're at the start of the timeseries
        if (it == values.begin()) {
            if (has_default) {
                return default_value;
            }
            Rcpp::stop("`time` is before our stored timeseries");
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
Rcpp::XPtr<Timeseries> create_timeseries(size_t size, double default_value) {
    return Rcpp::XPtr<Timeseries>(new Timeseries(size, default_value), true);
}

//[[Rcpp::export]]
double timeseries_at(Rcpp::XPtr<Timeseries> timeseries, double timestep, bool linear) {
    return timeseries->at(timestep, linear);
}

//[[Rcpp::export]]
void timeseries_push(Rcpp::XPtr<Timeseries> timeseries, double value, double timestep) {
    return timeseries->push(value, timestep);
}
