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
    _values.insert({time, value});
    if (max_size != -1) {
        while(_values.size() > max_size) {
            _values.erase(_values.begin());
        }
    }
}

double Timeseries::at(double time, bool linear) const {
    auto it = _values.lower_bound(time);
    if (it == _values.end()) {
        if (_values.size() > 0 && !linear) {
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
        if (it == _values.begin()) {
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

const std::map<double, double>& Timeseries::values() {
    return _values;
}

void Timeseries::set_values(std::map<double, double> values) {
    _values = std::move(values);
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

//[[Rcpp::export]]
Rcpp::List timeseries_save_state(Rcpp::XPtr<Timeseries> timeseries) {
    std::vector<double> timesteps;
    std::vector<double> values;
    for (const auto& entry: timeseries->values()) {
        timesteps.push_back(entry.first);
        values.push_back(entry.second);
    }
    return Rcpp::DataFrame::create(
        Rcpp::Named("timestep") = timesteps,
        Rcpp::Named("value") = values
    );
}

//[[Rcpp::export]]
void timeseries_restore_state(Rcpp::XPtr<Timeseries> timeseries, Rcpp::List state) {
    std::vector<double> timesteps = state["timestep"];
    std::vector<double> values = state["value"];
    if (timesteps.size() != values.size()) {
        Rcpp::stop("Bad size");
    }

    std::map<double, double> values_map;
    for (size_t i = 0; i < timesteps.size(); i++) {
        values_map.insert({timesteps[i], values[i]});
    }
    timeseries->set_values(std::move(values_map));
}
