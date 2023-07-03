/*
 * timeseries.h
 *
 *  Created on: 08 Apr 2021
 *      Author: gc1610
 */

#ifndef SRC_TIMESERIES_H_
#define SRC_TIMESERIES_H_

#include <map>

class Timeseries {
private:
    std::map<double,double> values;
    size_t max_size;
    void clean();
    bool has_default;
    double default_value;
public:
    Timeseries();
    Timeseries(size_t);
    Timeseries(size_t, double);
    void push(double, double);
    double at(double, bool = true) const;
};

#endif /* SRC_TIMESERIES_ */
