/*
 * history.h
 *
 *  Created on: 08 Apr 2021
 *      Author: gc1610
 */

#ifndef SRC_HISTORY_H_
#define SRC_HISTORY_H_

#include <map>

class History {
private:
    std::map<double,double> values;
    size_t max_size;
    void clean();
    bool has_default;
    double default_value;
public:
    History();
    History(size_t);
    History(size_t, double);
    void push(double, double);
    double at(double) const;
};

#endif /* SRC_HISTORY_H_ */
