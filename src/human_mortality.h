// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef SRC_HUMAN_MORTALITY_H_
#define SRC_HUMAN_MORTALITY_H_
#include <vector>
#include "Random.h"

std::vector<size_t> sample_mothers(const std::vector<bool>&,
                                   const std::vector<size_t>&,
                                   const std::vector<size_t>&,
                                   RandomInterface*);

#endif
