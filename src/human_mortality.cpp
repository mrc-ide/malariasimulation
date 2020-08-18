// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "human_mortality.h"

std::vector<size_t> sample_mothers(const std::vector<bool>& sampleable,
                                   const std::vector<size_t>& groups,
                                   const size_t n_groups,
                                   const std::vector<size_t>& died,
                                   RandomInterface* random) {
    std::vector<std::vector<size_t>> index(n_groups);
    for (auto i = 0u; i < sampleable.size(); ++i) {
        if (sampleable[i]) {
            index[groups[i] - 1].push_back(i);
        }
    }

    std::vector<size_t> n_died(n_groups);
    std::vector<std::vector<size_t>> index_dest(died.size());
    for (auto i = 0u; i < died.size(); ++i) {
        size_t group = groups[died[i]];
        n_died[group - 1]++;
        index_dest[group - 1].push_back(i);
    }

    std::vector<size_t> ret(died.size());
    for (auto i = 0u; i < n_groups; ++i) {
        std::vector<size_t> take = random->sample(index[i].size(), n_died[i]);
        const auto& index_dest_i = index_dest[i];
        for (auto j = 0u; j < take.size(); ++j) {
            ret[index_dest_i[j]] = index[i][take[j]];
        }
    }

    return ret;
}
