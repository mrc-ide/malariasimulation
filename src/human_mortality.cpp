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

std::vector<bool> sampleable(const variable_vector_t& birth, size_t timestep) {
    auto ret = std::vector<bool>(birth.size());
    auto upper = 35 * 365;
    auto lower = 15 * 365;
    for (auto i = 0u; i < birth.size(); ++i) {
        auto age = timestep - birth[i];
        ret[i] = age >= lower && age <= upper;
    }
    return ret;
}

std::vector<size_t> intersect(
    const variable_vector_t& a,
    const individual_index_t& b) {
    auto res = std::vector<size_t>();
    for (auto i : a) {
        if (b[i] == 1.) {
            res.push_back(i);
        }
    }
    return res;
}


std::vector<size_t> calculate_died_from_severe(
                const variable_vector_t& is_severe,
                const individual_index_t& diseased,
                double v,
                const individual_index_t& treated,
                double fvt
            ) {
    auto at_risk  = intersect(is_severe, diseased);
    auto treated_at_risk  = intersect(is_severe, treated);
    auto died_at_risk = random->bernoulli(
        at_risk.size(),
        v
    );
    auto treated_died_at_risk = random->bernoulli(
        treated_at_risk.size(),
        fvt
    );
    auto died_from_severe = std::vector<size_t>(
        died_at_risk.size() + treated_died_at_risk.size()
    );
    for (auto i = 0u; i < died_at_risk.size(); ++i) {
        died_from_severe[i] = at_risk[died_at_risk];
    }
    for (auto i = died_at_risk.size(); i < died_from_severe.size(); ++i) {
        died_from_severe[i] = treated_at_risk[treated_died_at_risk];
    }
    return died_from_severe;
}

//' @title Human mortality process
//' @description
//' This is the process for human mortality, it defines which humans die from
//' natural causes and severe infection and replaces dead individuals with
//' newborns.
//' @param random a pointer to an interface for RNG
Rcpp::XPtr<process_t> create_mortality_process(RandomInterface* random) {
    return Rcpp::XPtr<process_t>(new process_t([=] (ProcessAPI& api) {
        const auto& parameters = api.get_parameters();
        auto died = random->bernoulli(
            parameters.at("human_population")[0],
            parameters.at("mortality_rate")[0]
        );
        api.render("natural_deaths", died.size());

        if (parameters.at("severe_enabled")[0] == 1.) {
            auto died_from_severe = calculate_died_from_severe(
                api.get_variable("human", "is_severe"),
                api.get_state("human", "D"), //This is zero based because individual is clever
                parameters.at("v")[0],
                api.get_state("human", "Tr"),
                parameters.at("fvt")[0]
            );
            //api$render('severe_deaths', length(severe_deaths))
            //died <- union(died, severe_deaths)
        }


    }), true);
}
