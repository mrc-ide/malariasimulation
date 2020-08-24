// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "human_mortality.h"
#include <unordered_set>

std::vector<size_t> sample_mothers(
    const std::vector<bool>& sampleable,
    const variable_vector_t& groups,
    const size_t n_groups,
    const std::vector<size_t>& died,
    RandomInterface* random
    ) {
    std::vector<std::vector<size_t>> index(n_groups);
    for (auto i = 0u; i < sampleable.size(); ++i) {
        if (sampleable[i]) {
            index[groups[i] - 1].push_back(i);
        }
    }

    std::vector<size_t> n_died(n_groups);
    std::vector<std::vector<size_t>> index_dest(n_groups);
    for (auto i = 0u; i < died.size(); ++i) {
        size_t group = groups[died[i]];
        n_died[group - 1]++;
        index_dest[group - 1].push_back(i);
    }

    std::vector<size_t> ret(died.size());
    for (auto i = 0u; i < n_groups; ++i) {
        const auto take = random->sample(index[i].size(), n_died[i]);
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
    const variable_vector_t& variable, //{0, 1, 0, 0...}
    const individual_index_t& index) { //{3, 6, 9...}
    auto res = std::vector<size_t>();
    for (auto i : index) {
        if (variable[i] == 1.) {
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
        double fvt,
        RandomInterface *random
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
        died_from_severe[i] = at_risk[died_at_risk[i]];
    }
    for (auto i = died_at_risk.size(); i < died_from_severe.size(); ++i) {
        died_from_severe[i] = treated_at_risk[treated_died_at_risk[i]];
    }
    return died_from_severe;
}

void union_vectors(std::vector<size_t>& full, const std::vector<size_t> others) {
    auto seen = std::unordered_set<size_t>(full.begin(), full.end());
    for (auto i : others) {
        if (seen.find(i) == seen.end()) {
            full.push_back(i);
        }
    }
}


//' @title Human mortality process
//' @description
//' This is the process for human mortality, it defines which humans die from
//' natural causes and severe infection and replaces dead individuals with
//' newborns.
//[[Rcpp::export]]
Rcpp::XPtr<process_t> create_mortality_process_cpp() {
    return create_mortality_process(&Random::get_instance());
}

Rcpp::XPtr<process_t> create_mortality_process(RandomInterface* random) {
    return Rcpp::XPtr<process_t>(new process_t([=] (ProcessAPI& api) {
        const auto& parameters = api.get_parameters();
        auto died = random->bernoulli(
            parameters.at("human_population")[0],
            parameters.at("mortality_rate")[0]
        );

        if (parameters.at("severe_enabled")[0] == 1.) {
            auto died_from_severe = calculate_died_from_severe(
                api.get_variable("human", "is_severe"),
                api.get_state("human", "D"), //This is zero based because individual is clever
                parameters.at("v")[0],
                api.get_state("human", "Tr"),
                parameters.at("fvt")[0],
                random
            );

            union_vectors(died, died_from_severe);
        }

        if (died.size() > 0) {
            auto timestep = api.get_timestep();
            // Calculate new maternal immunities
            auto groups = api.get_variable("human", "zeta_group");
            auto is_sampleable = sampleable(api.get_variable("human", "birth"), timestep);
            auto ica = api.get_variable("human", "ICA");
            auto iva = api.get_variable("human", "IVA");
            auto mothers = sample_mothers(
                is_sampleable,
                groups,
                parameters.at("n_heterogeneity_groups")[0],
                died,
                random
            );

            auto birth_icm = variable_vector_t(mothers.size());
            auto birth_ivm = variable_vector_t(mothers.size());
            auto pcm = parameters.at("pcm")[0];
            auto pvm = parameters.at("pvm")[0];
            for (auto i = 0u; i < mothers.size(); ++i) {
                birth_icm[i] = ica[mothers[i]] * pcm;
                birth_ivm[i] = iva[mothers[i]] * pvm;
            }

            api.clear_schedule("infection", died);
            api.clear_schedule("asymptomatic_infection", died);

            api.queue_variable_update("human", "birth", died, variable_vector_t{static_cast<double>(timestep)});
            api.queue_variable_update("human", "last_boosted_ib", died, variable_vector_t{-1});
            api.queue_variable_update("human", "last_boosted_ica", died, variable_vector_t{-1});
            api.queue_variable_update("human", "last_boosted_iva", died, variable_vector_t{-1});
            api.queue_variable_update("human", "last_boosted_id", died, variable_vector_t{-1});
            api.queue_variable_update("human", "IB", died, variable_vector_t{0});
            api.queue_variable_update("human", "ICA", died, variable_vector_t{0});
            api.queue_variable_update("human", "IVA", died, variable_vector_t{0});
            api.queue_variable_update("human", "ID", died, variable_vector_t{0});
            api.queue_variable_update("human", "ICM", died, birth_icm);
            api.queue_variable_update("human", "IVM", died, birth_ivm);
            api.queue_variable_update("human", "drug", died, variable_vector_t{0});
            api.queue_variable_update("human", "drug_time", died, variable_vector_t{-1});
            api.queue_variable_update("human", "infectivity", died, variable_vector_t{0});
            //zeta and zeta group survive rebirth
        }
    }), true);
}
