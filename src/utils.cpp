
#include <Rcpp.h>
#include "Random.h"
#include <individual.h>

//[[Rcpp::export]]
void random_seed(size_t seed) {
    Random::get_instance().seed(seed);
}

//[[Rcpp::export]]
std::string random_save_state() {
    return Random::get_instance().save_state();
}

//[[Rcpp::export]]
void random_restore_state(std::string state) {
    return Random::get_instance().restore_state(state);
}

//[[Rcpp::export]]
std::vector<size_t> bernoulli_multi_p_cpp(const std::vector<double> p) {
    auto values = Random::get_instance().bernoulli_multi_p(p);
    for (auto i = 0u; i < values.size(); ++i) {
        values[i]++;
    }
    return values;
}

//[[Rcpp::export]]
std::vector<size_t> bitset_index_cpp(
    Rcpp::XPtr<individual_index_t> a,
    Rcpp::XPtr<individual_index_t> b
    ) {
    if (a->max_size() != b->max_size()) {
        Rcpp::stop("Incompatible bitmap sizes, %d vs %d",
                   a->max_size(), b->max_size());
    }

    auto values = std::vector<size_t>();
    auto i = 1u;
    for (const auto& v : *a) {
        if (b->find(v) != b->cend()) {
            values.push_back(i);
        }
        ++i;
    }
    return values;
}

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector fast_weighted_sample(
    size_t size,
    std::vector<double> probs
    ) {
    Rcpp::IntegerVector values(size);
    Random::get_instance().prop_sample_bucket(
        size,
        probs,
        INTEGER(values)
    );
    values = values + 1;
    return values;
}

//[[Rcpp::export]]
std::string individual_headers_version() {
    return INDIVIDUAL_VERSION;
}
