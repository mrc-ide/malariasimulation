
#include <Rcpp.h>
#include "Random.h"

//[[Rcpp::export]]
Rcpp::XPtr<individual_index_t> bernoulli_multi_p_cpp(const std::vector<double> p) {
    return Rcpp::XPtr<individual_index_t>(
        new individual_index_t(Random::get_instance().bernoulli_multi_p(p)),
        true
    );
}
