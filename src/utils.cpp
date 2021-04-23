#include "utils.h"

//[[Rcpp::export]]
std::vector<size_t> bernoulli_multi_p_cpp(const std::vector<double> p) {
    auto values = Random::get_instance().bernoulli_multi_p(p);
    for (auto i = 0u; i < values.size(); ++i) {
        values[i]++;
    }
    return values;
}

//[[Rcpp::export]]
Rcpp::IntegerVector simulate_competing_outcomes(Rcpp::NumericMatrix p) {
    return simulate_competing_outcomes(&Random::get_instance(), p);
}

Rcpp::IntegerVector simulate_competing_outcomes(
    RandomInterface* random,
    Rcpp::NumericMatrix p
    ) {
    auto nrow = p.nrow();
    auto ncol = p.ncol();
    std::vector<double> prob_occurs(nrow);
	for (auto i = 0; i < nrow; ++i) {
		double row_product = 1.;
		for (auto j = 0; j < ncol; ++j) {
			row_product *= 1 - p(i,j);
		}
		prob_occurs[i] = 1 - row_product;
	}

    auto rand = random->runif(nrow);

    Rcpp::IntegerVector outcomes(nrow);
    outcomes.fill(NA_INTEGER);
    for (auto i = 0; i < nrow; ++i) {
        if (rand[i] < prob_occurs[i]) {
            outcomes[i] = random->sample_int(
                ncol,
                Rcpp::as<std::vector<double>>(Rcpp::NumericVector(p.row(i)))
            );
        }
    }
    return outcomes;
}
