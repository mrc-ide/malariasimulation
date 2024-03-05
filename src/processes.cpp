#include <Rcpp.h>
#include <individual.h>

/**
 * An iterator adaptor which yields the same values as the underlying iterator,
 * but scaled by a pre-determined factor.
 *
 * This is used by the exponential_process below to scale an std::vector by a
 * constant.
 *
 * There are two straightforward ways of performing the operation. The first is
 * to create an empty vector, use `reserve(N)` to pre-allocate the vector and
 * then call `push_back` with each new value. The second way would be to create
 * a zero-initialised vector of size N and then use `operator[]` to fill in the
 * values.
 *
 * Unfortunately both approaches have significant overhead. In the former, the
 * use of `push_back` requires repeated checks as to whether the vector needs
 * growing, despite the prior reserve call. These calls inhibits optimizations
 * such as unrolling and auto-vectorization of the loop. The latter approach
 * requires an initial memset when zero-initializing the vector, even though the
 * vector then gets overwritten entirely. Sadly gcc fails to optimize out either
 * of those. Ideally we want a way to create a pre-sized but uninitialised
 * std::vector we can write to ourselves, but there is no API in the standard
 * library to do this. All existing workarounds end up with an std::vector with
 * non-default item type or allocators.
 *
 * There is however a way out! std::vector has a constructor which accepts a
 * pair of iterators and fills the vector with values from the iterators. Using
 * `std::distance` on the iterator pair it can even pre-allocate the vector to
 * the right size. No zero-initialisation or no capacity checks, just one
 * allocation and a straightforward easily optimizable loop. All we need is an
 * iterator yielding the right values, hence `scale_iterator`. In C++20 we would
 * probably be able to use the new ranges library as our iterators.
 *
 * How much does this matter? On microbenchmarks, for small and medium sized
 * vector (<= 1M doubles), this version is about 30% faster than the
 * zero-initialising implementation and 60% faster than the one which uses
 * push_back. For larger vector sizes the difference is less pronounced,
 * possibly because caches become saturated. At the time of writing, on a
 * real-word run of malariasimulation with a population size of 1M the overall
 * speedup is about 2-3%.
 *
 * https://wolchok.org/posts/cxx-trap-1-constant-size-vector/
 * https://codingnest.com/the-little-things-the-missing-performance-in-std-vector/
 * https://lemire.me/blog/2012/06/20/do-not-waste-time-with-stl-vectors/
 */
template<typename underlying_iterator>
struct scale_iterator  {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = typename std::iterator_traits<underlying_iterator>::difference_type;
    using value_type = typename std::iterator_traits<underlying_iterator>::value_type;
    using pointer = typename std::iterator_traits<underlying_iterator>::pointer;

    // We skirt the rules a bit by returning a prvalue from `operator*`, even
    // though the C++17 (and prior) standard says forward iterators are supposed
    // to return a reference type (ie. a glvalue). Because the scaling is
    // applied on the fly, there is no glvalue we could return a reference to.
    //
    // An input iterator would be allowed to return a prvalue, but the
    // std::vector constructor wouldn't be able to figure out the length ahead
    // of time if we were an input iterator.
    //
    // C++20 actually introduces parallel definitions of input and forward
    // iterators, which relax this requirement, so under that classification our
    // implementation in correct.
    //
    // In practice though, this does not really matter. We only use this
    // iterator in one specific context, and the vector constructor doesn't do
    // anything elaborate that we would be upsetting.
    using reference = value_type;

    scale_iterator(underlying_iterator it, value_type factor) : it(it), factor(factor) {}
    reference operator*() {
        return factor * (*it);
    }
    bool operator==(const scale_iterator& other) {
        return it == other.it;
    }
    bool operator!=(const scale_iterator& other) {
        return it != other.it;
    }
    scale_iterator& operator++() {
        it++;
        return *this;
    }
    scale_iterator operator++(int) {
        return scale_iterator(it++, factor);
    }

    private:
    underlying_iterator it;
    value_type factor;
};

template<typename T>
scale_iterator<T> make_scale_iterator(T&& it, typename std::iterator_traits<T>::value_type scale) {
    return scale_iterator<T>(std::forward<T>(it), scale);
}

//[[Rcpp::export]]
Rcpp::XPtr<process_t> exponential_process_cpp(
    Rcpp::XPtr<DoubleVariable> variable,
    const double rate
){
    return Rcpp::XPtr<process_t>(
        new process_t([=](size_t t){
            const std::vector<double>& values = variable->get_values();
            std::vector<double> new_values(
              make_scale_iterator(values.cbegin(), rate),
              make_scale_iterator(values.cend(), rate));

            variable->queue_update(std::move(new_values), std::vector<size_t>());
        }),
        true
    );
}
