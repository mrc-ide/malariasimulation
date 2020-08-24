#include "test-helpers.h"
#include <vector>
#include <unordered_set>
#include <individual.h>

individual_index_t individual_index(size_t size, std::vector<size_t> i) {
    return individual_index_t(size, i.cbegin(), i.cend());
}
