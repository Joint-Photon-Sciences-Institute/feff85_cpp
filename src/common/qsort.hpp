#pragma once
// Index-based sorting utility.
// Replaces Fortran qsorti() (ACM Algorithm 402) with std::sort.

#include <vector>
#include <algorithm>
#include <numeric>

namespace feff::common {

/// Return indices that would sort the array in ascending order.
/// Equivalent to Fortran qsorti(ord, n, a) but returns 0-based indices.
/// Example: argsort({3.0, 1.0, 2.0}) returns {1, 2, 0}.
inline std::vector<int> argsort(const std::vector<double>& a) {
    std::vector<int> idx(a.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&a](int i, int j) { return a[i] < a[j]; });
    return idx;
}

/// Overload for raw array with count.
inline std::vector<int> argsort(const double* a, int n) {
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [a](int i, int j) { return a[i] < a[j]; });
    return idx;
}

} // namespace feff::common
