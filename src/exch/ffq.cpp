#include "ffq.hpp"
#include <cmath>

namespace feff::exch {

double ffq(double q, double ef, double xk, double wp, double alph) {
    double wq = std::sqrt(wp * wp + alph * q * q + q * q * q * q);
    double result = (wp + wq) / (q * q) + alph / (2.0 * wp);
    result = ((ef * wp) / (4.0 * xk)) * std::log(result);
    return result;
}

} // namespace feff::exch
