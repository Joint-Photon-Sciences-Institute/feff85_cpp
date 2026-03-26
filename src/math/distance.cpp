#include "distance.hpp"
#include <cmath>

namespace feff::math {

double dist(const double r0[3], const double r1[3]) {
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
        double d = r0[i] - r1[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

float sdist(const float r0[3], const float r1[3]) {
    float sum = 0.0f;
    for (int i = 0; i < 3; ++i) {
        float d = r0[i] - r1[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

} // namespace feff::math
