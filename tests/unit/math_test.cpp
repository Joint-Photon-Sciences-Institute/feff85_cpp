#include <gtest/gtest.h>
#include <cmath>
#include <complex>

#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <feff/types.hpp>

#include "../../src/math/distance.hpp"
#include "../../src/math/bessel.hpp"
#include "../../src/math/legendre.hpp"
#include "../../src/math/wigner.hpp"
#include "../../src/math/integration.hpp"
#include "../../src/math/interpolation.hpp"
#include "../../src/math/phase_amplitude.hpp"
#include "../../src/math/polynomial_roots.hpp"
#include "../../src/math/determinant.hpp"

using namespace feff;
using namespace feff::math;

constexpr double TOL = 1.0e-12;

// ============== Distance tests ==============

TEST(DistanceTest, SamePoint) {
    double r0[] = {1.0, 2.0, 3.0};
    double r1[] = {1.0, 2.0, 3.0};
    EXPECT_NEAR(dist(r0, r1), 0.0, TOL);
}

TEST(DistanceTest, UnitDistance) {
    double r0[] = {0.0, 0.0, 0.0};
    double r1[] = {1.0, 0.0, 0.0};
    EXPECT_NEAR(dist(r0, r1), 1.0, TOL);
}

TEST(DistanceTest, Diagonal) {
    double r0[] = {0.0, 0.0, 0.0};
    double r1[] = {1.0, 1.0, 1.0};
    EXPECT_NEAR(dist(r0, r1), std::sqrt(3.0), TOL);
}

TEST(DistanceTest, SinglePrecision) {
    float r0[] = {0.0f, 0.0f, 0.0f};
    float r1[] = {3.0f, 4.0f, 0.0f};
    EXPECT_NEAR(sdist(r0, r1), 5.0f, 1.0e-6f);
}

// ============== Bessel function tests ==============

TEST(BesselTest, BjnserL0) {
    FeffComplex x(2.0, 0.0);
    FeffComplex jl, nl;
    bjnser(x, 0, jl, nl, 0);

    // j0(x) = sin(x)/x
    FeffComplex expected_jl = std::sin(x) / x;
    EXPECT_NEAR(std::abs(jl - expected_jl), 0.0, 1.0e-10);

    // n0(x) = -cos(x)/x
    FeffComplex expected_nl = -std::cos(x) / x;
    EXPECT_NEAR(std::abs(nl - expected_nl), 0.0, 1.0e-10);
}

TEST(BesselTest, ExjlnlL0) {
    FeffComplex z(5.0, 0.0);
    FeffComplex jl, nl;
    exjlnl(z, 0, jl, nl);

    FeffComplex expected_jl = std::sin(z) / z;
    EXPECT_NEAR(std::abs(jl - expected_jl), 0.0, 1.0e-10);
}

TEST(BesselTest, ExjlnlL1) {
    FeffComplex z(3.0, 0.0);
    FeffComplex jl, nl;
    exjlnl(z, 1, jl, nl);

    FeffComplex expected_jl = std::sin(z) / (z * z) - std::cos(z) / z;
    EXPECT_NEAR(std::abs(jl - expected_jl), 0.0, 1.0e-10);
}

TEST(BesselTest, BesjnComplex) {
    FeffComplex x(5.0, 0.5);
    FeffComplex jl[feff::ltot + 2], nl[feff::ltot + 2];
    besjn(x, jl, nl);

    // j0(x) = sin(x)/x
    FeffComplex expected_j0 = std::sin(x) / x;
    EXPECT_NEAR(std::abs(jl[0] - expected_j0), 0.0, 1.0e-8);
}

// ============== Legendre polynomial tests ==============

TEST(LegendreTest, P0) {
    double pl0[5];
    cpl0(0.5, pl0, 5);
    EXPECT_NEAR(pl0[0], 1.0, TOL);
}

TEST(LegendreTest, P1) {
    double pl0[5];
    cpl0(0.5, pl0, 5);
    EXPECT_NEAR(pl0[1], 0.5, TOL);
}

TEST(LegendreTest, P2) {
    double pl0[5];
    double x = 0.5;
    cpl0(x, pl0, 5);
    // P2(x) = (3x^2 - 1)/2
    double expected = (3.0 * x * x - 1.0) / 2.0;
    EXPECT_NEAR(pl0[2], expected, TOL);
}

// ============== Wigner 3j tests ==============

TEST(WignerTest, Cwig3jBasic) {
    // 3j(1, 1, 0; 0, 0, ient=1) should be nonzero
    double val = cwig3j(1, 1, 0, 0, 0, 1);
    // Known value: (-1)^(1) / sqrt(3)
    EXPECT_NEAR(std::abs(val), 1.0 / std::sqrt(3.0), 1.0e-10);
}

TEST(WignerTest, Cwig3jTriangleViolation) {
    // j1+j2 < j3 should give zero
    double val = cwig3j(1, 1, 5, 0, 0, 1);
    EXPECT_NEAR(val, 0.0, TOL);
}

// ============== Integration tests ==============

TEST(IntegrationTest, TrapLinear) {
    // Integrate y = x from 0 to 1
    const int n = 101;
    double x[n], y[n];
    for (int i = 0; i < n; ++i) {
        x[i] = static_cast<double>(i) / (n - 1);
        y[i] = x[i];
    }
    double sum;
    trap(x, y, n, sum);
    EXPECT_NEAR(sum, 0.5, 1.0e-6);
}

TEST(IntegrationTest, TrapQuadratic) {
    // Integrate y = x^2 from 0 to 1
    const int n = 1001;
    double x[n], y[n];
    for (int i = 0; i < n; ++i) {
        x[i] = static_cast<double>(i) / (n - 1);
        y[i] = x[i] * x[i];
    }
    double sum;
    trap(x, y, n, sum);
    EXPECT_NEAR(sum, 1.0 / 3.0, 1.0e-5);
}

// ============== Interpolation tests ==============

TEST(InterpolationTest, LocatMiddle) {
    double xx[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    int idx = locat(2.5, 5, xx);
    EXPECT_EQ(idx, 2);  // xx[1] = 2.0 <= 2.5 < xx[2] = 3.0
}

TEST(InterpolationTest, TerpLinear) {
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[] = {0.0, 1.0, 4.0, 9.0, 16.0};
    double y0;
    terp(x, y, 5, 2, 1.5, y0);
    // Quadratic interpolation of x^2 should be exact
    EXPECT_NEAR(y0, 2.25, 1.0e-10);
}

// ============== Determinant tests ==============

TEST(DeterminantTest, Identity2x2) {
    double mat[] = {1.0, 0.0, 0.0, 1.0};
    double det = determ(mat, 2, 2);
    EXPECT_NEAR(det, 1.0, TOL);
}

TEST(DeterminantTest, Known2x2) {
    // |1 2| = 1*4 - 2*3 = -2
    // |3 4|
    double mat[] = {1.0, 2.0, 3.0, 4.0};
    double det = determ(mat, 2, 2);
    EXPECT_NEAR(det, -2.0, TOL);
}

// ============== Polynomial roots tests ==============

TEST(PolynomialTest, QuadraticRoots) {
    // x^2 - 3x + 2 = 0 => x = 1, 2
    FeffComplex coef[] = {1.0, -3.0, 2.0};
    FeffComplex sol[2];
    int nsol;
    cqdrtc(coef, sol, nsol);
    EXPECT_EQ(nsol, 2);
    // Sort by real part
    if (sol[0].real() > sol[1].real()) std::swap(sol[0], sol[1]);
    EXPECT_NEAR(sol[0].real(), 1.0, 1.0e-10);
    EXPECT_NEAR(sol[1].real(), 2.0, 1.0e-10);
}

// ============== Phase amplitude tests ==============

TEST(PhaseAmplitudeTest, AtanccReal) {
    // atan(1) = pi/4
    FeffComplex temp(1.0, 0.0);
    FeffComplex phx;
    atancc(temp, phx);
    EXPECT_NEAR(phx.real(), feff::pi / 4.0, 1.0e-10);
    EXPECT_NEAR(phx.imag(), 0.0, 1.0e-10);
}
