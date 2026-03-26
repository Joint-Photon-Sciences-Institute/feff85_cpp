// Unit tests for EXCH module
#include <gtest/gtest.h>
#include "../../src/exch/vbh.hpp"
#include "../../src/exch/edp.hpp"
#include "../../src/exch/rhl.hpp"
#include "../../src/exch/quinn.hpp"
#include "../../src/exch/ffq.hpp"
#include "../../src/exch/cubic_exch.hpp"
#include "../../src/exch/imhl.hpp"
#include "../../src/exch/fndsng.hpp"
#include "../../src/exch/csigma.hpp"
#include <feff/constants.hpp>
#include <cmath>

using namespace feff;
using namespace feff::exch;

constexpr double TOL = 1.0e-8;

// ===== Von Barth-Hedin =====

TEST(ExchVbhTest, TypicalMetal) {
    // rs ~ 2 for typical metal, unpolarized (xmag=1)
    double vxc;
    vbh(2.0, 1.0, vxc);
    // VBH should give negative xc potential (attractive)
    EXPECT_LT(vxc, 0.0);
    // Typical range: -0.5 to -0.1 Hartrees
    EXPECT_GT(vxc, -1.0);
}

TEST(ExchVbhTest, LargeRs) {
    // Large rs = low density → small vxc
    double vxc;
    vbh(10.0, 1.0, vxc);
    EXPECT_LT(vxc, 0.0);
    // Should be smaller magnitude than rs=2
    double vxc2;
    vbh(2.0, 1.0, vxc2);
    EXPECT_GT(vxc, vxc2);  // less negative (smaller magnitude)
}

TEST(ExchVbhTest, SmallRs) {
    // Small rs = high density → large |vxc|
    double vxc;
    vbh(0.5, 1.0, vxc);
    EXPECT_LT(vxc, -0.5);
}

// ===== Dirac-Hara =====

TEST(ExchEdpTest, AtFermiLevel) {
    // At Fermi level, xk = xf = fa/rs
    double rs = 2.0;
    double xf = feff::fa / rs;
    double vr;
    edp(rs, xf, vr);
    // Should give exchange potential (negative)
    EXPECT_LT(vr, 0.0);
}

TEST(ExchEdpTest, ZeroMomentum) {
    // At zero momentum
    double vr;
    edp(2.0, 0.001, vr);
    // Should give finite negative value
    EXPECT_LT(vr, 0.0);
    EXPECT_GT(vr, -10.0);
}

// ===== Hedin-Lundqvist =====

TEST(ExchRhlTest, BelowFermi) {
    // Below Fermi level (xk < kf)
    double erl, eim;
    rhl(2.0, 0.5, erl, eim);
    // Real part should be negative (correlation hole)
    EXPECT_LT(erl, 0.0);
    // Imaginary part should be zero or very small below Fermi
    EXPECT_NEAR(eim, 0.0, 0.1);
}

TEST(ExchRhlTest, AboveFermi) {
    // Well above Fermi level (xk > kf)
    double erl, eim;
    rhl(2.0, 2.0, erl, eim);
    // Imaginary part should be negative (damping)
    EXPECT_LE(eim, 0.0);
}

// ===== Quinn =====

TEST(ExchQuinnTest, BelowThreshold) {
    // Below plasmon threshold, ei should be zero or small
    double ei;
    double rs = 2.0;
    double wp = std::sqrt(3.0 / (rs * rs * rs));
    double ef = 0.5 * feff::fa * feff::fa / (rs * rs);
    quinn(0.5, rs, wp, ef, ei);
    // ei should be non-positive (loss)
    EXPECT_LE(ei, 0.001);
}

// ===== FFQ =====

TEST(ExchFfqTest, BasicValue) {
    double result = ffq(1.0, 1.0, 1.0, 1.0, 4.0 / 3.0);
    // Just check it returns finite value
    EXPECT_TRUE(std::isfinite(result));
}

// ===== Cubic Solver =====

TEST(ExchCubicTest, BasicSolve) {
    double rad, qplus, qminus;
    cubic_exch(1.0, 1.0, 4.0 / 3.0, rad, qplus, qminus);
    // Check outputs are finite
    EXPECT_TRUE(std::isfinite(rad));
    EXPECT_TRUE(std::isfinite(qplus));
    EXPECT_TRUE(std::isfinite(qminus));
}

// ===== Imaginary HL =====

TEST(ExchImhlTest, AtFermi) {
    double eim;
    int icusp;
    // At Fermi level, imaginary part should be zero
    imhl(2.0, 1.0, eim, icusp);
    // eim near Fermi should be small
    EXPECT_TRUE(std::isfinite(eim));
}

// ===== Singularity Finder =====

TEST(ExchFndsngTest, BasicCall) {
    FeffComplex limit1(0.0, 0.0);
    FeffComplex limit2(3.0, 0.0);
    int nsing = 0;
    FeffComplex xsing[20] = {};
    double dppar[10] = {1.0, 1.0, 1.0};  // kF, wp, energy
    FeffComplex cpar[10] = {{1.0, 0.0}};
    int ifcn = 1;

    fndsng(limit1, limit2, nsing, xsing, dppar, cpar, ifcn);
    // Should find 0 or more singularities
    EXPECT_GE(nsing, 0);
    EXPECT_LE(nsing, 20);
}

// ===== CSigma internal functions =====

TEST(ExchCsigmaTest, HFExchange) {
    FeffComplex ck(1.0, 0.0);  // at Fermi momentum
    double ef = 1.0;
    double kf = 1.0;
    FeffComplex result = hfexc(ck, ef, kf);
    // HF exchange should be real and negative
    EXPECT_LT(result.real(), 0.0);
    EXPECT_NEAR(result.imag(), 0.0, 0.01);
}
