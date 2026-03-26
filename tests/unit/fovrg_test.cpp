// Unit tests for FOVRG module
#include <gtest/gtest.h>
#include "../../src/atom/atom_types.hpp"
#include "../../src/fovrg/nucdec.hpp"
#include "../../src/fovrg/radial_integrals_c.hpp"
#include <feff/constants.hpp>
#include <feff/dimensions.hpp>
#include <cmath>
#include <complex>

using namespace feff;
using namespace feff::atom;

constexpr double TOL = 1.0e-10;

// ===== Nuclear Potential (Complex version) =====

TEST(FovrgNucdecTest, PointCharge) {
    double av[10] = {};
    double dr[nrptx] = {};
    double dv[nrptx] = {};
    double dz = 29.0;  // Copper
    double hx = 0.05;
    int nuc = 0;
    int np = nrptx;
    int ndor = 10;
    double dr1 = dz * 1.0e-5;

    feff::fovrg::nucdec(av, dr, dv, dz, hx, nuc, np, ndor, dr1);

    // Mesh should be exponential
    EXPECT_GT(dr[0], 0.0);
    EXPECT_GT(dr[1], dr[0]);
    double ratio = dr[1] / dr[0];
    EXPECT_NEAR(ratio, std::exp(hx), 1.0e-6);

    // Nuclear potential should be negative
    EXPECT_LT(dv[10], 0.0);
}

// ===== Complex Polynomial Product =====

TEST(FovrgRadialTest, Aprdec) {
    FeffComplex a[10] = {{2.0, 1.0}, {3.0, 0.5}};
    double b[10] = {4.0, 5.0};

    // l=1: a[0]*b[0] = (2+i)*(4) = (8+4i)
    FeffComplex result = feff::fovrg::aprdec(a, b, 1);
    EXPECT_NEAR(result.real(), 8.0, TOL);
    EXPECT_NEAR(result.imag(), 4.0, TOL);

    // l=2: a[0]*b[1] + a[1]*b[0] = (2+i)*5 + (3+0.5i)*4 = (10+5i) + (12+2i) = (22+7i)
    result = feff::fovrg::aprdec(a, b, 2);
    EXPECT_NEAR(result.real(), 22.0, TOL);
    EXPECT_NEAR(result.imag(), 7.0, TOL);
}

TEST(FovrgRadialTest, AprdepReal) {
    double a[10] = {2.0, 3.0};
    double b[10] = {4.0, 5.0};

    double result = feff::fovrg::aprdep(a, b, 1);
    EXPECT_NEAR(result, 8.0, TOL);

    result = feff::fovrg::aprdep(a, b, 2);
    EXPECT_NEAR(result, 22.0, TOL);
}

// ===== Struct Initialization =====

TEST(FovrgTypesTest, FovrgStateInit) {
    FovrgState state;
    EXPECT_EQ(state.scf.norb, 0);
    EXPECT_EQ(state.scf.nz, 0);
    EXPECT_NEAR(state.work.cl, 0.0, TOL);
    EXPECT_EQ(state.mesh.idim, nrptx);
}

TEST(FovrgTypesTest, DiracWorkspaceComplex) {
    DiracWorkspaceComplex work;
    // All zeros by default
    EXPECT_NEAR(work.gg[0].real(), 0.0, TOL);
    EXPECT_NEAR(work.gg[0].imag(), 0.0, TOL);
    EXPECT_NEAR(work.cl, 0.0, TOL);
}

TEST(FovrgTypesTest, AngularCoefficientsC) {
    AngularCoefficientsC ang;
    // Test offset indexing: (-ltot-1..ltot, 30, 0..3)
    ang(-5, 0, 0) = 3.14;
    EXPECT_NEAR(ang(-5, 0, 0), 3.14, TOL);

    ang(ltot, 0, 0) = 2.71;
    EXPECT_NEAR(ang(ltot, 0, 0), 2.71, TOL);
}
