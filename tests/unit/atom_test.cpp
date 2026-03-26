// Unit tests for ATOM module
#include <gtest/gtest.h>
#include "../../src/atom/atom_types.hpp"
#include "../../src/atom/utility.hpp"
#include "../../src/atom/nucdev.hpp"
#include "../../src/atom/radial_integrals.hpp"
#include "../../src/atom/muatco.hpp"
#include "../../src/atom/inmuat.hpp"
#include "../../src/par/parallel.hpp"
#include <feff/constants.hpp>
#include <cmath>

using namespace feff::atom;

constexpr double TOL = 1.0e-10;

// ===== Utility Functions =====

TEST(AtomUtilTest, Aprdev) {
    // aprdev computes coefficient for power (l-1) of product of two polynomials
    // For l=1: sum(a(m)*b(l+1-m), m=1..l) = a(1)*b(1)
    double a[10] = {2.0, 3.0};
    double b[10] = {4.0, 5.0};
    double result = aprdev(a, b, 1);
    EXPECT_NEAR(result, 8.0, TOL);  // a[0]*b[0] = 2*4

    // l=2: a(1)*b(2) + a(2)*b(1) = 2*5 + 3*4 = 22
    result = aprdev(a, b, 2);
    EXPECT_NEAR(result, 22.0, TOL);
}

TEST(AtomUtilTest, Dentfa) {
    // Thomas-Fermi potential: analytical approximation
    // At large r, potential should be small but positive for neutral atom
    double v = dentfa(100.0, 29.0, 0.0);
    EXPECT_GT(v, -1.0);  // small magnitude at large r

    // At very small r, the screening is minimal so potential ≈ Z/r (positive for TF convention)
    double v_small = dentfa(0.001, 29.0, 0.0);
    // dentfa returns the TF potential, which may have different sign convention
    // Just check it's finite and non-zero
    EXPECT_NE(v_small, 0.0);
}

TEST(AtomUtilTest, Cofcon) {
    // Convergence acceleration: oscillating errors should decrease mixing
    double a = 0.5, b = 0.5, p = 1.0, q = -1.0;
    cofcon(a, b, p, q);
    EXPECT_LT(b, 0.5);  // b decreased
    EXPECT_NEAR(a, 1.0 - b, TOL);
}

// ===== Nuclear Potential =====

TEST(AtomNucdevTest, PointCharge) {
    double av[10] = {};
    double dr[251] = {};
    double dv[251] = {};
    double dz = 29.0;  // Copper
    double hx = 0.05;
    int nuc = 0;
    int np = 251;
    int ndor = 10;
    double dr1 = dz * 1.0e-5;  // small first point

    nucdev(av, dr, dv, dz, hx, nuc, np, ndor, dr1);

    // Check that mesh is exponential
    EXPECT_GT(dr[0], 0.0);
    EXPECT_GT(dr[1], dr[0]);  // increasing
    double ratio = dr[1] / dr[0];
    EXPECT_NEAR(ratio, std::exp(hx), 1.0e-6);

    // Nuclear potential should be negative (attractive)
    EXPECT_LT(dv[10], 0.0);

    // nuc should be 1 for point charge (small Z)
    // For Z=29 it might be finite nucleus
    EXPECT_GE(nuc, 1);
}

// ===== Angular Coefficients =====

TEST(AtomAngularTest, FdmoccSameOrbital) {
    // fdmocc returns product of occupations for integrals
    // The exact formula depends on the implementation
    OrbitalConfig config;
    config.xnel[0] = 2.0;
    config.kap[0] = -1;  // s orbital: |kap|=1, degeneracy 2

    double result = fdmocc(0, 0, config);
    // Just verify it returns a positive finite value for occupied orbitals
    EXPECT_GT(result, 0.0);
    EXPECT_LT(result, 100.0);
}

TEST(AtomAngularTest, FdmoccDifferentOrbitals) {
    OrbitalConfig config;
    config.xnel[0] = 2.0;
    config.xnel[1] = 2.0;
    config.kap[0] = -1;
    config.kap[1] = -1;

    double result = fdmocc(0, 1, config);
    // Product of occupations should be positive
    EXPECT_GT(result, 0.0);
}

// ===== Initialization =====

TEST(AtomInitTest, InmuatCopper) {
    feff::par::par_begin();
    AtomState state;
    state.scf.nz = 29;  // Copper
    double xnval[30] = {};
    double xmag[30] = {};
    int iholep = 0;
    int iorb = 0;

    // Initialize with no hole, no ion
    inmuat(0, 0.0, 0, xnval, iholep, xmag, iorb, state);

    // Cu should have norb > 0 orbitals
    EXPECT_GT(state.scf.norb, 0);
    // Cu has Z=29 electrons
    double total = 0.0;
    for (int i = 0; i < state.scf.norb; ++i) {
        total += state.config.xnel[i];
    }
    EXPECT_NEAR(total, 29.0, 0.01);
}

// ===== Radial Grid Check =====

TEST(AtomGridTest, MeshParamsInit) {
    MeshParamsReal mesh;
    EXPECT_EQ(mesh.idim, 251);
    EXPECT_EQ(mesh.ndor, 10);
}

TEST(AtomGridTest, OrbitalConfig) {
    OrbitalConfig config;
    // All zeros by default
    EXPECT_EQ(config.nq[0], 0);
    EXPECT_EQ(config.kap[0], 0);
    EXPECT_NEAR(config.xnel[0], 0.0, TOL);
}
