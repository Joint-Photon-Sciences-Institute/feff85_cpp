// Integration tests: verify that feff.inp parsing (RDINP) works for each
// material by reading the baseline feff.inp and checking key parameters.

#include "test_helpers.hpp"

#include <filesystem>
#include <string>

// RDINP module headers
#include "../../src/rdinp/iniall.hpp"
#include "../../src/rdinp/rdinp.hpp"
#include <feff/feff_input.hpp>

namespace fs = std::filesystem;

// ============================================================================
// Expected parameters for each material (from the POTENTIALS / EDGE cards)
// ============================================================================
struct MaterialExpected {
    std::string name;
    int nph;        // number of unique potentials (0..nph)
    int ihole;      // core-hole type: 1=K, 4=L3, etc.
    int iz0;        // atomic number of absorbing atom (pot 0)
    int min_natt;   // minimum number of atoms expected
};

static const MaterialExpected materials[] = {
    // name             nph  ihole  iz0  min_natt
    {"Copper",            1,    1,   29,     50},
    {"NiO",               2,    1,   28,    100},
    {"UO2",               2,    4,   92,     50},   // L3 edge -> ihole=4
    {"Zircon",            3,    1,   14,     50},
    {"ferrocene",         6,    1,   26,     50},
    {"bromoadamantane",   2,    1,   35,     20},
    {"LCO-para",          3,    1,   29,     50},
    {"LCO-perp",          3,    1,   29,     50},
};

// ============================================================================
// Parameterized fixture
// ============================================================================
class RdinpTest : public ::testing::TestWithParam<MaterialExpected> {
protected:
    fs::path original_cwd;
    std::string work_dir;

    void SetUp() override {
        auto& mat = GetParam();

        // Set up work dir mirroring Fortran test structure
        feff::test::setup_material_work_dir(mat.name);
        work_dir = feff::test::get_material_work_dir(mat.name);

        ASSERT_TRUE(fs::exists(work_dir + "/feff.inp"))
            << "feff.inp missing in: " << work_dir;

        original_cwd = fs::current_path();
        fs::current_path(work_dir);
    }

    void TearDown() override {
        fs::current_path(original_cwd);
        // Keep working dirs for inspection
    }
};

// ---------------------------------------------------------------------------
// Test: parse feff.inp and check key fields
// ---------------------------------------------------------------------------
TEST_P(RdinpTest, ParseFeffInp) {
    auto& expected = GetParam();

    // Initialize input structure with defaults
    feff::FeffInput inp;
    feff::rdinp::iniall(inp);

    // Parse feff.inp
    int nabs = feff::rdinp::rdinp(inp, "feff.inp");

    // nabs should be >= 1 (number of absorbers)
    EXPECT_GE(nabs, 1)
        << expected.name << ": rdinp returned nabs=" << nabs;

    // Check number of unique potentials
    EXPECT_EQ(inp.nph, expected.nph)
        << expected.name << ": nph mismatch";

    // Check core-hole type
    EXPECT_EQ(inp.ihole, expected.ihole)
        << expected.name << ": ihole mismatch"
        << " (1=K, 2=L1, 3=L2, 4=L3)";

    // Check atomic number of the absorbing atom (potential 0)
    EXPECT_EQ(inp.iz[0], expected.iz0)
        << expected.name << ": iz[0] (absorber Z) mismatch";

    // Check that atoms were read (natt = number of atoms in ATOMS list)
    EXPECT_GE(inp.natt, expected.min_natt)
        << expected.name << ": too few atoms parsed (natt=" << inp.natt << ")";

    // Sanity: rmax should be positive if set
    EXPECT_GT(inp.rmax, 0.0f)
        << expected.name << ": rmax should be positive";

    // Sanity: at least some potential labels should be non-empty
    EXPECT_FALSE(inp.potlbl[0].empty())
        << expected.name << ": potlbl[0] should not be empty";
}

// ---------------------------------------------------------------------------
// Register test cases
// ---------------------------------------------------------------------------
INSTANTIATE_TEST_SUITE_P(Materials, RdinpTest,
    ::testing::ValuesIn(materials),
    [](const ::testing::TestParamInfo<MaterialExpected>& info) {
        std::string name = info.param.name;
        std::replace(name.begin(), name.end(), '-', '_');
        return name;
    });
