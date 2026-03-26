// Unit tests for COMMON module
#include <gtest/gtest.h>
#include "../../src/common/string_utils.hpp"
#include "../../src/common/periodic_table.hpp"
#include "../../src/common/radial_grid.hpp"
#include "../../src/common/physics_utils.hpp"
#include "../../src/common/qsort.hpp"
#include "../../src/common/pad_io.hpp"
#include "../../src/common/itoken.hpp"
#include "../../src/common/logging.hpp"
#include "../../src/par/parallel.hpp"
#include <feff/constants.hpp>
#include <cmath>
#include <sstream>

using namespace feff;
using namespace feff::common;

constexpr double TOL = 1.0e-12;

// ===== String Utilities =====

TEST(StringUtilsTest, Trim) {
    EXPECT_EQ(ltrim("  hello"), "hello");
    EXPECT_EQ(ltrim("\thello"), "hello");
    EXPECT_EQ(rtrim("hello  "), "hello");
    EXPECT_EQ(trim("  hello  "), "hello");
    EXPECT_EQ(trim(""), "");
}

TEST(StringUtilsTest, Case) {
    EXPECT_EQ(to_upper("hello"), "HELLO");
    EXPECT_EQ(to_lower("HELLO"), "hello");
    EXPECT_EQ(to_upper("Hello123"), "HELLO123");
}

TEST(StringUtilsTest, SplitWords) {
    auto words = split_words("  one  two  three  ");
    ASSERT_EQ(words.size(), 3u);
    EXPECT_EQ(words[0], "one");
    EXPECT_EQ(words[1], "two");
    EXPECT_EQ(words[2], "three");

    // Comma-separated
    words = split_words("a, b,c");
    ASSERT_EQ(words.size(), 3u);
    EXPECT_EQ(words[0], "a");
    EXPECT_EQ(words[1], "b");
    EXPECT_EQ(words[2], "c");
}

TEST(StringUtilsTest, IsComment) {
    EXPECT_TRUE(is_comment(""));
    EXPECT_TRUE(is_comment("   "));
    EXPECT_TRUE(is_comment("* comment"));
    EXPECT_TRUE(is_comment("# comment"));
    EXPECT_TRUE(is_comment("; comment"));
    EXPECT_TRUE(is_comment("% comment"));
    EXPECT_FALSE(is_comment("not a comment"));
}

TEST(StringUtilsTest, IsNumeric) {
    EXPECT_TRUE(is_numeric("123"));
    EXPECT_TRUE(is_numeric("-3.14"));
    EXPECT_TRUE(is_numeric("1.5e10"));
    EXPECT_TRUE(is_numeric("2.5d-3"));
    EXPECT_FALSE(is_numeric("abc"));
    EXPECT_FALSE(is_numeric("1e2e3"));  // two exponents
    EXPECT_FALSE(is_numeric(""));
}

TEST(StringUtilsTest, ParseDouble) {
    double val;
    EXPECT_TRUE(parse_double("3.14", val));
    EXPECT_NEAR(val, 3.14, TOL);

    EXPECT_TRUE(parse_double("2.5d-3", val));
    EXPECT_NEAR(val, 0.0025, TOL);

    EXPECT_FALSE(parse_double("abc", val));
}

TEST(StringUtilsTest, ParseInt) {
    int val;
    EXPECT_TRUE(parse_int("42", val));
    EXPECT_EQ(val, 42);

    EXPECT_TRUE(parse_int("-7", val));
    EXPECT_EQ(val, -7);

    EXPECT_FALSE(parse_int("3.14", val));
}

// ===== Periodic Table =====

TEST(PeriodicTableTest, Iron) {
    EXPECT_NEAR(atomic_weight(26), 55.85, 0.01);
    EXPECT_STREQ(atomic_symbol(26), "Fe");
}

TEST(PeriodicTableTest, Copper) {
    EXPECT_NEAR(atomic_weight(29), 63.55, 0.01);
    EXPECT_STREQ(atomic_symbol(29), "Cu");
}

TEST(PeriodicTableTest, Hydrogen) {
    EXPECT_NEAR(atomic_weight(1), 1.0079, 0.001);
    EXPECT_STREQ(atomic_symbol(1), "H");
}

TEST(PeriodicTableTest, Lookup) {
    EXPECT_EQ(atomic_number("Fe"), 26);
    EXPECT_EQ(atomic_number("fe"), 26);
    EXPECT_EQ(atomic_number("Cu"), 29);
    EXPECT_EQ(atomic_number("Xx"), 0);  // not found
}

TEST(PeriodicTableTest, BoundsCheck) {
    EXPECT_THROW(atomic_weight(0), std::out_of_range);
    EXPECT_THROW(atomic_weight(104), std::out_of_range);
}

// ===== Radial Grid =====

TEST(RadialGridTest, GridValues) {
    // xx(1) = -8.8
    EXPECT_NEAR(xx(1), -8.8, TOL);
    // xx(2) = -8.8 + 0.05 = -8.75
    EXPECT_NEAR(xx(2), -8.75, TOL);
    // rr(1) = exp(-8.8)
    EXPECT_NEAR(rr(1), std::exp(-8.8), TOL);
}

TEST(RadialGridTest, RoundTrip) {
    // ii(rr(j)) should give j (or j-1 due to floor)
    for (int j = 10; j <= 200; ++j) {
        double r = rr(j);
        int idx = ii(r);
        EXPECT_TRUE(idx == j || idx == j - 1)
            << "j=" << j << " r=" << r << " ii(r)=" << idx;
    }
}

// ===== Physics Utils =====

TEST(PhysicsTest, GetXK) {
    // k = sqrt(2*0.5) = 1.0 for E > 0
    EXPECT_NEAR(getxk(0.5), 1.0, TOL);
    // k = -sqrt(2*0.5) = -1.0 for E < 0
    EXPECT_NEAR(getxk(-0.5), -1.0, TOL);
    // k = 0 for E = 0
    EXPECT_NEAR(getxk(0.0), 0.0, TOL);
}

TEST(PhysicsTest, Setkap) {
    int kinit, linit;

    // K edge (1s): l=0, kappa=-1
    setkap(1, kinit, linit);
    EXPECT_EQ(linit, 0);
    EXPECT_EQ(kinit, -1);

    // LIII edge (2p3/2): l=1, kappa=-2
    setkap(4, kinit, linit);
    EXPECT_EQ(linit, 1);
    EXPECT_EQ(kinit, -2);

    // LII edge (2p1/2): l=1, kappa=1
    setkap(3, kinit, linit);
    EXPECT_EQ(linit, 1);
    EXPECT_EQ(kinit, 1);
}

TEST(PhysicsTest, Pijump) {
    double ph = 3.0;
    double old = 3.0 + 2.0 * feff::pi;  // old is 2pi ahead
    pijump(ph, old);
    // ph should be adjusted to be close to old
    EXPECT_NEAR(ph - old, 0.0, 0.1);
}

TEST(PhysicsTest, SetgamCuK) {
    // Cu (Z=29) K-edge core-hole width: approximately 1.5 eV
    feff::par::par_begin();  // needed for logger
    double gamach;
    setgam(29, 1, gamach);
    EXPECT_GT(gamach, 0.5);
    EXPECT_LT(gamach, 5.0);
}

// ===== Quicksort =====

TEST(QsortTest, Argsort) {
    std::vector<double> a = {3.0, 1.0, 2.0};
    auto idx = argsort(a);
    ASSERT_EQ(idx.size(), 3u);
    EXPECT_EQ(idx[0], 1);  // smallest is a[1]=1.0
    EXPECT_EQ(idx[1], 2);  // next is a[2]=2.0
    EXPECT_EQ(idx[2], 0);  // largest is a[0]=3.0
}

TEST(QsortTest, ArgsortAlreadySorted) {
    std::vector<double> a = {1.0, 2.0, 3.0};
    auto idx = argsort(a);
    EXPECT_EQ(idx[0], 0);
    EXPECT_EQ(idx[1], 1);
    EXPECT_EQ(idx[2], 2);
}

// ===== PAD I/O =====

TEST(PadIoTest, RoundTrip) {
    // Encode and decode a few values
    constexpr int npack = 8;
    std::vector<double> test_values = {0.0, 1.0, -1.0, 3.14159265358979,
                                        1.0e-20, -2.5e15, 1.0e37};
    for (double val : test_values) {
        std::string encoded;
        pad_encode(val, npack, encoded);
        double decoded = pad_decode(encoded, npack);
        if (val == 0.0) {
            EXPECT_NEAR(decoded, 0.0, 1.0e-15) << "val=" << val;
        } else {
            double rel_err = std::abs((decoded - val) / val);
            EXPECT_LT(rel_err, 1.0e-10) << "val=" << val << " decoded=" << decoded;
        }
    }
}

TEST(PadIoTest, WriteReadDouble) {
    constexpr int npack = 8;
    double arr[] = {1.5, -2.7, 3.14, 0.001, 1e20};
    int n = 5;

    // Write to stringstream
    std::ostringstream oss;
    write_pad_double(oss, npack, arr, n);

    // Read back
    std::istringstream iss(oss.str());
    double result[5];
    int nread = read_pad_double(iss, npack, result, n);
    EXPECT_EQ(nread, n);

    for (int i = 0; i < n; ++i) {
        if (arr[i] == 0.0) {
            EXPECT_NEAR(result[i], 0.0, 1e-15);
        } else {
            double rel = std::abs((result[i] - arr[i]) / arr[i]);
            EXPECT_LT(rel, 1e-10) << "i=" << i;
        }
    }
}

TEST(PadIoTest, WriteReadComplex) {
    constexpr int npack = 8;
    FeffComplex arr[] = {{1.5, -0.3}, {-2.7, 4.1}, {0.0, 0.0}};
    int n = 3;

    std::ostringstream oss;
    write_pad_complex(oss, npack, arr, n);

    std::istringstream iss(oss.str());
    FeffComplex result[3];
    int nread = read_pad_complex(iss, npack, result, n);
    EXPECT_EQ(nread, n);

    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(result[i].real(), arr[i].real(), std::abs(arr[i].real()) * 1e-10 + 1e-15);
        EXPECT_NEAR(result[i].imag(), arr[i].imag(), std::abs(arr[i].imag()) * 1e-10 + 1e-15);
    }
}

// ===== IToken =====

TEST(ITokenTest, FeffTokens) {
    // Token IDs match Fortran itoken.f exactly
    EXPECT_EQ(itoken("ATOMS", "feff.inp"), 1);   // ATOM → 1
    EXPECT_EQ(itoken("HOLE", "feff.inp"), 2);
    EXPECT_EQ(itoken("POTENTIALS", "feff.inp"), 14);  // POTE → 14
    EXPECT_EQ(itoken("EDGE", "feff.inp"), 35);
    EXPECT_EQ(itoken("TITLE", "feff.inp"), 7);   // TITL → 7
    EXPECT_EQ(itoken("EXCH", "feff.inp"), 5);
}

TEST(ITokenTest, CaseInsensitive) {
    EXPECT_EQ(itoken("atoms", "feff.inp"), 1);
    EXPECT_EQ(itoken("Hole", "feff.inp"), 2);
}

TEST(ITokenTest, Unknown) {
    // Fortran returns 0 for unknown tokens
    EXPECT_EQ(itoken("ZZZZ", "feff.inp"), 0);
}
