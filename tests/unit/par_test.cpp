// Unit tests for PAR module (sequential mode)
#include <gtest/gtest.h>
#include "../../src/par/parallel.hpp"

using namespace feff::par;

TEST(ParTest, InitSequential) {
    par_begin();
    auto& s = state();
    EXPECT_EQ(s.numprocs, 1);
    EXPECT_EQ(s.my_rank, 0);
    EXPECT_TRUE(s.master);
    EXPECT_FALSE(s.worker);
    EXPECT_FALSE(s.parallel_run);
    EXPECT_EQ(s.par_type, 1);
    par_end();
}

TEST(ParTest, DecompSingle) {
    // Single processor gets everything
    int start, end;
    mpe_decomp1d(100, 1, 0, start, end);
    EXPECT_EQ(start, 1);
    EXPECT_EQ(end, 100);
}

TEST(ParTest, DecompEven) {
    // 100 elements across 4 processors = 25 each
    int start, end;
    mpe_decomp1d(100, 4, 0, start, end);
    EXPECT_EQ(start, 1);
    EXPECT_EQ(end, 25);

    mpe_decomp1d(100, 4, 1, start, end);
    EXPECT_EQ(start, 26);
    EXPECT_EQ(end, 50);

    mpe_decomp1d(100, 4, 3, start, end);
    EXPECT_EQ(start, 76);
    EXPECT_EQ(end, 100);
}

TEST(ParTest, DecompUneven) {
    // 10 elements across 3 processors: 4, 3, 3
    int start, end;
    mpe_decomp1d(10, 3, 0, start, end);
    EXPECT_EQ(start, 1);
    EXPECT_EQ(end, 4);  // gets extra

    mpe_decomp1d(10, 3, 1, start, end);
    EXPECT_EQ(start, 5);
    EXPECT_EQ(end, 7);

    mpe_decomp1d(10, 3, 2, start, end);
    EXPECT_EQ(start, 8);
    EXPECT_EQ(end, 10);
}

TEST(ParTest, Seconds) {
    double t1 = seconds();
    double t2 = seconds();
    EXPECT_GE(t2, t1);
}
