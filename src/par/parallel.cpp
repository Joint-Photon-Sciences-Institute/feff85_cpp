// Parallel/sequential abstraction layer — sequential implementation.
// Converted from: src/PAR/sequential.f

#include "parallel.hpp"
#include <iostream>
#include <cstdlib>

namespace feff::par {

static ParallelState g_state;

ParallelState& state() {
    return g_state;
}

void par_begin() {
    g_state.numprocs = 1;
    g_state.my_rank = 0;
    g_state.this_process = 0;
    g_state.master = true;
    g_state.worker = false;
    g_state.parallel_run = false;
    g_state.par_type = 1;
}

void par_end() {
    // No-op in sequential mode
}

[[noreturn]] void par_stop(const std::string& msg) {
    std::cerr << " Fatal error: " << msg << std::endl;
    std::exit(1);
}

void par_barrier() {
    // No-op in sequential mode
}

void mpe_decomp1d(int n, int num_procs, int myid, int& start, int& end) {
    // Balanced decomposition of [1..n] across num_procs processors.
    // Identical to Fortran MPE_DECOMP1D.
    int nlocal = n / num_procs;
    int deficit = n % num_procs;

    // Processors with rank < deficit get one extra element
    if (myid < deficit) {
        start = myid * (nlocal + 1) + 1;
        end = start + nlocal;
    } else {
        start = myid * nlocal + deficit + 1;
        end = start + nlocal - 1;
    }

    // Clamp to valid range
    if (end > n) end = n;
    if (start > end) start = end + 1;  // empty range
}

double seconds() {
    using clock = std::chrono::high_resolution_clock;
    static auto t0 = clock::now();
    auto now = clock::now();
    std::chrono::duration<double> elapsed = now - t0;
    return elapsed.count();
}

} // namespace feff::par
