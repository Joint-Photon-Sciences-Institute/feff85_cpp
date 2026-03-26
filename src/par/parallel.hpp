#pragma once
// Parallel/sequential abstraction layer.
// Converted from: src/PAR/parallel.f, src/PAR/sequential.f
//
// Currently implements sequential-only mode. MPI support can be added
// by compiling with FEFF_USE_MPI and linking against MPI::MPI_CXX.

#include <string>
#include <chrono>
#include <stdexcept>

namespace feff::par {

/// Parallel execution state (replaces Fortran COMMON /parallel/)
struct ParallelState {
    int numprocs = 1;
    int my_rank = 0;
    int this_process = 0;
    bool master = true;
    bool worker = false;
    bool parallel_run = false;
    // par_type: 0=not initialized, 1=normal, 2=suppress output (seq loop), 3=no log file
    int par_type = 0;
};

/// Global parallel state accessor
ParallelState& state();

/// Initialize parallel environment (MPI_Init in parallel mode)
void par_begin();

/// Normal shutdown (MPI_Finalize in parallel mode)
void par_end();

/// Abnormal termination with error message
[[noreturn]] void par_stop(const std::string& msg);

/// Synchronization barrier (MPI_Barrier in parallel mode)
void par_barrier();

/// Decompose 1D array [1..n] across num_procs processors.
/// Returns the start and end indices for processor myid (1-based, inclusive).
void mpe_decomp1d(int n, int num_procs, int myid, int& start, int& end);

/// Wall clock time in seconds
double seconds();

} // namespace feff::par
