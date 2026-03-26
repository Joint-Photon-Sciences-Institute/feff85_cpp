PAR Module
==========

**Namespace:** ``feff::par``

**Header:** ``src/par/parallel.hpp``

The PAR module provides an abstraction layer for parallel execution. In the
current implementation, all calculations run sequentially, but the module
provides hooks for future MPI parallelization.

Data Structures
---------------

``ParallelState``
   Holds the parallel execution state:

   - ``par_type`` --- parallel mode (0 = sequential)
   - ``this_process`` --- process rank (0 for sequential)
   - ``numprocs`` --- number of processes (1 for sequential)

Functions
---------

``void par_begin(ParallelState& state, int argc, char* argv[])``
   Initialize the parallel execution environment.

``void par_end(ParallelState& state)``
   Shut down the parallel execution environment.

``void par_barrier(ParallelState& state)``
   Synchronization barrier (no-op in sequential mode).

``void mpe_decomp1d(int n, int numprocs, int myid, int& s, int& e)``
   Decompose a 1D range ``[1, n]`` across processes. Returns the start ``s``
   and end ``e`` indices for process ``myid``.
