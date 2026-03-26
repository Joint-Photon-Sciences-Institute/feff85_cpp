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

   - ``numprocs`` --- number of processes (1 for sequential)
   - ``my_rank`` --- MPI rank (0 for sequential)
   - ``this_process`` --- process index (0 for sequential)
   - ``master`` --- true if this is the master process
   - ``worker`` --- true if this is a worker process
   - ``parallel_run`` --- true if running in parallel mode
   - ``par_type`` --- output mode (0=not initialized, 1=normal, 2=suppress output, 3=no log file)

Functions
---------

``ParallelState& state()``
   Global accessor for the parallel state singleton.

``void par_begin()``
   Initialize the parallel execution environment (``MPI_Init`` in parallel mode).

``void par_end()``
   Normal shutdown (``MPI_Finalize`` in parallel mode).

``void par_stop(const std::string& msg)``
   Abnormal termination with error message. Marked ``[[noreturn]]``.

``void par_barrier()``
   Synchronization barrier (no-op in sequential mode).

``void mpe_decomp1d(int n, int num_procs, int myid, int& start, int& end)``
   Decompose a 1D range ``[1, n]`` across processes. Returns the ``start``
   and ``end`` indices (1-based, inclusive) for process ``myid``.

``double seconds()``
   Wall clock time in seconds.
