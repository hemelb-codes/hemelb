#include "benchmark.h"
#ifndef NOMPI
#include <mpi.h>
#endif

// Returns true if elapsed time ? target amount of time.
bool BenchmarkTimer::IsBenchSectionFinished (double minutes, double elapsed_time_in_seconds)
{
  int is_bench_section_finished = 0;
  
  if (elapsed_time_in_seconds > minutes * 60.)
    {
      is_bench_section_finished = 1;
    }

#ifndef NOMPI
  MPI_Bcast (&is_bench_section_finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  return is_bench_section_finished > 0;
}

