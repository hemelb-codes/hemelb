#include "benchmark.h"
#include <mpi.h>

int IsBenchSectionFinished (double minutes, double elapsed_time)
{
  int is_bench_section_finished = 0;
  
  if (elapsed_time > minutes * 60.)
    {
      is_bench_section_finished = 1;
    }
#ifndef NOMPI
  MPI_Bcast (&is_bench_section_finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  
  if (is_bench_section_finished)
    {
      return 1;
    }
  return 0;
}

