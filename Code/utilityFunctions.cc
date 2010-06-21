#ifndef NOMPI
#ifdef XT3
#include <mpi.h>
#else
#include "mpi.h"
#endif
#endif

#include "utilityFunctions.h"

// Return the smaller of the two ints
int UtilityFunctions::min (int a, int b)
{
  if (a < b)
    {
      return a;
    }
  else
    {
      return b;
    }
}

//Return the larger of the two ints
int UtilityFunctions::max (int a, int b)
{
  if (a > b)
    {
      return a;
    }
  else
    {
      return b;
    }
}

// If number < lowerBound, returns lowerBound. If number > upperBound, returns upperBound
// If number between bounds, returns number.
// Consider the behaviour undefined if lowerBound > upperBound, so don't try it!
int UtilityFunctions::enforceBounds(int number, int lowerBound, int upperBound)
{
  return max(lowerBound, min(number, upperBound));
}

// Returns the number of seconds to 6dp elapsed since the Epoch
double UtilityFunctions::myClock ()
{
#ifdef NOMPI
  struct timeval time_data;
  
  gettimeofday (&time_data, NULL);
  
  return (double)time_data.tv_sec + (double)time_data.tv_usec / 1.e6;
#else
  return MPI_Wtime();
#endif
}