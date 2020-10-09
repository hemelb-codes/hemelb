
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_MPIERROR_H
#define HEMELB_COMM_MPIERROR_H

#include <mpi.h>
#include "Exception.h"

namespace hemelb
{
  namespace comm
  {
    /**
     * Indicate an error to do with MPI.
     *
     * Will be thrown by the HEMELB_MPI_CALL macro, as in:
     *   HEMELB_MPI_CALL(MPI_Send, (buffer, count, dtype, toRank, tag, communicator) );
     *
     */
    class MpiError : public ::hemelb::Exception
    {
      public:
        MpiError(const char* mpiFunc_, int errorCode_, const char* fileName_, const int lineNo);

        // Print this->what() to std::cerr without having to put
        // <iostream> in this header
        void LogError() const;
      private:
        const char* mpiFunc;
        const int errorCode;
        const char* fileName;
        const int lineNo;
    };
  }
}

// Wrap MPI calls with error checking. If not successful, will throw
// an exception which will print a nice error message.
#define HEMELB_MPI_CALL( mpiFunc, args ) \
{ \
  int _check_result = mpiFunc args; \
  if (_check_result != MPI_SUCCESS) \
    throw ::hemelb::comm::MpiError(#mpiFunc, _check_result, __FILE__, __LINE__); \
}

// For use in destructors and other nothrow functions
#define HEMELB_MPI_CALL_NOTHROW( mpiFunc, args ) {			\
  int _check_result = mpiFunc args;					\
  if (_check_result != MPI_SUCCESS) {					\
    ::hemelb::comm::MpiError(#mpiFunc, _check_result, __FILE__, __LINE__).LogError(); \
    MPI_Abort(MPI_COMM_WORLD, -1);					\
  }									\
}

#endif // HEMELB_COMM_MPIERROR_H
