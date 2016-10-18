
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/MpiError.h"
#include <sstream>

namespace hemelb
{
  namespace net
  {
    MpiError::MpiError(const char* mpiFunc_, const int errorCode_, const char* fileName_,
                       const int lineNo_) :
      mpiFunc(mpiFunc_), errorCode(errorCode_), fileName(fileName_), lineNo(lineNo_)
    {
      // Construct the first bit of error message
      *this << fileName << ":" << lineNo << ": " << mpiFunc << ": ";

      // Query the MPI implementation for its reason for failure
      char buffer[MPI_MAX_ERROR_STRING];
      int len;
      MPI_Error_string(errorCode, buffer, &len);

      *this << buffer;
    }

  }
}
