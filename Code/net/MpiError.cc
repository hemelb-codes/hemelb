//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

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
      // Query the MPI implementation for its reason for failure
      char buffer[MPI_MAX_ERROR_STRING];
      int len;
      MPI_Error_string(errorCode_, buffer, &len);

      // Construct the complete error message
      {
        std::stringstream msg;
        msg << fileName << ":" << lineNo << ": ";
        msg << mpiFunc_ << ": ";
        message = msg.str();
      }
      message.append(buffer, len);
    }

    MpiError::~MpiError() throw ()
    {
    }

    const char* MpiError::what() const throw ()
    {
      return message.c_str();
    }
  }
}
