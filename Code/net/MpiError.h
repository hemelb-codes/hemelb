//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_MPIERROR_H
#define HEMELB_NET_MPIERROR_H

#include <mpi.h>
#include <string>

namespace hemelb
{
  namespace net
  {
    class MpiError : public std::exception
    {
      public:
        MpiError(const char* mpiFunc_, int errorCode_, const char* fileName_, const int lineNo);
        virtual ~MpiError() throw ();

        virtual const char * what() const throw ();
      private:
        const char* mpiFunc;
        const int errorCode;
        const char* fileName;
        const int lineNo;
        std::string message;
    };
  }
}

#define HEMELB_MPI_CALL( mpiFunc, args ) \
{ \
  int _check_result = mpiFunc args; \
  if (_check_result != MPI_SUCCESS) \
    throw ::hemelb::net::MpiError(#mpiFunc, _check_result, __FILE__, __LINE__); \
}

#endif // HEMELB_NET_MPIERROR_H
