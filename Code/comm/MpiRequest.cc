//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "comm/MpiRequest.h"

namespace hemelb
{
  namespace comm
  {
    MpiRequest::MpiRequest() :
        req(MPI_REQUEST_NULL)
    {
    }
    MpiRequest::MpiRequest(MPI_Request req_) :
        req(req_)
    {
    }

    MpiRequest::operator bool() const
    {
      return req != MPI_REQUEST_NULL;
    }

    void MpiRequest::Wait()
    {
      HEMELB_MPI_CALL(MPI_Wait, (&req, MPI_STATUS_IGNORE));
    }

    void MpiRequest::WaitAll(ReqVec& reqs)
    {
      size_t n = reqs.size();
      MPI_Request* rawreqs = new MPI_Request[n];
      try
      {
        for (size_t i = 0; i < n; ++i)
        {
          rawreqs[i] = reqs[i];
        }

        HEMELB_MPI_CALL(
            MPI_Waitall,
            (n, rawreqs, MPI_STATUSES_IGNORE)
        );
      }
      catch (const std::exception& e)
      {
        delete[] rawreqs;
        throw;
      }

      delete[] rawreqs;
    }

    bool MpiRequest::Test()
    {
      int flag;
      HEMELB_MPI_CALL(MPI_Test, (&req, &flag, MPI_STATUS_IGNORE));
      return flag;
    }

    bool MpiRequest::TestAll(ReqVec& reqs)
    {
      int flag;
      size_t n = reqs.size();
      if (n == 0)
        return true;

      MPI_Request* rawreqs = new MPI_Request[n];
      try
      {
        for (size_t i = 0; i < n; ++i)
        {
          rawreqs[i] = reqs[i];
        }

        HEMELB_MPI_CALL(
            MPI_Testall,
            (n, rawreqs, &flag, MPI_STATUSES_IGNORE)
        );
      }
      catch (const std::exception& e)
      {
        delete[] rawreqs;
        throw;
      }

      delete[] rawreqs;
      return flag;
    }

  }
}
