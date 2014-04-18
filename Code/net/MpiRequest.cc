//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "net/MpiRequest.h"
#include "net/MpiStatus.h"

namespace hemelb
{
  namespace net
  {
    MpiRequest::MpiRequest() :
        reqPtr()
    {
    }
    MpiRequest::MpiRequest(MPI_Request req) :
        reqPtr()
    {
      reqPtr.reset(new MPI_Request(req));
    }

    void MpiRequest::Wait()
    {
      HEMELB_MPI_CALL(MPI_Wait, (reqPtr.get(), MPI_STATUS_IGNORE));
    }
    void MpiRequest::Wait(MpiStatus& stat)
    {
      HEMELB_MPI_CALL(MPI_Wait, (reqPtr.get(), stat.statPtr.get()));
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

    void MpiRequest::WaitAll(ReqVec& reqs, StatVec& stats)
    {
      size_t n = reqs.size();
      MPI_Request* rawreqs = new MPI_Request[n];
      MPI_Status* rawstats = new MPI_Status[n];
      try
      {
        for (size_t i = 0; i < n; ++i)
        {
          rawreqs[i] = reqs[i];
        }

        HEMELB_MPI_CALL(
            MPI_Waitall,
            (n, rawreqs, rawstats)
        );

        for (size_t i = 0; i < n; ++i)
        {
          stats[i] = MpiStatus(rawstats[i]);
        }
      }
      catch (const std::exception& e)
      {
        delete[] rawreqs;
        delete[] rawstats;
        throw;
      }

      delete[] rawreqs;
      delete[] rawstats;
    }
    /*
    void MpiRequest::WaitAll(ReqVec& reqs, StatVec& stats)
    {
      stats.resize(reqs.size());
      ReqVec::iterator reqIt = reqs.begin();
      StatVec::iterator statIt = stats.begin();
      for (; reqIt != reqs.end(); ++reqIt, ++statIt)
      {
        *statIt = reqIt ->Wait();
      }
    }
    */

    bool MpiRequest::Test()
    {
      int flag;
      HEMELB_MPI_CALL(MPI_Test, (reqPtr.get(), &flag, MPI_STATUS_IGNORE));
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
