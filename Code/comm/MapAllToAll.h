//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_MAPALLTOALL_H
#define HEMELB_NET_MAPALLTOALL_H

#include <map>
#include "comm/Communicator.h"
#include "comm/Request.h"

namespace hemelb
{
  namespace comm
  {
    template <typename T>
    void MapAllToAll(Communicator::ConstPtr comm,
                     const std::map<int, T>& valsToSend,
                     std::map<int, T>& receivedVals,
                     const int tag = 1234)
    {
      receivedVals.clear();

      auto sendReqs = comm->MakeRequestList();
      int localRank = comm->Rank();
      int nSends = valsToSend.size();
      // Is localRank in the map of send destinations?
      if (valsToSend.find(localRank) != valsToSend.end())
        // If so, reduce the number of sends by 1.
        nSends -= 1;
      // Set up a container for all the Status objs
      sendReqs->resize(nSends);

      typename std::map<int, T>::const_iterator iter = valsToSend.begin(),
          end = valsToSend.end();
      int i = 0;
      for (; iter != end; ++iter) {
        int rank = iter->first;
        const T& val = iter->second;
        if (rank == localRank)
        {
          receivedVals[localRank] = val;
        }
        else
        {
          // Synchronous to ensure we know when this is matched
          sendReqs->set(i, std::move(*comm->Issend(val, rank, tag)));
          i++;
        }
      }

      bool allProcsHaveHadAllSendsMatched = false;
      bool mySendsMatched = false;
      Request::Ptr barrierReq;

      // Allocate outside of loop to ensure compiler isn't reinitialising this
      // all the time.
      MPI_Status status;
      while(!mySendsMatched || !allProcsHaveHadAllSendsMatched)
      {
        if (comm->Iprobe(MPI_ANY_SOURCE, tag, &status)) {
          // There is a message waiting
          comm->Recv(receivedVals[status.MPI_SOURCE], status.MPI_SOURCE, status.MPI_TAG);
        }

        if (!mySendsMatched)
        {
          if (sendReqs->TestAll())
          {
            mySendsMatched = true;
            barrierReq = comm->Ibarrier();
          }
        }
        else
        {
          allProcsHaveHadAllSendsMatched = barrierReq->Test();
        }
      }
    }
  }
}

#endif /* HEMELB_NET_MAPALLTOALL_H */
