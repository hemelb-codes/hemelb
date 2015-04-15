//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "net/MpiCommunicator.h"
#include "net/MpiGroup.h"

namespace hemelb
{
  namespace net
  {
    namespace
    {
      void Deleter(MPI_Comm* comm)
      {
        int finalized;
        HEMELB_MPI_CALL(MPI_Finalized, (&finalized));
        if (!finalized)
          HEMELB_MPI_CALL(MPI_Comm_free, (comm));
        delete comm;
      }
    }

    MpiCommunicator MpiCommunicator::World()
    {
      return MpiCommunicator(MPI_COMM_WORLD, false);
    }

    MpiCommunicator::MpiCommunicator() :
        commPtr()
    {
    }

    MpiCommunicator::MpiCommunicator(MPI_Comm communicator, bool owner) :
        commPtr()
    {
      if (communicator == MPI_COMM_NULL)
        return;

      if (owner)
      {
        commPtr.reset(new MPI_Comm(communicator), Deleter);
      }
      else
      {
        commPtr.reset(new MPI_Comm(communicator));
      }
    }

    MpiCommunicator::~MpiCommunicator()
    {
    }

    bool operator==(const MpiCommunicator& comm1, const MpiCommunicator& comm2)
    {
      if (comm1)
      {
        if (comm2)
        {
          int result;
          HEMELB_MPI_CALL(MPI_Comm_compare, (comm1, comm2, &result));
          return result == MPI_IDENT;
        }
        return false;
      }
      return (!comm2);
    }

    bool operator!=(const MpiCommunicator& comm1, const MpiCommunicator& comm2)
    {
      return ! (comm1 == comm2);
    }

    int MpiCommunicator::Rank() const
    {
      int rank;
      HEMELB_MPI_CALL(MPI_Comm_rank, (*commPtr, &rank));
      return rank;
    }

    int MpiCommunicator::Size() const
    {
      int size;
      HEMELB_MPI_CALL(MPI_Comm_size, (*commPtr, &size));
      return size;
    }

    MpiGroup MpiCommunicator::Group() const
    {
      MPI_Group grp;
      HEMELB_MPI_CALL(MPI_Comm_group, (*commPtr, &grp));
      return MpiGroup(grp, true);
    }

    MpiCommunicator MpiCommunicator::Create(const MpiGroup& grp) const
    {
      MPI_Comm newComm;
      HEMELB_MPI_CALL(MPI_Comm_create, (*commPtr, grp, &newComm));
      return MpiCommunicator(newComm, true);
    }

    void MpiCommunicator::Abort(int errCode) const
    {
      HEMELB_MPI_CALL(MPI_Abort, (*commPtr, errCode));
    }

    MpiCommunicator MpiCommunicator::Duplicate() const
    {
      MPI_Comm newComm;
      HEMELB_MPI_CALL(MPI_Comm_dup, (*commPtr, &newComm));
      return MpiCommunicator(newComm, true);
    }
  }
}
