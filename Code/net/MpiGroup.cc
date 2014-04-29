//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "net/MpiGroup.h"
#include "net/MpiConstness.h"

namespace hemelb
{
  namespace net
  {
    namespace
    {
      void Deleter(MPI_Group* comm)
      {
        int finalized;
        HEMELB_MPI_CALL(MPI_Finalized, (&finalized));
        if (!finalized)
          HEMELB_MPI_CALL(MPI_Group_free, (comm));
        delete comm;
      }
    }

    MpiGroup::MpiGroup(): groupPtr()
    {
    }

    int MpiGroup::Rank() const
    {
      int rank;
      HEMELB_MPI_CALL(MPI_Group_rank, (*groupPtr, &rank));
      return rank;
    }

    int MpiGroup::Size() const
    {
      int size;
      HEMELB_MPI_CALL(MPI_Group_size, (*groupPtr, &size));
      return size;
    }

    MpiGroup MpiGroup::Exclude(const std::vector<int>& ranksToExclude)
    {
      MPI_Group ans;
      HEMELB_MPI_CALL(MPI_Group_excl,
          (*groupPtr, ranksToExclude.size(), MpiConstCast(&ranksToExclude.front()), &ans))
      return MpiGroup(ans, true);
    }

    MpiGroup MpiGroup::Include(const std::vector<int>& ranksToInclude)
    {
      MPI_Group ans;
      HEMELB_MPI_CALL(MPI_Group_incl,
                      (*groupPtr, ranksToInclude.size(), MpiConstCast(&ranksToInclude.front()), &ans))
      return MpiGroup(ans, true);
    }

    MpiGroup::MpiGroup(MPI_Group grp, bool own)
    {
      if (grp != MPI_GROUP_EMPTY)
      {
        if (own)
        {
          groupPtr.reset(new MPI_Group(grp), Deleter);
        }
        else
        {
          groupPtr.reset(new MPI_Group(grp));
        }
      }
    }

  }
}
