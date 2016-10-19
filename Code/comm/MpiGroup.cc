
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "comm/MpiGroup.h"
#include "comm/MpiError.h"

namespace hemelb
{
  namespace comm
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

    Group::Ptr MpiGroup::Exclude(const std::vector<int>& ranksToExclude) const
    {
      MPI_Group ans;
      HEMELB_MPI_CALL(MPI_Group_excl,
		      (*groupPtr, ranksToExclude.size(), ranksToExclude.data(), &ans));
      return std::make_shared<MpiGroup>(ans, true);
    }

    Group::Ptr MpiGroup::Include(const std::vector<int>& ranksToInclude) const
    {
      MPI_Group ans;
      HEMELB_MPI_CALL(MPI_Group_incl,
                      (*groupPtr, ranksToInclude.size(), ranksToInclude.data(), &ans));
      return std::make_shared<MpiGroup>(ans, true);
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
