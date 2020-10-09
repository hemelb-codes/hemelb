
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_MPIGROUP_H
#define HEMELB_COMM_MPIGROUP_H

#include "comm/Group.h"
#include <mpi.h>
#include <memory>

namespace hemelb
{
  namespace comm
  {
    class MpiGroup : public Group
    {
      public:
        /**
         * Default c'tor - initialises equivalent to MPI_GROUP_NULL
         */
        MpiGroup();
        // Direct construction.
        /**
         * Construct from an MPI_Group
         * @param grp
         * @param own - Whether this instance is responsible for deleting the
         * group when all copies are destroyed.
         */
        MpiGroup(MPI_Group grp, bool own);

	virtual ~MpiGroup() = default;
        /**
         * Returns the local rank within the group
         * @return
         */
        virtual int Rank() const;

        /**
         * Returns the size of the group
         * @return
         */
        virtual int Size() const;

        /**
         * Exclude the provided ranks
         * @param ranksToExclude
         * @return
         */
        virtual Group::Ptr Exclude(const std::vector<proc_t>& ranksToExclude) const;
        /**
         * Include the provided ranks
         * @param ranksToExclude
         * @return
         */
        virtual Group::Ptr Include(const std::vector<proc_t>& ranksToInclude) const;

        /**
         * Implicit cast to the underlying MPI_group
         * @return
         */
        operator MPI_Group() const
        {
          return *groupPtr;
        }

      private:
        std::shared_ptr<MPI_Group> groupPtr;
    };
  }
}

#endif /* HEMELB_COMM_MPIGROUP_H */
