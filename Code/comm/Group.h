
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_GROUP_H
#define HEMELB_COMM_GROUP_H

#include "units.h"
#include <vector>
#include <memory>

namespace hemelb
{
  namespace comm
  {
    class Group
    {
      public:
      typedef std::shared_ptr<Group> Ptr;
        /**
         * Returns the local rank within the group
         * @return
         */
        virtual int Rank() const = 0;

        /**
         * Returns the size of the group
         * @return
         */
        virtual int Size() const = 0;

        /**
         * Exclude the provided ranks
         * @param ranksToExclude
         * @return
         */
        virtual Group::Ptr Exclude(const std::vector<proc_t>& ranksToExclude) const = 0;
        /**
         * Include the provided ranks
         * @param ranksToExclude
         * @return
         */
        virtual Group::Ptr Include(const std::vector<proc_t>& ranksToInclude) const = 0;
    };
  }
}

#endif /* HEMELB_COMM_GROUP_H */
