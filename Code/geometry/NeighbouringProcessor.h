
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_NEIGHBOURINGPROCESSOR_H
#define HEMELB_GEOMETRY_NEIGHBOURINGPROCESSOR_H

#include "units.h"

namespace hemelb
{
  namespace geometry
  {
    /**
     * Encapsulates data about a processor which has at least 1 site which
     * neighbours at least one site local to this core.
     */
    struct NeighbouringProcessor
    {
      public:
        //! Rank of the neighbouring processor.
        proc_t Rank;

        //! The number of distributions shared between this neighbour and the current processor.
        //! Note that this is not equal to the number of interfacing sites * lattice connectivity
        //! because we only send the distributions that point towards the neighbour.
        site_t SharedDistributionCount;

        //! Index on this processor of the first distribution shared between this
        //! neighbour and the current processor.
        site_t FirstSharedDistribution;
    };
  }
}

#endif /* HEMELB_GEOMETRY_NEIGHBOURINGPROCESSOR_H */
