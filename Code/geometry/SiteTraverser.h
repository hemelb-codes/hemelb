
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_SITETRAVERSER_H
#define HEMELB_GEOMETRY_SITETRAVERSER_H

#include "geometry/VolumeTraverser.h"

namespace hemelb
{
  namespace geometry
  {
    class LatticeData;

    /**
     * SiteTraverser is used to traverse the sites in a speficied block
     * within the lattice data
     */
    class SiteTraverser : public VolumeTraverser
    {
      public:
        SiteTraverser(const LatticeData & iLatticeDat);

        virtual ~SiteTraverser();

        site_t GetXCount() const;

        site_t GetYCount() const;

        site_t GetZCount() const;

      private:
        //Returns the block size in number of sites
        site_t GetBlockSize() const;

        const LatticeData & mLatticeData;
    };

  }
}

#endif // HEMELB_GEOMETRY_SITETRAVERSER_H
