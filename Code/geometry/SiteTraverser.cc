
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/LatticeData.h"
#include "geometry/SiteTraverser.h"

namespace hemelb
{
  namespace geometry
  {
    SiteTraverser::SiteTraverser(const geometry::LatticeData& iLatticeData) :
      mLatticeData(iLatticeData)

    {
    }

    SiteTraverser::~SiteTraverser()
    {
    }

    site_t SiteTraverser::GetXCount() const
    {
      return GetBlockSize();
    }

    site_t SiteTraverser::GetYCount() const
    {
      return GetBlockSize();
    }

    site_t SiteTraverser::GetZCount() const
    {
      return GetBlockSize();
    }

    site_t SiteTraverser::GetBlockSize() const
    {
      return mLatticeData.GetBlockSize();
    }
  }
}
