// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/Domain.h"
#include "geometry/SiteTraverser.h"

namespace hemelb::geometry
{
    SiteTraverser::SiteTraverser(const geometry::Domain& iLatticeData) :
        mLatticeData(iLatticeData)

    {
    }

    U16 SiteTraverser::GetXCount() const
    {
      return GetBlockSize();
    }

    U16 SiteTraverser::GetYCount() const
    {
      return GetBlockSize();
    }

    U16 SiteTraverser::GetZCount() const
    {
      return GetBlockSize();
    }

    U16 SiteTraverser::GetBlockSize() const
    {
      return mLatticeData.GetBlockSize();
    }
  }
