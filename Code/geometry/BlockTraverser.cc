// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/BlockTraverser.h"
#include "geometry/Domain.h"
#include "util/Vector3D.h"

namespace hemelb::geometry
{
    BlockTraverser::BlockTraverser(const geometry::Domain& iLatDat) :
        VolumeTraverser(), mLatticeData(iLatDat)
    {
    }

    site_t BlockTraverser::CurrentBlockNumber() const
    {
      return GetCurrentIndex();
    }

    util::Vector3D<site_t> BlockTraverser::GetSiteCoordinatesOfLowestSiteInCurrentBlock()
    {
      return GetCurrentLocation().as<site_t>() * mLatticeData.GetBlockSize();
    }

    const Block& BlockTraverser::GetCurrentBlockData()
    {
      return mLatticeData.GetBlock(GetCurrentLocation());
    }

    const Block& BlockTraverser::GetBlockDataForLocation(const Vec16& iLocation)
    {
      return mLatticeData.GetBlock(iLocation);
    }

    site_t BlockTraverser::GetBlockSize()
    {
      return mLatticeData.GetBlockSize();
    }

    SiteTraverser BlockTraverser::GetSiteTraverser()
    {
      return SiteTraverser(mLatticeData);
    }

    bool BlockTraverser::IsValidLocation(Vec16 const& iBlock)
    {
      return mLatticeData.IsValidBlock(iBlock);
    }

    bool BlockTraverser::GoToNextBlock()
    {
      return TraverseOne();
    }

    U16 BlockTraverser::GetXCount() const
    {
      return mLatticeData.GetBlockDimensions().x();
    }

    U16 BlockTraverser::GetYCount() const
    {
      return mLatticeData.GetBlockDimensions().y();
    }

    U16 BlockTraverser::GetZCount() const
    {
      return mLatticeData.GetBlockDimensions().z();
    }
}
