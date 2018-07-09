
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/BlockTraverser.h"
#include "geometry/LatticeData.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
    BlockTraverser::BlockTraverser(const geometry::LatticeData& iLatDat) :
        VolumeTraverser(), mLatticeData(iLatDat)
    {
    }

    BlockTraverser::~BlockTraverser()
    {
    }

    site_t BlockTraverser::CurrentBlockNumber() const
    {
      return GetCurrentIndex();
    }

    util::Vector3D<site_t> BlockTraverser::GetSiteCoordinatesOfLowestSiteInCurrentBlock()
    {
      return GetCurrentLocation() * mLatticeData.GetBlockSize();
    }

    const Block& BlockTraverser::GetCurrentBlockData()
    {
      return mLatticeData.GetBlock(GetCurrentIndex());
    }

    const Block& BlockTraverser::GetBlockDataForLocation(const util::Vector3D<site_t>& iLocation)
    {
      return mLatticeData.GetBlock(GetIndexFromLocation(iLocation));
    }

    site_t BlockTraverser::GetBlockSize()
    {
      return mLatticeData.GetBlockSize();
    }

    SiteTraverser BlockTraverser::GetSiteTraverser()
    {
      return SiteTraverser(mLatticeData);
    }

    bool BlockTraverser::IsValidLocation(util::Vector3D<site_t> iBlock)
    {
      return mLatticeData.IsValidBlock(iBlock.x, iBlock.y, iBlock.z);
    }

    bool BlockTraverser::GoToNextBlock()
    {
      return TraverseOne();
    }

    site_t BlockTraverser::GetXCount() const
    {
      return mLatticeData.GetBlockDimensions().x;
    }

    site_t BlockTraverser::GetYCount() const
    {
      return mLatticeData.GetBlockDimensions().y;
    }

    site_t BlockTraverser::GetZCount() const
    {
      return mLatticeData.GetBlockDimensions().z;
    }
  }
}
