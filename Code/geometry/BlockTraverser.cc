// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
