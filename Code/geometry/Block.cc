
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "constants.h"
#include "geometry/Block.h"

namespace hemelb
{
  namespace geometry
  {
    const site_t Block::SOLID_SITE_ID = 1U << 31;

    Block::Block()
    {
    }

    Block::Block(site_t sitesPerBlock) :
        processorRankForEachBlockSite(sitesPerBlock, SITE_OR_BLOCK_SOLID), localContiguousIndex(sitesPerBlock, SOLID_SITE_ID)
    {
    }

    Block::~Block()
    {
    }

    bool Block::IsEmpty() const
    {
      return localContiguousIndex.empty();
    }

    proc_t Block::GetProcessorRankForSite(site_t localSiteIndex) const
    {
      return processorRankForEachBlockSite[localSiteIndex];
    }

    site_t Block::GetLocalContiguousIndexForSite(site_t localSiteIndex) const
    {
      return localContiguousIndex[localSiteIndex];
    }

    bool Block::SiteIsSolid(site_t localSiteIndex) const
    {
      return localContiguousIndex[localSiteIndex] == SOLID_SITE_ID;
    }

    void Block::SetProcessorRankForSite(site_t localSiteIndex, proc_t rank)
    {
      processorRankForEachBlockSite[localSiteIndex] = rank;
    }

    void Block::SetLocalContiguousIndexForSite(site_t localSiteIndex, site_t contiguousIndex)
    {
      localContiguousIndex[localSiteIndex] = contiguousIndex;
    }

  }
}
