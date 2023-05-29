// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/Block.h"
#include "hassert.h"
#include "constants.h"

namespace hemelb::geometry
{

    Block::Block(site_t sitesPerBlock) :
        processorRankForEachBlockSite(sitesPerBlock, SITE_OR_BLOCK_SOLID),
        localContiguousIndex(sitesPerBlock, SOLID_SITE_ID)
    {
    }

    bool Block::IsEmpty() const
    {
      return localContiguousIndex.empty();
    }

    proc_t Block::GetProcessorRankForSite(site_t localSiteIndex) const
    {
      HASSERT(localSiteIndex < processorRankForEachBlockSite.size());
      return processorRankForEachBlockSite[localSiteIndex];
    }

    site_t Block::GetLocalContiguousIndexForSite(site_t localSiteIndex) const
    {
      HASSERT(localSiteIndex < localContiguousIndex.size());
      return localContiguousIndex[localSiteIndex];
    }

    bool Block::SiteIsSolid(site_t localSiteIndex) const
    {
      HASSERT(localSiteIndex < localContiguousIndex.size());
      return localContiguousIndex[localSiteIndex] == SOLID_SITE_ID;
    }

    void Block::SetProcessorRankForSite(site_t localSiteIndex, proc_t rank)
    {
      HASSERT(localSiteIndex < processorRankForEachBlockSite.size());
      processorRankForEachBlockSite[localSiteIndex] = rank;
    }

    void Block::SetLocalContiguousIndexForSite(site_t localSiteIndex, site_t contiguousIndex)
    {
      HASSERT(localSiteIndex < localContiguousIndex.size());
      localContiguousIndex[localSiteIndex] = contiguousIndex;
    }

}
