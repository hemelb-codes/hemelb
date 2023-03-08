// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/GmyReadResult.h"
#include "geometry/LookupTree.h"

namespace hemelb::geometry {
    GmyReadResult::GmyReadResult(const util::Vector3D<site_t>& dimensionsInBlocks, site_t blockSize) :
            dimensionsInBlocks(dimensionsInBlocks), blockSize(blockSize),
            blockCount(dimensionsInBlocks.x() * dimensionsInBlocks.y() * dimensionsInBlocks.z()),
            sitesPerBlock(util::NumericalFunctions::IntegerPower(blockSize, 3)),
            Blocks(blockCount)
    {
    }

    GmyReadResult::~GmyReadResult() {

    }
    site_t GmyReadResult::FindSiteIndexInBlock(site_t fluidSiteBlock, site_t fluidSitesToPass) const
    {
        site_t siteIndex = 0;
        while (true)
        {
            // We keep going through the sites on the block until we've passed as many fluid
            // sites as we need to.
            if (Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor != SITE_OR_BLOCK_SOLID)
            {
                fluidSitesToPass--;
            }
            if (fluidSitesToPass < 0)
            {
                break;
            }
            siteIndex++;
        }
        return siteIndex;
    }

    site_t GmyReadResult::FindFluidSiteIndexInBlock(site_t fluidSiteBlock, site_t neighbourSiteId) const
    {
        site_t SiteId = 0;
        // Calculate the site's id over the whole geometry,
        for (site_t neighSite = 0; neighSite < GetSitesPerBlock(); neighSite++)
        {
            if (neighSite == neighbourSiteId)
            {
                break;
            }
            else if (Blocks[fluidSiteBlock].Sites[neighSite].targetProcessor != SITE_OR_BLOCK_SOLID)
            {
                SiteId++;
            }
        }

        return SiteId;
    }
}