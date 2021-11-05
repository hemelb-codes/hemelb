// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_GEOMETRY_H
#define HEMELB_GEOMETRY_GEOMETRY_H

#include <vector>
#include "units.h"
#include "constants.h"
#include "geometry/GeometryBlock.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
    /***
     * Model of the information in a geometry file
     */
    class Geometry
    {
      public:
        /**
         * Default constructor initialises internal variables
         */
        Geometry(const util::Vector3D<site_t>& dimensionsInBlocks, site_t blockSize) :
            dimensionsInBlocks(dimensionsInBlocks), blockSize(blockSize),
                blockCount(dimensionsInBlocks.x * dimensionsInBlocks.y * dimensionsInBlocks.z),
                sitesPerBlock(util::NumericalFunctions::IntegerPower(blockSize, 3)),
                Blocks(blockCount)
        {

        }

        /***
         * Returns the total count of blocks in the bounding box of the geometry.
         * @return count of blocks in the geometry.
         */
        site_t GetBlockCount() const
        {
          return blockCount;
        }

        /**
         * Returns the number of sites in each cubic block.
         * @return Number of sites in each block of the geometry.
         */
        site_t GetSitesPerBlock() const
        {
          return sitesPerBlock;
        }

        /***
         * Get the i.d. of a block, i.e. the one-d coordinate, from the three-d coordinate.
         * @todo Use this to replace LatticeData::GetBlockIdFromBlockCoords
         */
        site_t GetBlockIdFromBlockCoordinates(site_t blockI, site_t blockJ, site_t blockK) const
        {
          return (blockI * dimensionsInBlocks.y + blockJ) * dimensionsInBlocks.z + blockK;
        }

        /***
         * Get the coordinates of a block from a block id.
         */
        util::Vector3D<site_t> GetBlockCoordinatesFromBlockId(site_t blockId) const
        {
          site_t blockZ = blockId % dimensionsInBlocks.z;
          site_t remainder = blockId / dimensionsInBlocks.z;
          site_t blockY = remainder % dimensionsInBlocks.y;
          site_t blockX = remainder / dimensionsInBlocks.y;
          return util::Vector3D<site_t>(blockX, blockY, blockZ);
        }

        /***
         * Get the i.d. of a site, i.e. the one-d coordinate, from the three-d coordinate.
         */
        site_t GetSiteIdFromSiteCoordinates(site_t siteI, site_t siteJ, site_t siteK) const
        {
          return (siteI * blockSize + siteJ) * blockSize + siteK;
        }

        /**
         * True if the given block coordinates are within the geometry bounding-box.
         */
        bool AreBlockCoordinatesValid(const util::Vector3D<site_t>& blockCoords) const
        {
          return blockCoords.x >= 0 && blockCoords.y >= 0 && blockCoords.z >= 0
              && blockCoords.x < dimensionsInBlocks.x && blockCoords.y < dimensionsInBlocks.y
              && blockCoords.z < dimensionsInBlocks.z;
        }

        /**
         * True if the given site coordinates are within a bounding-box for a block.
         */
        bool AreLocalSiteCoordinatesValid(const util::Vector3D<site_t>& siteCoords) const
        {
          return siteCoords.x >= 0 && siteCoords.y >= 0 && siteCoords.z >= 0
              && siteCoords.x < blockSize && siteCoords.y < blockSize && siteCoords.z < blockSize;
        }

        /**
         * Get the dimensions of the bounding box in terms of blocks.
         * @return Dimensions of the bounding box in blocks.
         */
        const util::Vector3D<site_t>& GetBlockDimensions() const
        {
          return dimensionsInBlocks;
        }

        /**
         * Get the number of sites along one edge of a block.
         * @return Number of sites that make up one block-length.
         */
        site_t GetBlockSize() const
        {
          return blockSize;
        }

        /* Find a site index, taking into account ALL lattice sites. */
        site_t FindSiteIndexInBlock(site_t fluidSiteBlock, site_t fluidSitesToPass) const
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

        /* Find a site index, taking into account ONLY fluid sites. */
        site_t FindFluidSiteIndexInBlock(site_t fluidSiteBlock, site_t neighbourSiteId) const
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

      private:
        const util::Vector3D<site_t> dimensionsInBlocks; //! The count of blocks in each direction
        const site_t blockSize; //! Size of a block, in sites.

        const site_t blockCount;
        const site_t sitesPerBlock;

      public:
        std::vector<BlockReadResult> Blocks; //! Array of Block models

    };

  }
}
#endif // HEMELB_GEOMETRY_GEOMETRYREADRESULT_H
