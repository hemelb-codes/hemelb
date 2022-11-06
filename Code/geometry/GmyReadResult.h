// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_GMYREADRESULT_H
#define HEMELB_GEOMETRY_GMYREADRESULT_H

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
    class GmyReadResult
    {
      public:
        /**
         * Default constructor initialises internal variables
         */
        GmyReadResult(const util::Vector3D<site_t>& dimensionsInBlocks, site_t blockSize);

        /***
         * Returns the total count of blocks in the bounding box of the geometry.
         * @return count of blocks in the geometry.
         */
        inline site_t GetBlockCount() const
        {
          return blockCount;
        }

        /**
         * Returns the number of sites in each cubic block.
         * @return Number of sites in each block of the geometry.
         */
        inline site_t GetSitesPerBlock() const
        {
          return sitesPerBlock;
        }

        /***
         * Get the i.d. of a block, i.e. the one-d coordinate, from the three-d coordinate.
         * @todo Use this to replace domain_type::GetBlockIdFromBlockCoords
         */
        inline site_t GetBlockIdFromBlockCoordinates(site_t blockI, site_t blockJ, site_t blockK) const
        {
          return (blockI * dimensionsInBlocks.y + blockJ) * dimensionsInBlocks.z + blockK;
        }

        /***
         * Get the coordinates of a block from a block id.
         */
        inline util::Vector3D<site_t> GetBlockCoordinatesFromBlockId(site_t blockId) const
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
        inline site_t GetSiteIdFromSiteCoordinates(site_t siteI, site_t siteJ, site_t siteK) const
        {
          return (siteI * blockSize + siteJ) * blockSize + siteK;
        }

        /**
         * True if the given block coordinates are within the geometry bounding-box.
         */
        inline bool AreBlockCoordinatesValid(const util::Vector3D<site_t>& blockCoords) const
        {
          return blockCoords.x >= 0 && blockCoords.y >= 0 && blockCoords.z >= 0
              && blockCoords.x < dimensionsInBlocks.x && blockCoords.y < dimensionsInBlocks.y
              && blockCoords.z < dimensionsInBlocks.z;
        }

        /**
         * True if the given site coordinates are within a bounding-box for a block.
         */
        inline bool AreLocalSiteCoordinatesValid(const util::Vector3D<site_t>& siteCoords) const
        {
          return siteCoords.x >= 0 && siteCoords.y >= 0 && siteCoords.z >= 0
              && siteCoords.x < blockSize && siteCoords.y < blockSize && siteCoords.z < blockSize;
        }

        /**
         * Get the dimensions of the bounding box in terms of blocks.
         * @return Dimensions of the bounding box in blocks.
         */
        inline util::Vector3D<site_t> const& GetBlockDimensions() const
        {
          return dimensionsInBlocks;
        }

        /**
         * Get the number of sites along one edge of a block.
         * @return Number of sites that make up one block-length.
         */
        inline site_t GetBlockSize() const
        {
          return blockSize;
        }

        /* Find a site index, taking into account ALL lattice sites. */
        site_t FindSiteIndexInBlock(site_t fluidSiteBlock, site_t fluidSitesToPass) const;

        /* Find a site index, taking into account ONLY fluid sites. */
        site_t FindFluidSiteIndexInBlock(site_t fluidSiteBlock, site_t neighbourSiteId) const;

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
#endif // HEMELB_GEOMETRY_GMYREADRESULT_H
