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
        util::Vector3D<site_t> blocks; //! The count of blocks in each direction
        site_t blockSize; //! Size of a block, in sites.
        double voxelSize; //! Size of a block, in real-world units.
        util::Vector3D<double> origin;

        std::vector<BlockReadResult> Blocks; //! Array of Block models

        /***
         * Total count of blocks in the site.
         * @return count of blocks in the site.
         */
        site_t GetBlockCount() const
        {
          return blocks.x * blocks.y * blocks.z;
        }

        site_t GetSitesPerBlock() const
        {
          return util::NumericalFunctions::IntegerPower(blockSize, 3);
        }

        /***
         * Get the i.d. of a block, i.e. the one-d coordinate, from the three-d coordinate.
         * @todo Use this to replace LatticeData::GetBlockIdFromBlockCoords
         */
        site_t GetBlockIdFromBlockCoordinates(site_t blockI, site_t blockJ, site_t blockK) const
        {
          return (blockI * blocks.y + blockJ) * blocks.z + blockK;
        }

        /***
         * Get the coordinates of a block from a block id.
         */
        util::Vector3D<site_t> GetBlockCoordinatesFromBlockId(site_t blockId) const
        {
          site_t blockZ = blockId % blocks.z;
          site_t remainder = blockId / blocks.z;
          site_t blockY = remainder % blocks.y;
          site_t blockX = remainder / blocks.y;
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
          return blockCoords.x >= 0 && blockCoords.y >= 0 && blockCoords.z >= 0 && blockCoords.x < blocks.x
              && blockCoords.y < blocks.y && blockCoords.z < blocks.z;
        }
    };

  }
}
#endif // HEMELB_GEOMETRY_GEOMETRYREADRESULT_H
