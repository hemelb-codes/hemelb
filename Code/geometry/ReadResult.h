#ifndef HEMELB_GEOMETRY_READRESULT_H
#define HEMELB_GEOMETRY_READRESULT_H

#include <vector>
#include "units.h"
#include "constants.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
    //! Model of the data read in about a link from one site in one direction.
    struct LinkReadResult
    {
      public:
        //! Enumeration of the different intersections the link might make between the current
        //! site and the next lattice point in this direction: no intersection,
        //! intersection with a vessel wall and intersection with an inlet or outlet.
        enum IntersectionType
        {
          NO_INTERSECTION = 0,
          WALL_INTERSECTION = 1,
          INLET_INTERSECTION = 2,
          OUTLET_INTERSECTION = 3
        } type;

        //! Default constructor. Has no intersection, nonsense values for intersection distance
        //! and iolet id.
        LinkReadResult() :
            type(NO_INTERSECTION), distanceToIntersection(-1.0), ioletId(-1)
        {
        }

        //! The proportion of the lattice vector traversed before an intersection is hit.
        float distanceToIntersection;

        //! The identifier of the inlet or outlet hit along the lattice vector (if one is hit).
        int ioletId;
    };

    /***
     * Model of the data for a site, as contained within a geometry file.
     * this data will be broken up and placed in various arrays in hemelb::Geometry::LatticeData
     */
    struct SiteReadResult
    {
      public:
        //! Basic constructor for solid and fluid sites.
        SiteReadResult(bool siteIsFluid) :
            targetProcessor(siteIsFluid ?
              -1 :
              BIG_NUMBER2), isFluid(siteIsFluid)
        {
        }

        SiteReadResult(const SiteReadResult& other) :
            targetProcessor(other.targetProcessor), isFluid(other.isFluid), links(other.links)
        {

        }

        //! Processor on which to perform lattice-Boltzmann for the site.
        proc_t targetProcessor;

        //! True iff the site is fluid, i.e. it is within the geometry and we will be doing
        //! lattice-Boltzmann with it.
        bool isFluid;

        //! A vector of the link data for each direction in the lattice currently being used
        //! (NOT necessarily the same as the lattice used by the geometry file).
        std::vector<LinkReadResult> links;
    };

    /***
     * Model of the information stored for a block in a geometry file.
     * Just gives the array of sites
     */
    struct BlockReadResult
    {
      public:
        std::vector<SiteReadResult> Sites;
    };

    /***
     * Model of the information in a geometry file
     */
    struct GeometryReadResult
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
#endif // HEMELB_GEOMETRY_READRESULT_H
