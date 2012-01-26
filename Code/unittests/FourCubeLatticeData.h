#ifndef HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H
#define HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H

#include <cstdlib>
#include "geometry/LatticeData.h"
#include "io/formats/geometry.h"

namespace hemelb
{
  namespace unittests
  {
    class FourCubeLatticeData : public geometry::LatticeData
    {
      public:

        /**
         * The create function makes a 4 x 4 x 4 cube of sites from (0,0,0) to (3,3,3).
         * The plane (x,y,0) is an inlet (boundary 0).
         * The plane (x,y,3) is an outlet (boundary 1).
         * The planes (0,y,z), (3,y,z), (x,0,z) and (x,3,z) are all walls.
         *
         * @return
         */
        static LatticeData* Create(site_t sitesPerBlockUnit = 4, proc_t rankCount = 1)
        {
          hemelb::geometry::GeometryReadResult readResult;

          readResult.voxelSize = 0.01;
          readResult.origin = util::Vector3D<distribn_t>::Zero();
          readResult.blockSize = sitesPerBlockUnit;
          readResult.blocks = util::Vector3D<site_t>::Ones();

          readResult.Blocks = std::vector < hemelb::geometry::BlockReadResult > (1);

          hemelb::geometry::BlockReadResult& block = readResult.Blocks[0];
          block.Sites.resize(readResult.GetSitesPerBlock(), geometry::SiteReadResult(false));

          site_t index = -1;
          for (site_t i = 0; i < sitesPerBlockUnit; ++i)
          {
            for (site_t j = 0; j < sitesPerBlockUnit; ++j)
            {
              for (site_t k = 0; k < sitesPerBlockUnit; ++k)
              {
                ++index;

                hemelb::geometry::SiteReadResult& site = block.Sites[index];

                site.isFluid = true;
                site.targetProcessor = 0;

                for (Direction direction = 1; direction < D3Q15::NUMVECTORS; ++direction)
                {
                  site_t neighI = i + D3Q15::CX[direction];
                  site_t neighJ = j + D3Q15::CY[direction];
                  site_t neighK = k + D3Q15::CZ[direction];

                  hemelb::geometry::LinkReadResult link;

                  float randomDistance = (float(std::rand() % 10000) / 10000.0);

                  // The inlet is by the minimal z value.
                  if (neighK < 0)
                  {
                    link.ioletId = 0;
                    link.type = geometry::LinkReadResult::INLET_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }
                  // The outlet is by the maximal z value.
                  else if (neighK >= sitesPerBlockUnit)
                  {
                    link.ioletId = 0;
                    link.type = geometry::LinkReadResult::OUTLET_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }
                  // Walls are by extremes of x and y.
                  else if (neighI < 0 || neighJ < 0 || neighI >= sitesPerBlockUnit || neighJ >= sitesPerBlockUnit)
                  {
                    link.type = geometry::LinkReadResult::WALL_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }

                  site.links.push_back(link);
                }
              }
            }
          }

          FourCubeLatticeData* returnable = new FourCubeLatticeData(readResult);

          // First, fiddle with the fluid site count, for tests that require this set.
          returnable->fluidSitesOnEachProcessor.resize(rankCount);
          returnable->fluidSitesOnEachProcessor[0] = sitesPerBlockUnit * sitesPerBlockUnit * sitesPerBlockUnit;
          for (proc_t rank = 1; rank < rankCount; ++rank)
          {
            returnable->fluidSitesOnEachProcessor[rank] = rank * 1000;
          }

          return returnable;
        }

      protected:
        FourCubeLatticeData(hemelb::geometry::GeometryReadResult& readResult) :
            hemelb::geometry::LatticeData(D3Q15::GetLatticeInfo(), readResult)
        {

        }

        FourCubeLatticeData() :
            hemelb::geometry::LatticeData(D3Q15::GetLatticeInfo())
        {

        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H */
