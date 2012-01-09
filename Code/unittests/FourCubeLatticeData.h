#ifndef HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H
#define HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H

#include <cstdlib>
#include "geometry/LatticeData.h"

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
          readResult.blocks = util::Vector3D<site_t>::Unity();

          readResult.Blocks = std::vector<hemelb::geometry::BlockReadResult>(1);

          hemelb::geometry::BlockReadResult& block = readResult.Blocks[0];
          block.Sites.resize(readResult.GetSitesPerBlock());

          site_t index = -1;
          for (site_t i = 0; i < sitesPerBlockUnit; ++i)
          {
            for (site_t j = 0; j < sitesPerBlockUnit; ++j)
            {
              for (site_t k = 0; k < sitesPerBlockUnit; ++k)
              {
                ++index;

                hemelb::geometry::SiteReadResult& site = block.Sites[index];

                bool xMin = i == 0;
                bool yMin = j == 0;
                bool zMin = k == 0;
                bool xMax = i == (sitesPerBlockUnit - 1);
                bool yMax = j == (sitesPerBlockUnit - 1);
                bool zMax = k == (sitesPerBlockUnit - 1);

                bool nearWall = xMin || xMax || yMin || yMax;

                unsigned siteData = 0;

                // Near inlet
                if (zMin)
                {
                  siteData |= hemelb::geometry::INLET_TYPE;
                  siteData |= 0 << hemelb::geometry::SiteData::BOUNDARY_ID_SHIFT;

                  // Also near wall
                  if (nearWall)
                  {
                    siteData |= hemelb::geometry::SiteData::PRESSURE_EDGE_MASK;
                  }
                }
                // Near outlet

                else if (zMax)
                {
                  siteData |= hemelb::geometry::OUTLET_TYPE;
                  siteData |= 1 << hemelb::geometry::SiteData::BOUNDARY_ID_SHIFT;

                  if (nearWall)
                  {
                    siteData |= hemelb::geometry::SiteData::PRESSURE_EDGE_MASK;
                  }
                }
                // Not near in/outlet

                else
                {
                  siteData |= hemelb::geometry::FLUID_TYPE;
                  if (nearWall)
                  {
                    siteData |= hemelb::geometry::SiteData::PRESSURE_EDGE_MASK;
                  }
                }

                site.targetProcessor = 0;

                for (unsigned int ll = 1; ll < D3Q15::NUMVECTORS; ++ll)
                {
                  site_t neighI = i + D3Q15::CX[ll];
                  site_t neighJ = j + D3Q15::CY[ll];
                  site_t neighK = k + D3Q15::CZ[ll];

                  if (neighI < 0 || neighJ < 0 || neighK < 0 || neighI >= sitesPerBlockUnit
                      || neighJ >= sitesPerBlockUnit || neighK >= sitesPerBlockUnit)
                  {
                    siteData |= 1U << (hemelb::geometry::SiteData::BOUNDARY_CONFIG_SHIFT + ll - 1);
                    site.cutDistance[ll - 1] = (double(std::rand() % 10000) / 10000.0);
                  }
                }

                site.siteData = hemelb::geometry::SiteData(siteData);
              }
            }
          }

          FourCubeLatticeData* returnable = new FourCubeLatticeData(readResult);

          // First, fiddle with the fluid site count, for tests that require this set.
          returnable->fluidSitesOnEachProcessor.resize(rankCount);
          returnable->fluidSitesOnEachProcessor[0] = sitesPerBlockUnit * sitesPerBlockUnit
              * sitesPerBlockUnit;
          for (proc_t rank = 1; rank < rankCount; ++rank)
          {
            returnable->fluidSitesOnEachProcessor[rank] = rank * 1000;
          }

          return returnable;
        }

      protected:
        FourCubeLatticeData(hemelb::geometry::GeometryReadResult& readResult) :
          hemelb::geometry::LatticeData(readResult)
        {

        }

        FourCubeLatticeData() :
          hemelb::geometry::LatticeData()
        {

        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H */
