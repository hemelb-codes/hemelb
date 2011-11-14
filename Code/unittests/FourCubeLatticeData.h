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
         * The constructor makes a 4 x 4 x 4 cube of sites from (0,0,0) to (3,3,3).
         * The plane (x,y,0) is an inlet (boundary 0).
         * The plane (x,y,3) is an outlet (boundary 1).
         * The planes (0,y,z), (3,y,z), (x,0,z) and (x,3,z) are all walls.
         *
         * @return
         */
        FourCubeLatticeData() :
          LatticeData()
        {
          globLatDat.SetBasicDetails(1, 1, 1, 4, 0.01, 0.0, 0.0, 0.0);

          // All sites are dealt with on the local processor.
          localLatDat.my_inner_sites = globLatDat.GetSitesPerBlockVolumeUnit()
              * globLatDat.GetBlockCount();
          localLatDat.Initialise(localLatDat.my_inner_sites);
          localLatDat.SetSharedSiteCount(0);

          geometry::BlockData* block = &globLatDat.Blocks[0];

          block->ProcessorRankForEachBlockSite
              = new proc_t[globLatDat.GetSitesPerBlockVolumeUnit()];
          block ->site_data = new unsigned int[globLatDat.GetSitesPerBlockVolumeUnit()];
          block->wall_data = new geometry::WallData[globLatDat.GetSitesPerBlockVolumeUnit()];

          // Iterate through the fluid sides and assign variables as necessary.
          for (unsigned int collisionType = 0; collisionType < COLLISION_TYPES; ++collisionType)
          {
            localLatDat.my_inner_collisions[collisionType] = 0;
            localLatDat.my_inter_collisions[collisionType] = 0;
          }

          site_t index = -1;
          for (site_t i = 0; i < globLatDat.GetXSiteCount(); ++i)
          {
            for (site_t j = 0; j < globLatDat.GetYSiteCount(); ++j)
            {
              for (site_t k = 0; k < globLatDat.GetZSiteCount(); ++k)
              {
                ++index;

                bool xMin = i == 0;
                bool yMin = j == 0;
                bool zMin = k == 0;
                bool xMax = i == (globLatDat.GetXSiteCount() - 1);
                bool yMax = j == (globLatDat.GetYSiteCount() - 1);
                bool zMax = k == (globLatDat.GetZSiteCount() - 1);

                bool nearWall = xMin || xMax || yMin || yMax;

                int collType = -1;
                localLatDat.mSiteData[index] = 0;

                // Near inlet
                if (zMin)
                {
                  localLatDat.mSiteData[index] |= INLET_TYPE;
                  localLatDat.mSiteData[index] |= 0 << BOUNDARY_ID_SHIFT;

                  // Also near wall
                  if (nearWall)
                  {
                    collType = 4;
                  }
                  else
                  {
                    collType = 2;
                  }
                }
                // Near outlet
                else if (zMax)
                {
                  localLatDat.mSiteData[index] |= OUTLET_TYPE;
                  localLatDat.mSiteData[index] |= 0 << BOUNDARY_ID_SHIFT;

                  if (nearWall)
                  {
                    collType = 5;
                  }
                  else
                  {
                    collType = 3;
                  }
                }
                // Not near in/outlet
                else
                {
                  localLatDat.mSiteData[index] |= FLUID_TYPE;
                  if (nearWall)
                  {
                    collType = 1;
                  }
                  else
                  {
                    collType = 0;
                  }
                }

                for (unsigned int ll = 1; ll < D3Q15::NUMVECTORS; ++ll)
                {
                  if (!globLatDat.IsValidLatticeSite(i + D3Q15::CX[ll],
                                                     j + D3Q15::CY[ll],
                                                     k + D3Q15::CZ[ll]))
                  {
                    localLatDat.mSiteData[index] |= 1U << (BOUNDARY_CONFIG_SHIFT + ll - 1);
                    block->wall_data[index].cut_dist[ll - 1] = double(std::rand() % 10000)
                        / 10000.0;
                  }
                  localLatDat.SetDistanceToWall(index, block->wall_data[index].cut_dist);
                }

                localLatDat.my_inner_collisions[collType]++;
                block->ProcessorRankForEachBlockSite[index] = 0;
                block->site_data[index] = (unsigned int) index;
              }
            }
          }

          InitialiseNeighbourLookup(NULL, 0, localLatDat.mSiteData);
        }

        ~FourCubeLatticeData()
        {

        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H */
