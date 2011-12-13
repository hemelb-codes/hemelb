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

        class CubeGlobalLatticeData : public geometry::LatticeData::GlobalLatticeData
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
            static CubeGlobalLatticeData* Create(site_t sitesPerBlockUnit = 4, proc_t rankCount = 1)
            {
              CubeGlobalLatticeData* globalLattice = new CubeGlobalLatticeData();
              globalLattice->SetBasicDetails(1, 1, 1, sitesPerBlockUnit, 0.01, 0.0, 0.0, 0.0);

              // All sites are dealt with on the local processor.
              globalLattice->fluidSitesOnEachProcessor = new site_t[rankCount];
              globalLattice->fluidSitesOnEachProcessor[0] = sitesPerBlockUnit * sitesPerBlockUnit
                  * sitesPerBlockUnit;
              for (proc_t rank = 1; rank < rankCount; ++rank)
              {
                globalLattice->fluidSitesOnEachProcessor[rank] = rank * 1000;
              }

              geometry::BlockData* block = &globalLattice->Blocks[0];

              block->ProcessorRankForEachBlockSite
                  = new proc_t[globalLattice->GetSitesPerBlockVolumeUnit()];
              block ->site_data = new unsigned int[globalLattice->GetSitesPerBlockVolumeUnit()];
              block->wall_data
                  = new geometry::WallData[globalLattice->GetSitesPerBlockVolumeUnit()];

              site_t index = -1;
              for (site_t i = 0; i < globalLattice->GetXSiteCount(); ++i)
              {
                for (site_t j = 0; j < globalLattice->GetYSiteCount(); ++j)
                {
                  for (site_t k = 0; k < globalLattice->GetZSiteCount(); ++k)
                  {
                    ++index;

                    bool xMin = i == 0;
                    bool yMin = j == 0;
                    bool zMin = k == 0;
                    bool xMax = i == (globalLattice->GetXSiteCount() - 1);
                    bool yMax = j == (globalLattice->GetYSiteCount() - 1);
                    bool zMax = k == (globalLattice->GetZSiteCount() - 1);

                    bool nearWall = xMin || xMax || yMin || yMax;

                    block->site_data[index] = 0;

                    // Near inlet
                    if (zMin)
                    {
                      block->site_data[index] |= INLET_TYPE;
                      block->site_data[index] |= 0 << BOUNDARY_ID_SHIFT;

                      // Also near wall
                      if (nearWall)
                      {
                        block->site_data[index] |= PRESSURE_EDGE_MASK;
                      }
                    }
                    // Near outlet

                    else if (zMax)
                    {
                      block->site_data[index] |= OUTLET_TYPE;
                      block->site_data[index] |= 1 << BOUNDARY_ID_SHIFT;

                      if (nearWall)
                      {
                        block->site_data[index] |= PRESSURE_EDGE_MASK;
                      }
                    }
                    // Not near in/outlet

                    else
                    {
                      block->site_data[index] |= FLUID_TYPE;
                      if (nearWall)
                      {
                        block->site_data[index] |= PRESSURE_EDGE_MASK;
                      }
                    }

                    block->ProcessorRankForEachBlockSite[index] = 0;

                    for (unsigned int ll = 1; ll < D3Q15::NUMVECTORS; ++ll)
                    {
                      if (!globalLattice->IsValidLatticeSite(i + D3Q15::CX[ll],
                                                             j + D3Q15::CY[ll],
                                                             k + D3Q15::CZ[ll]))
                      {
                        block->site_data[index] |= 1U << (BOUNDARY_CONFIG_SHIFT + ll - 1);
                        block->wall_data[index].cut_dist[ll - 1] = double(std::rand() % 10000)
                            / 10000.0;
                      }
                    }
                  }
                }
              }

              return globalLattice;
            }
        };

        /**
         * The create function makes a 4 x 4 x 4 cube of sites from (0,0,0) to (3,3,3).
         * The plane (x,y,0) is an inlet (boundary 0).
         * The plane (x,y,3) is an outlet (boundary 1).
         * The planes (0,y,z), (3,y,z), (x,0,z) and (x,3,z) are all walls.
         *
         * @return
         */
        static FourCubeLatticeData* Create(site_t rankCount = 1)
        {
          LatticeData::GlobalLatticeData* globalLattice = CubeGlobalLatticeData::Create(4,
                                                                                        rankCount);

          LatticeData::LocalLatticeData* localLattice =
              new LocalLatticeData(globalLattice->GetSitesPerBlockVolumeUnit()
                  * globalLattice->GetBlockCount());

          localLattice->SetSharedSiteCount(0);

          // Iterate through the fluid sides and assign variables as necessary.
          for (unsigned int collisionType = 0; collisionType < COLLISION_TYPES; ++collisionType)
          {
            localLattice->my_inner_collisions[collisionType] = 0;
            localLattice->my_inter_collisions[collisionType] = 0;
          }

          site_t index = -1;
          for (site_t i = 0; i < globalLattice->GetXSiteCount(); ++i)
          {
            for (site_t j = 0; j < globalLattice->GetYSiteCount(); ++j)
            {
              for (site_t k = 0; k < globalLattice->GetZSiteCount(); ++k)
              {
                ++index;

                bool xMin = i == 0;
                bool yMin = j == 0;
                bool zMin = k == 0;
                bool xMax = i == (globalLattice->GetXSiteCount() - 1);
                bool yMax = j == (globalLattice->GetYSiteCount() - 1);
                bool zMax = k == (globalLattice->GetZSiteCount() - 1);

                bool nearWall = xMin || xMax || yMin || yMax;

                int collType = -1;
                localLattice->mSiteData[index] = 0;

                // Near inlet
                if (zMin)
                {
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
                  localLattice->SetDistanceToWall(index,
                                                  globalLattice->Blocks[0].wall_data[index].cut_dist);
                }

                localLattice->my_inner_collisions[collType]++;
                localLattice->mSiteData[index] = globalLattice->Blocks[0].site_data[index];
              }
            }
          }

          return new FourCubeLatticeData(globalLattice, localLattice);
        }

        ~FourCubeLatticeData()
        {

        }

      protected:
        FourCubeLatticeData(GlobalLatticeData* globalLattice, LocalLatticeData* localLattice) :
          LatticeData(localLattice, globalLattice)
        {
        }

    };
  }
}

#endif /* HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H */
