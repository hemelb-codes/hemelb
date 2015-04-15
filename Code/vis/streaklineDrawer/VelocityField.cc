// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <cmath>
#include <vector>
#include <cassert>

#include "debug/Debugger.h"
#include "geometry/BlockTraverser.h"
#include "geometry/SiteTraverser.h"
#include "vis/streaklineDrawer/VelocityField.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      VelocityField::VelocityField(
          proc_t localRank_, std::map<proc_t, NeighbouringProcessor>& neighbouringProcessorsIn,
          const lb::MacroscopicPropertyCache& propertyCache) :
          counter(0), localRank(localRank_), neighbouringProcessors(neighbouringProcessorsIn),
              propertyCache(propertyCache)
      {
      }

      void VelocityField::BuildVelocityField(const geometry::LatticeData& latDat)
      {
        velocityField.resize(latDat.GetBlockCount());

        // Iterate over each block with some sites on this rank.
        geometry::BlockTraverser blockTraverser(latDat);
        do
        {
          const geometry::Block& block = blockTraverser.GetCurrentBlockData();

          if (block.IsEmpty())
          {
            continue;
          }

          geometry::SiteTraverser siteTraverser(latDat);
          do
          {
            // Only interested if the site lives on this rank.
            if (localRank != block.GetProcessorRankForSite(siteTraverser.GetCurrentIndex()))
            {
              continue;
            }

            // Calculate the bounds of the unit cube around the current site (within the lattice)
            const site_t startI =
                util::NumericalFunctions::max<site_t>(0,
                                                      blockTraverser.GetX()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetX() - 1);

            const site_t startJ =
                util::NumericalFunctions::max<site_t>(0,
                                                      blockTraverser.GetY()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetY() - 1);

            const site_t startK =
                util::NumericalFunctions::max<site_t>(0,
                                                      blockTraverser.GetZ()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetZ() - 1);

            const site_t endI =
                util::NumericalFunctions::min<site_t>(latDat.GetSiteDimensions().x - 1,
                                                      blockTraverser.GetX()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetX() + 1);

            const site_t endJ =
                util::NumericalFunctions::min<site_t>(latDat.GetSiteDimensions().y - 1,
                                                      blockTraverser.GetY()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetY() + 1);

            const site_t endK =
                util::NumericalFunctions::min<site_t>(latDat.GetSiteDimensions().z - 1,
                                                      blockTraverser.GetZ()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetZ() + 1);

            // Iterate over the sites in the unit cube.
            for (site_t neighbourI = startI; neighbourI <= endI; neighbourI++)
            {
              for (site_t neighbourJ = startJ; neighbourJ <= endJ; neighbourJ++)
              {
                for (site_t neighbourK = startK; neighbourK <= endK; neighbourK++)
                {
                  // Get the rank that the neighbour lives on.
                  const proc_t neigh_proc_id = latDat.GetProcIdFromGlobalCoords(util::Vector3D<
                      site_t>(neighbourI, neighbourJ, neighbourK));

                  // If we have data for it, we should initialise a block in the velocity field
                  // for the neighbour site.
                  if (neigh_proc_id == BIG_NUMBER2)
                  {
                    continue;
                  }

                  InitializeVelocityFieldBlock(latDat,
                                               util::Vector3D<site_t>(neighbourI,
                                                                      neighbourJ,
                                                                      neighbourK),
                                               neigh_proc_id);

                  // If the neighbour is on this rank, ignore it.
                  if (localRank == neigh_proc_id)
                  {
                    continue;
                  }

                  if (neighbouringProcessors.count(neigh_proc_id) == 0)
                  {
                    NeighbouringProcessor newProc(neigh_proc_id);
                    neighbouringProcessors[neigh_proc_id] = newProc;
                  }
                }
              }
            }
          }
          while (siteTraverser.TraverseOne());
        }
        while (blockTraverser.TraverseOne());

        // Iterate over the blocks, updating the value of the counter variable wherever
        // there is velocity field data.
        for (site_t block = 0; block < latDat.GetBlockCount(); block++)
        {
          if (velocityField[block].empty())
          {
            continue;
          }

          if (latDat.GetBlock(block).IsEmpty())
          {
            continue;
          }

          // Update the site id on each velocity field unit as required.
          for (site_t localSiteId = 0; localSiteId < latDat.GetSitesPerBlockVolumeUnit();
              localSiteId++)
          {
            velocityField[block][localSiteId].site_id =
                latDat.GetBlock(block).GetLocalContiguousIndexForSite(localSiteId);
          }
        }
      }

      bool VelocityField::BlockContainsData(size_t blockNumber) const
      {
        return !velocityField[blockNumber].empty();
      }

      VelocitySiteData& VelocityField::GetSiteData(site_t blockNumber, site_t siteNumber)
      {
        return velocityField[blockNumber][siteNumber];
      }

      // Returns the velocity site data for a given index, or nullptr if the index isn't valid / has
      // no data.
      VelocitySiteData* VelocityField::GetVelocitySiteData(const geometry::LatticeData& latDat,
                                                           const util::Vector3D<site_t>& location)
      {
        if (!latDat.IsValidLatticeSite(location))
        {
          return nullptr;
        }

        util::Vector3D<site_t> blockCoords, siteCoords;
        latDat.GetBlockAndLocalSiteCoords(location, blockCoords, siteCoords);

        site_t block_id = latDat.GetBlockIdFromBlockCoords(blockCoords);

        if (!BlockContainsData(static_cast<size_t>(block_id)))
        {
          return nullptr;
        }

        site_t site_id = latDat.GetLocalSiteIdFromLocalSiteCoords(siteCoords);

        return &GetSiteData(block_id, site_id);
      }

      // Function to initialise the velocity field at given coordinates.
      void VelocityField::InitializeVelocityFieldBlock(const geometry::LatticeData& latDat,
                                                       const util::Vector3D<site_t> location,
                                                       const proc_t proc_id)
      {
        util::Vector3D<site_t> blockCoords, siteCoords;
        latDat.GetBlockAndLocalSiteCoords(location, blockCoords, siteCoords);

        site_t blockId = latDat.GetBlockIdFromBlockCoords(blockCoords);

        if (!BlockContainsData(blockId))
        {
          for (site_t localSiteId = 0; localSiteId < latDat.GetSitesPerBlockVolumeUnit();
              ++localSiteId)
          {
            velocityField[blockId].push_back(VelocitySiteData());
          }
        }

        site_t localSiteId = latDat.GetLocalSiteIdFromLocalSiteCoords(siteCoords);
        velocityField[blockId][localSiteId].proc_id = proc_id;
      }

      // Populate the matrix v with all the velocity field data at each index.
      // Returns true if this the area resides entirely on this core.
      void VelocityField::GetVelocityFieldAroundPoint(
          const util::Vector3D<site_t> location, const geometry::LatticeData& latDat,
          util::Vector3D<float> localVelocityField[2][2][2])
      {
        for (int unitGridI = 0; unitGridI <= 1; ++unitGridI)
        {
          site_t neighbourI = location.x + unitGridI;

          for (int unitGridJ = 0; unitGridJ <= 1; ++unitGridJ)
          {
            site_t neighbourJ = location.y + unitGridJ;

            for (int unitGridK = 0; unitGridK <= 1; ++unitGridK)
            {
              site_t neighbourK = location.z + unitGridK;

              util::Vector3D<site_t> neighbour(neighbourI, neighbourJ, neighbourK);

              if (!latDat.IsValidLatticeSite(neighbour))
              {
                // it is a solid site and the velocity is
                // assumed to be zero
                localVelocityField[unitGridI][unitGridJ][unitGridK] = util::Vector3D<float>::Zero();
                continue;
              }

              VelocitySiteData *vel_site_data_p =
                  GetVelocitySiteData(latDat,
                                      util::Vector3D<site_t>(neighbourI, neighbourJ, neighbourK));

              if (vel_site_data_p == nullptr || vel_site_data_p->proc_id == -1)
              {
                // it is a solid site and the velocity is
                // assumed to be zero
                localVelocityField[unitGridI][unitGridJ][unitGridK] = util::Vector3D<float>::Zero();
                continue;
              }

              if (vel_site_data_p->counter != counter)
              {
                UpdateLocalField(vel_site_data_p, latDat);
              }

              localVelocityField[unitGridI][unitGridJ][unitGridK] = vel_site_data_p->velocity;
            }
          }
        }
      }

      bool VelocityField::NeededFromNeighbour(const util::Vector3D<site_t> location,
                                              const geometry::LatticeData& latDat,
                                              proc_t* sourceProcessor)
      {
        if (!latDat.IsValidLatticeSite(location))
        {
          return false;
        }

        VelocitySiteData *vel_site_data_p = GetVelocitySiteData(latDat, location);

        if (vel_site_data_p == nullptr || vel_site_data_p->proc_id == -1
            || vel_site_data_p->proc_id == localRank || vel_site_data_p->counter == counter)
        {
          return false;
        }

        vel_site_data_p->counter = counter;
        *sourceProcessor = vel_site_data_p->proc_id;

        return true;
      }

      void VelocityField::UpdateLocalField(const util::Vector3D<site_t>& position,
                                           const geometry::LatticeData& latDat)
      {
        VelocitySiteData *localVelocitySiteData = GetVelocitySiteData(latDat, position);

        if (log::Logger::ShouldDisplay<log::Debug>())
        {
          if (localRank != localVelocitySiteData->proc_id)
          {
            log::Logger::Log<log::Warning, log::OnePerCore>("Got a request for velocity data "
                                                            "that actually seems to be on rank %i",
                                                            localVelocitySiteData->proc_id);
          }
        }

        UpdateLocalField(localVelocitySiteData, latDat);
      }

      void VelocityField::UpdateLocalField(VelocitySiteData* localVelocitySiteData,
                                           const geometry::LatticeData& latDat)
      {
        // the local counter is set equal to the global one
        // and the local velocity is calculated
        localVelocitySiteData->counter = counter;

        const util::Vector3D<distribn_t>& velocity =
            propertyCache.velocityCache.Get(localVelocitySiteData->site_id);
        localVelocitySiteData->velocity.x = (float) velocity.x;
        localVelocitySiteData->velocity.y = (float) velocity.y;
        localVelocitySiteData->velocity.z = (float) velocity.z;
      }

      // Interpolates a velocity field to get the velocity at the position of a particle.
      util::Vector3D<float> VelocityField::InterpolateVelocityForPoint(
          const util::Vector3D<float> point,
          const util::Vector3D<float> localVelocityField[2][2][2]) const
      {
        float dummy;

        // Get the fractional parts of each of x, y, z
        float dx = modff(point.x, &dummy);
        float dy = modff(point.y, &dummy);
        float dz = modff(point.z, &dummy);

        util::Vector3D<float> yInterpolatedVelocity[2];

        for (int unitX = 0; unitX <= 1; unitX++)
        {
          util::Vector3D<float> zInterpolatedVelocityForThisX[2];

          for (int unitY = 0; unitY <= 1; unitY++)
          {
            zInterpolatedVelocityForThisX[unitY] = localVelocityField[unitX][unitY][0] * (1.F - dz)
                + localVelocityField[unitX][unitY][1] * dz;
          }

          yInterpolatedVelocity[unitX] = zInterpolatedVelocityForThisX[0] * (1.F - dy)
              + zInterpolatedVelocityForThisX[1] * dy;
        }

        return yInterpolatedVelocity[0] * (1.F - dx) + yInterpolatedVelocity[1] * dx;
      }

      void VelocityField::InvalidateAllCalculatedVelocities()
      {
        ++counter;
      }

    }
  }
}
