#include <math.h>
#include <vector>
#include <cassert>

#include "geometry/BlockTraverser.h"
#include "geometry/SiteTraverser.h"
#include "vis/streaklineDrawer/StreaklineDrawer.h"
#include "vis/streaklineDrawer/VelocityField.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      VelocityField::VelocityField(std::map<proc_t, NeighbouringProcessor>& neighbouringProcessorsIn) :
        neighbouringProcessors(neighbouringProcessorsIn)
      {
        counter = 1;
      }

      void VelocityField::BuildVelocityField(const geometry::LatticeData& latDat,
                                             StreaklineDrawer* streaklineDrawer)
      {
        velocityField.resize(latDat.GetBlockCount());

        site_t inlet_sites = 0;

        geometry::BlockTraverser blockTraverser(latDat);
        do
        {
          geometry::BlockData* block = blockTraverser.GetCurrentBlockData();

          if (block->site_data == NULL)
          {
            continue;
          }

          geometry::SiteTraverser siteTraverser(latDat);
          do
          {
            if (topology::NetworkTopology::Instance()->GetLocalRank()
                != block->ProcessorRankForEachBlockSite[siteTraverser.GetCurrentIndex()])
            {
              continue;
            }

            const site_t startI = util::NumericalFunctions::max<site_t>(0, blockTraverser.GetX()
                * blockTraverser.GetBlockSize() + siteTraverser.GetX() - 1);

            const site_t startJ = util::NumericalFunctions::max<site_t>(0, blockTraverser.GetY()
                * blockTraverser.GetBlockSize() + siteTraverser.GetY() - 1);

            const site_t startK = util::NumericalFunctions::max<site_t>(0, blockTraverser.GetZ()
                * blockTraverser.GetBlockSize() + siteTraverser.GetZ() - 1);

            const site_t endI =
                util::NumericalFunctions::min<site_t>(latDat.GetXSiteCount() - 1,
                                                      blockTraverser.GetX()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetX() + 1);

            const site_t endJ =
                util::NumericalFunctions::min<site_t>(latDat.GetYSiteCount() - 1,
                                                      blockTraverser.GetY()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetY() + 1);

            const site_t endK =
                util::NumericalFunctions::min<site_t>(latDat.GetZSiteCount() - 1,
                                                      blockTraverser.GetZ()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetZ() + 1);

            for (site_t neigh_i = startI; neigh_i <= endI; neigh_i++)
            {
              for (site_t neigh_j = startJ; neigh_j <= endJ; neigh_j++)
              {
                for (site_t neigh_k = startK; neigh_k <= endK; neigh_k++)
                {
                  const proc_t* neigh_proc_id = latDat.GetProcIdFromGlobalCoords(util::Vector3D<
                      site_t>(neigh_i, neigh_j, neigh_k));

                  if (neigh_proc_id == NULL || *neigh_proc_id == BIG_NUMBER2)
                  {
                    continue;
                  }

                  InitializeVelFieldBlock(latDat,
                                          util::Vector3D<site_t>(neigh_i, neigh_j, neigh_k),
                                          *neigh_proc_id);

                  if (topology::NetworkTopology::Instance()->GetLocalRank() == *neigh_proc_id)
                  {
                    continue;
                  }

                  VelocitySiteData* vel_site_data_p = GetVelocitySiteData(latDat, util::Vector3D<
                      site_t>(neigh_i, neigh_j, neigh_k));

                  vel_site_data_p->counter = counter;

                  if (neighbouringProcessors.count(*neigh_proc_id) == 0)
                  {
                    NeighbouringProcessor newProc(*neigh_proc_id);
                    neighbouringProcessors[*neigh_proc_id] = newProc;
                  }
                }
              }
            }

            site_t siteIndex = block->site_data[siteTraverser.GetCurrentIndex()];

            // if the lattice site is not an inlet
            if (latDat.GetSiteType(siteIndex) != geometry::LatticeData::INLET_TYPE)
            {
              continue;
            }
            ++inlet_sites;

            // TODO this is a problem on multiple cores.
            if (inlet_sites % 50 != 0)
            {
              continue;
            }

            streaklineDrawer->particleSeeds.push_back(Particle(static_cast<float> (blockTraverser.GetX()
                                                                   * blockTraverser.GetBlockSize()
                                                                   + siteTraverser.GetX()),
                                                               static_cast<float> (blockTraverser.GetY()
                                                                   * blockTraverser.GetBlockSize()
                                                                   + siteTraverser.GetY()),
                                                               static_cast<float> (blockTraverser.GetZ()
                                                                   * blockTraverser.GetBlockSize()
                                                                   + siteTraverser.GetZ()),
                                                               latDat.GetBoundaryId(siteIndex)));

          }
          while (siteTraverser.TraverseOne());
        }
        while (blockTraverser.TraverseOne());

        for (site_t n = 0; n < latDat.GetBlockCount(); n++)
        {
          if (velocityField[n].empty())
          {
            continue;
          }

          for (site_t m = 0; m < latDat.GetSitesPerBlockVolumeUnit(); m++)
          {
            velocityField[n][m].counter = counter;
          }
          if (latDat.GetBlock(n)->site_data == NULL)
            continue;

          for (site_t m = 0; m < latDat.GetSitesPerBlockVolumeUnit(); m++)
          {
            velocityField[n][m].site_id = latDat.GetBlock(n)->site_data[m];
          }
        }

        counter = 0;
      }

      bool VelocityField::BlockContainsData(size_t blockNumber) const
      {
        return !velocityField[blockNumber].empty();
      }

      VelocitySiteData& VelocityField::GetSiteData(site_t blockNumber, site_t siteNumber)
      {
        return velocityField[blockNumber][siteNumber];
      }

      // Returns the velocity site data for a given index, or NULL if the index isn't valid / has
      // no data.
      VelocitySiteData* VelocityField::GetVelocitySiteData(const geometry::LatticeData& latDat,
                                                           const util::Vector3D<site_t>& location)
      {
        if (!latDat.IsValidLatticeSite(location.x, location.y, location.z))
        {
          return NULL;
        }

        site_t i = location.x >> latDat.GetLog2BlockSize();
        site_t j = location.y >> latDat.GetLog2BlockSize();
        site_t k = location.z >> latDat.GetLog2BlockSize();

        site_t block_id = latDat.GetBlockIdFromBlockCoords(i, j, k);

        if (!BlockContainsData(static_cast<size_t> (block_id)))
        {
          return NULL;
        }

        site_t ii = location.x - (i << latDat.GetLog2BlockSize());
        site_t jj = location.y - (j << latDat.GetLog2BlockSize());
        site_t kk = location.z - (k << latDat.GetLog2BlockSize());

        site_t site_id = ( ( (ii << latDat.GetLog2BlockSize()) + jj) << latDat.GetLog2BlockSize())
            + kk;

        return &GetSiteData(block_id, site_id);
      }

      // Function to initialise the velocity field at given coordinates.
      void VelocityField::InitializeVelFieldBlock(const geometry::LatticeData& latDat,
                                                  const util::Vector3D<site_t> location,
                                                  const proc_t proc_id)
      {
        site_t i = location.x >> latDat.GetLog2BlockSize();
        site_t j = location.y >> latDat.GetLog2BlockSize();
        site_t k = location.z >> latDat.GetLog2BlockSize();

        site_t block_id = latDat.GetBlockIdFromBlockCoords(i, j, k);

        if (!BlockContainsData(block_id))
        {
          velocityField[block_id]
              = std::vector<VelocitySiteData>(latDat.GetSitesPerBlockVolumeUnit());
        }

        site_t ii = location.x - (i << latDat.GetLog2BlockSize());
        site_t jj = location.y - (j << latDat.GetLog2BlockSize());
        site_t kk = location.z - (k << latDat.GetLog2BlockSize());

        site_t site_id = ( ( (ii << latDat.GetLog2BlockSize()) + jj) << latDat.GetLog2BlockSize())
            + kk;
        velocityField[block_id][site_id].proc_id = proc_id;
      }

      // Populate the matrix v with all the velocity field data at each index.
      // Returns true if this the area resides entirely on this core.
      void VelocityField::GetVelocityFieldAroundPoint(const util::Vector3D<site_t> location,
                                                      const geometry::LatticeData& latDat,
                                                      float v[2][2][2][3])
      {
        const proc_t thisRank = topology::NetworkTopology::Instance()->GetLocalRank();

        for (int unitGridI = 0; unitGridI <= 1; ++unitGridI)
        {
          site_t neighbourI = location.x + unitGridI;

          for (int unitGridJ = 0; unitGridJ <= 1; ++unitGridJ)
          {
            site_t neighbourJ = location.y + unitGridJ;

            for (int unitGridK = 0; unitGridK <= 1; ++unitGridK)
            {
              site_t neighbourK = location.z + unitGridK;

              if (!latDat.IsValidLatticeSite(neighbourI, neighbourJ, neighbourK))
              {
                // it is a solid site and the velocity is
                // assumed to be zero
                v[unitGridI][unitGridJ][unitGridK][0] = v[unitGridI][unitGridJ][unitGridK][1]
                    = v[unitGridI][unitGridJ][unitGridK][2] = 0.0F;
                continue;
              }

              VelocitySiteData *vel_site_data_p =
                  GetVelocitySiteData(latDat, util::Vector3D<site_t>(neighbourI,
                                                                     neighbourJ,
                                                                     neighbourK));

              if (vel_site_data_p->proc_id == -1)
              {
                // it is a solid site and the velocity is
                // assumed to be zero
                v[unitGridI][unitGridJ][unitGridK][0] = v[unitGridI][unitGridJ][unitGridK][1]
                    = v[unitGridI][unitGridJ][unitGridK][2] = 0.0F;
                continue;
              }

              if (vel_site_data_p->counter == counter)
              {
                // This means that the local velocity has already been
                // calculated at the current time step if the site
                // belongs to the current processor; if not, the
                // following instructions have no effect
                v[unitGridI][unitGridJ][unitGridK][0] = vel_site_data_p->vx;
                v[unitGridI][unitGridJ][unitGridK][1] = vel_site_data_p->vy;
                v[unitGridI][unitGridJ][unitGridK][2] = vel_site_data_p->vz;
              }
              else if (thisRank == vel_site_data_p->proc_id)
              {
                // the local counter is set equal to the global one
                // and the local velocity is calculated
                vel_site_data_p->counter = counter;
                distribn_t density, vx, vy, vz;

                D3Q15::CalculateDensityAndVelocity(latDat.GetFOld(vel_site_data_p->site_id
                    * D3Q15::NUMVECTORS), density, vx, vy, vz);

                v[unitGridI][unitGridJ][unitGridK][0] = vel_site_data_p->vx
                    = (float) (vx / density);
                v[unitGridI][unitGridJ][unitGridK][1] = vel_site_data_p->vy
                    = (float) (vy / density);
                v[unitGridI][unitGridJ][unitGridK][2] = vel_site_data_p->vz
                    = (float) (vz / density);
              }
              else
              {
                vel_site_data_p->counter = counter;

                neighbouringProcessors[vel_site_data_p->proc_id].AddSiteToRequestVelocityDataFor(neighbourI,
                                                                                                 neighbourJ,
                                                                                                 neighbourK);
              }
            }
          }
        }
      }

      // Interpolates a velocity field to get the velocity at the position of a particle.
      void VelocityField::InterpolateVelocityForPoint(const float x,
                                                      const float y,
                                                      const float z,
                                                      const float v[2][2][2][3],
                                                      float interp_v[3]) const
      {
        float dummy;

        // Get the fractional parts of each of x, y, z
        float dx = modff(x, &dummy);
        float dy = modff(y, &dummy);
        float dz = modff(z, &dummy);

        for (int l = 0; l < 3; l++)
        {
          float v_00z = (1.F - dz) * v[0][0][0][l] + dz * v[0][0][1][l];
          float v_01z = (1.F - dz) * v[0][1][0][l] + dz * v[0][1][1][l];
          float v_10z = (1.F - dz) * v[1][0][0][l] + dz * v[1][0][1][l];
          float v_11z = (1.F - dz) * v[1][1][0][l] + dz * v[1][1][1][l];

          float v_0y = (1.F - dy) * v_00z + dy * v_01z;
          float v_1y = (1.F - dy) * v_10z + dy * v_11z;

          interp_v[l] = (1.F - dx) * v_0y + dx * v_1y;
        }
      }

      void VelocityField::InvalidateAllCalculatedVelocities()
      {
        ++counter;
      }

    }
  }
}
