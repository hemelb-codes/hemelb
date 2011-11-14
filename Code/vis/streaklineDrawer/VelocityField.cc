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
      VelocityField::VelocityField(std::vector<NeighbouringProcessor>& neighbouringProcessorsIn) :
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
                  const proc_t* neigh_proc_id = latDat.GetProcIdFromGlobalCoords(neigh_i,
                                                                                 neigh_j,
                                                                                 neigh_k);

                  if (neigh_proc_id == NULL || *neigh_proc_id == BIG_NUMBER2)
                  {
                    continue;
                  }

                  initializeVelFieldBlock(latDat,
                                          neigh_i,
                                          neigh_j,
                                          neigh_k,
                                          *neigh_proc_id,
                                          streaklineDrawer);

                  if (topology::NetworkTopology::Instance()->GetLocalRank() == *neigh_proc_id)
                  {
                    continue;
                  }

                  VelocitySiteData* vel_site_data_p = velSiteDataPointer(latDat,
                                                                         neigh_i,
                                                                         neigh_j,
                                                                         neigh_k);

                  if (vel_site_data_p->counter == counter)
                  {
                    continue;
                  }

                  vel_site_data_p->counter = counter;

                  bool seenSelf = false;
                  for (size_t mm = 0; mm < neighbouringProcessors.size() && !seenSelf; mm++)
                  {
                    if (*neigh_proc_id == neighbouringProcessors[mm].mID)
                    {
                      seenSelf = true;
                      ++neighbouringProcessors[mm].send_vs;
                    }
                  }
                  if (seenSelf)
                  {
                    continue;
                  }

                  NeighbouringProcessor newParticle(*neigh_proc_id);

                  newParticle.send_vs = 1;
                  neighbouringProcessors.push_back(newParticle);
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

      bool VelocityField::BlockContainsData(size_t blockNumber)
      {
        return !velocityField[blockNumber].empty();
      }

      VelocitySiteData& VelocityField::GetSiteData(size_t blockNumber, size_t siteNumber)
      {
        return velocityField[blockNumber][siteNumber];
      }

      // Returns the velocity site data for a given index, or NULL if the index isn't valid / has
      // no data.
      VelocitySiteData* VelocityField::velSiteDataPointer(const geometry::LatticeData& latDat,
                                                          site_t site_i,
                                                          site_t site_j,
                                                          site_t site_k)
      {
        if (site_i >= latDat.GetXSiteCount() || site_j >= latDat.GetYSiteCount() || site_k
            >= latDat.GetZSiteCount())
        {
          return NULL;
        }
        site_t i = site_i >> latDat.GetLog2BlockSize();
        site_t j = site_j >> latDat.GetLog2BlockSize();
        site_t k = site_k >> latDat.GetLog2BlockSize();

        site_t block_id = latDat.GetBlockIdFromBlockCoords(i, j, k);

        if (!BlockContainsData(static_cast<size_t> (block_id)))
        {
          return NULL;
        }
        site_t ii = site_i - (i << latDat.GetLog2BlockSize());
        site_t jj = site_j - (j << latDat.GetLog2BlockSize());
        site_t kk = site_k - (k << latDat.GetLog2BlockSize());

        site_t site_id = ( ( (ii << latDat.GetLog2BlockSize()) + jj) << latDat.GetLog2BlockSize())
            + kk;

        return &GetSiteData(block_id, site_id);
      }

      // Function to initialise the velocity field at given coordinates.
      void VelocityField::initializeVelFieldBlock(const geometry::LatticeData& latDat,
                                                  site_t site_i,
                                                  site_t site_j,
                                                  site_t site_k,
                                                  proc_t proc_id,
                                                  StreaklineDrawer* streaklineDrawer)
      {
        site_t i = site_i >> latDat.GetLog2BlockSize();
        site_t j = site_j >> latDat.GetLog2BlockSize();
        site_t k = site_k >> latDat.GetLog2BlockSize();

        site_t block_id = latDat.GetBlockIdFromBlockCoords(i, j, k);

        if (!BlockContainsData(block_id))
        {
          velocityField[block_id]
              = std::vector<VelocitySiteData>(latDat.GetSitesPerBlockVolumeUnit());
        }

        site_t ii = site_i - (i << latDat.GetLog2BlockSize());
        site_t jj = site_j - (j << latDat.GetLog2BlockSize());
        site_t kk = site_k - (k << latDat.GetLog2BlockSize());

        site_t site_id = ( ( (ii << latDat.GetLog2BlockSize()) + jj) << latDat.GetLog2BlockSize())
            + kk;
        velocityField[block_id][site_id].proc_id = proc_id;
      }

      // Populate the matrix v with all the velocity field data at each index.
      void VelocityField::localVelField(site_t x,
                                        site_t y,
                                        site_t z,
                                        float v[2][2][2][3],
                                        int *is_interior,
                                        const geometry::LatticeData& latDat,
                                        StreaklineDrawer* streaklineDrawer)
      {
        /*
         site_t site_i = (site_t) iStreaklineDrawer->mParticleVec[p_index].x;
         site_t site_j = (site_t) iStreaklineDrawer->mParticleVec[p_index].y;
         site_t site_k = (site_t) iStreaklineDrawer->mParticleVec[p_index].z;*/

        *is_interior = 1;

        const proc_t thisRank = topology::NetworkTopology::Instance()->GetLocalRank();

        for (unsigned int i = 0; i < 2; i++)
        {
          site_t neigh_i = x + i;

          for (unsigned int j = 0; j < 2; j++)
          {
            site_t neigh_j = y + j;

            for (unsigned int k = 0; k < 2; k++)
            {
              site_t neigh_k = z + k;

              VelocitySiteData *vel_site_data_p = velSiteDataPointer(latDat,
                                                                     neigh_i,
                                                                     neigh_j,
                                                                     neigh_k);

              if (vel_site_data_p == NULL || vel_site_data_p->proc_id == -1)
              {
                // it is a solid site and the velocity is
                // assumed to be zero
                v[i][j][k][0] = v[i][j][k][1] = v[i][j][k][2] = 0.0F;
                continue;
              }
              if (thisRank != vel_site_data_p->proc_id)
              {
                *is_interior = 0;
              }
              if (vel_site_data_p->counter == counter)
              {
                // This means that the local velocity has already been
                // calculated at the current time step if the site
                // belongs to the current processor; if not, the
                // following instructions have no effect
                v[i][j][k][0] = vel_site_data_p->vx;
                v[i][j][k][1] = vel_site_data_p->vy;
                v[i][j][k][2] = vel_site_data_p->vz;
              }
              else if (thisRank == vel_site_data_p->proc_id)
              {
                // the local counter is set equal to the global one
                // and the local velocity is calculated
                vel_site_data_p->counter = counter;
                distribn_t density, vx, vy, vz;

                D3Q15::CalculateDensityAndVelocity(latDat.GetFOld(vel_site_data_p->site_id
                                                       * D3Q15::NUMVECTORS),
                                                   density,
                                                   vx,
                                                   vy,
                                                   vz);

                v[i][j][k][0] = vel_site_data_p->vx = (float) (vx / density);
                v[i][j][k][1] = vel_site_data_p->vy = (float) (vy / density);
                v[i][j][k][2] = vel_site_data_p->vz = (float) (vz / density);
              }
              else
              {
                vel_site_data_p->counter = counter;

                proc_t m =
                    streaklineDrawer->from_proc_id_to_neigh_proc_index[vel_site_data_p->proc_id];

                neighbouringProcessors[m].s_to_send[3 * neighbouringProcessors[m].send_vs + 0]
                    = neigh_i;
                neighbouringProcessors[m].s_to_send[3 * neighbouringProcessors[m].send_vs + 1]
                    = neigh_j;
                neighbouringProcessors[m].s_to_send[3 * neighbouringProcessors[m].send_vs + 2]
                    = neigh_k;
                ++ (neighbouringProcessors[m].send_vs);
              }
            }
          }
        }
      }

      // Interpolates a velocity field to get the velocity at the position of a particle.
      void VelocityField::GetVelocityAtPoint(float x,
                                             float y,
                                             float z,
                                             float v[2][2][2][3],
                                             float interp_v[3])
      {
        float v_00z, v_01z, v_10z, v_11z, v_0y, v_1y;

        float dummy;

        float dx = modff(x, &dummy);
        float dy = modff(y, &dummy);
        float dz = modff(z, &dummy);

        for (int l = 0; l < 3; l++)
        {
          v_00z = (1.F - dz) * v[0][0][0][l] + dz * v[0][0][1][l];
          v_01z = (1.F - dz) * v[0][1][0][l] + dz * v[0][1][1][l];
          v_10z = (1.F - dz) * v[1][0][0][l] + dz * v[1][0][1][l];
          v_11z = (1.F - dz) * v[1][1][0][l] + dz * v[1][1][1][l];

          v_0y = (1.F - dy) * v_00z + dy * v_01z;
          v_1y = (1.F - dy) * v_10z + dy * v_11z;

          interp_v[l] = (1.F - dx) * v_0y + dx * v_1y;
        }
      }

    }
  }
}
