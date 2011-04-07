//#include <stdlib.h>
#include <math.h>
#include <vector>

#include "util/utilityFunctions.h"

#include "vis/StreaklineDrawer.h"
#include "vis/Control.h"
#include "vis/ColPixel.h"

namespace hemelb
{
  namespace vis
  {

    // TODO the streaker needs to be fixed. The lines drawn differ depending on how many proccesors
    // are used to do the visualisation. This is a bug.

    // Function to initialise the velocity field at given coordinates.
    void StreaklineDrawer::initializeVelFieldBlock(const geometry::LatticeData* iLatDat,
                                                   unsigned int site_i,
                                                   unsigned int site_j,
                                                   unsigned int site_k,
                                                   int proc_id)
    {
      unsigned int i = site_i >> iLatDat->GetLog2BlockSize();
      unsigned int j = site_j >> iLatDat->GetLog2BlockSize();
      unsigned int k = site_k >> iLatDat->GetLog2BlockSize();

      unsigned int block_id = iLatDat->GetBlockIdFromBlockCoords(i, j, k);
      if (velocity_field[block_id].vel_site_data == NULL)
      {
        velocity_field[block_id].vel_site_data
            = new VelSiteData[iLatDat->GetSitesPerBlockVolumeUnit()];

        for (unsigned int site_id = 0; site_id < iLatDat->GetSitesPerBlockVolumeUnit(); site_id++)
        {
          velocity_field[block_id].vel_site_data[site_id].proc_id = -1;
          velocity_field[block_id].vel_site_data[site_id].counter = 0;
        }
      }

      unsigned int ii = site_i - (i << iLatDat->GetLog2BlockSize());
      unsigned int jj = site_j - (j << iLatDat->GetLog2BlockSize());
      unsigned int kk = site_k - (k << iLatDat->GetLog2BlockSize());

      unsigned int site_id = ( ( (ii << iLatDat->GetLog2BlockSize()) + jj)
          << iLatDat->GetLog2BlockSize()) + kk;
      velocity_field[block_id].vel_site_data[site_id].proc_id = proc_id;
    }

    // Returns the velocity site data for a given index, or NULL if the index isn't valid / has
    // no data.
    StreaklineDrawer::VelSiteData *StreaklineDrawer::velSiteDataPointer(geometry::LatticeData* iLatDat,
                                                                        unsigned int site_i,
                                                                        unsigned int site_j,
                                                                        unsigned int site_k)
    {
      if (site_i >= iLatDat->GetXSiteCount() || site_j >= iLatDat->GetYSiteCount() || site_k
          >= iLatDat->GetZSiteCount())
      {
        return NULL;
      }
      unsigned int i = site_i >> iLatDat->GetLog2BlockSize();
      unsigned int j = site_j >> iLatDat->GetLog2BlockSize();
      unsigned int k = site_k >> iLatDat->GetLog2BlockSize();

      unsigned int block_id = iLatDat->GetBlockIdFromBlockCoords(i, j, k);

      if (velocity_field[block_id].vel_site_data == NULL)
      {
        return NULL;
      }
      unsigned int ii = site_i - (i << iLatDat->GetLog2BlockSize());
      unsigned int jj = site_j - (j << iLatDat->GetLog2BlockSize());
      unsigned int kk = site_k - (k << iLatDat->GetLog2BlockSize());

      unsigned int site_id = ( ( (ii << iLatDat->GetLog2BlockSize()) + jj)
          << iLatDat->GetLog2BlockSize()) + kk;

      return &velocity_field[block_id].vel_site_data[site_id];
    }

    // Interpolates a velocity field to get the velocity at the position of a particle.
    void StreaklineDrawer::particleVelocity(Particle *particle_p,
                                            float v[2][2][2][3],
                                            float interp_v[3])
    {
      float dx, dy, dz;
      float v_00z, v_01z, v_10z, v_11z, v_0y, v_1y;

      dx = particle_p->x - (int) particle_p->x;
      dy = particle_p->y - (int) particle_p->y;
      dz = particle_p->z - (int) particle_p->z;

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

    // Create a particle with given position, velocity and inlet_id.
    void StreaklineDrawer::createParticle(float x, float y, float z, float vel, int inlet_id)
    {

      if (nParticles == particleVec.capacity())
      {
        particleVec.reserve(2 * particleVec.capacity());
      }

      particleVec[nParticles].x = x;
      particleVec[nParticles].y = y;
      particleVec[nParticles].z = z;
      particleVec[nParticles].vel = vel;
      particleVec[nParticles].inlet_id = inlet_id;
      ++nParticles;
    }

    // Delete the particle at given index. Do something a bit budget to ensure that
    // the particles remain in the first <particles> elements of an array,
    void StreaklineDrawer::deleteParticle(unsigned int p_index)
    {
      if (nParticles == 0)
        return;

      nParticles--;
      if (nParticles == 0)
        return;

      // its data are replaced with those of the last particle;
      if (p_index != nParticles)
      {
        particleVec[p_index].x = particleVec[nParticles].x;
        particleVec[p_index].y = particleVec[nParticles].y;
        particleVec[p_index].z = particleVec[nParticles].z;
        particleVec[p_index].vx = particleVec[nParticles].vx;
        particleVec[p_index].vy = particleVec[nParticles].vy;
        particleVec[p_index].vz = particleVec[nParticles].vz;
        particleVec[p_index].vel = particleVec[nParticles].vel;
        particleVec[p_index].inlet_id = particleVec[nParticles].inlet_id;
      }

    }

    // Create seed particles to begin the streaklines.
    void StreaklineDrawer::createSeedParticles()
    {
      for (unsigned int n = 0; n < nParticleSeeds; n++)
      {
        createParticle(particleSeedVec[n].x, particleSeedVec[n].y, particleSeedVec[n].z, 0.0F,
                       particleSeedVec[n].inlet_id);
      }
    }

    // Populate the matrix v with all the velocity field data at each index.
    void StreaklineDrawer::localVelField(int p_index,
                                         float v[2][2][2][3],
                                         int *is_interior,
                                         geometry::LatticeData* iLatDat)
    {
      unsigned int site_i = (unsigned int) particleVec[p_index].x;
      unsigned int site_j = (unsigned int) particleVec[p_index].y;
      unsigned int site_k = (unsigned int) particleVec[p_index].z;

      *is_interior = 1;

      for (unsigned int i = 0; i < 2; i++)
      {
        unsigned int neigh_i = site_i + i;

        for (unsigned int j = 0; j < 2; j++)
        {
          unsigned int neigh_j = site_j + j;

          for (unsigned int k = 0; k < 2; k++)
          {
            unsigned int neigh_k = site_k + k;

            VelSiteData *vel_site_data_p = velSiteDataPointer(iLatDat, neigh_i, neigh_j, neigh_k);

            if (vel_site_data_p == NULL || vel_site_data_p->proc_id == -1)
            {
              // it is a solid site and the velocity is
              // assumed to be zero
              v[i][j][k][0] = v[i][j][k][1] = v[i][j][k][2] = 0.0F;
              continue;
            }
            if ((int) mNetworkTopology->GetLocalRank() != vel_site_data_p->proc_id)
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
            else if ((int) mNetworkTopology->GetLocalRank() == vel_site_data_p->proc_id)
            {
              // the local counter is set equal to the global one
              // and the local velocity is calculated
              vel_site_data_p->counter = counter;
              double density, vx, vy, vz;

              D3Q15::CalculateDensityAndVelocity(iLatDat->GetFOld(vel_site_data_p->site_id
                  * D3Q15::NUMVECTORS), density, vx, vy, vz);

              v[i][j][k][0] = vel_site_data_p->vx = vx / density;
              v[i][j][k][1] = vel_site_data_p->vy = vy / density;
              v[i][j][k][2] = vel_site_data_p->vz = vz / density;
            }
            else
            {
              vel_site_data_p->counter = counter;

              unsigned int m = from_proc_id_to_neigh_proc_index[vel_site_data_p->proc_id];

              mNeighProcs[m].s_to_send[3 * mNeighProcs[m].send_vs + 0] = neigh_i;
              mNeighProcs[m].s_to_send[3 * mNeighProcs[m].send_vs + 1] = neigh_j;
              mNeighProcs[m].s_to_send[3 * mNeighProcs[m].send_vs + 2] = neigh_k;
              ++mNeighProcs[m].send_vs;
            }
          }
        }
      }
    }

    // Constructor, populating fields from lattice data objects.
    StreaklineDrawer::StreaklineDrawer(const topology::NetworkTopology * iNetworkTopology,
                                       geometry::LatticeData* iLatDat,
                                       Screen* iScreen,
                                       Viewpoint* iViewpoint,
                                       VisSettings* iVisSettings) :
      mScreen(iScreen), mViewpoint(iViewpoint), mVisSettings(iVisSettings)
    {
      mNetworkTopology = iNetworkTopology;

      particleVec.reserve(10000);
      nParticles = 0;
      particleSeedVec.reserve(100);
      nParticleSeeds = 0;

      num_blocks = iLatDat->GetBlockCount();
      velocity_field = new VelocityField[iLatDat->GetBlockCount()];

      for (unsigned int n = 0; n < iLatDat->GetBlockCount(); n++)
      {
        velocity_field[n].vel_site_data = NULL;
      }

      counter = 1;
      unsigned int inlet_sites = 0;
      unsigned int n = 0;

      for (unsigned int i = 0; i < iLatDat->GetXSiteCount(); i += iLatDat->GetBlockSize())
      {
        for (unsigned int j = 0; j < iLatDat->GetYSiteCount(); j += iLatDat->GetBlockSize())
        {
          for (unsigned int k = 0; k < iLatDat->GetZSiteCount(); k += iLatDat->GetBlockSize())
          {
            geometry::LatticeData::BlockData* lBlock = iLatDat->GetBlock(n);

            ++n;

            if (lBlock->site_data == NULL)
            {
              continue;
            }

            int m = -1;

            for (unsigned int site_i = i; site_i < i + iLatDat->GetBlockSize(); site_i++)
            {
              for (unsigned int site_j = j; site_j < j + iLatDat->GetBlockSize(); site_j++)
              {
                for (unsigned int site_k = k; site_k < k + iLatDat->GetBlockSize(); site_k++)
                {

                  m++;
                  if ((int) mNetworkTopology->GetLocalRank()
                      != lBlock->ProcessorRankForEachBlockSite[m])
                    continue;

                  const unsigned int startI = util::NumericalFunctions::max<int>(0, (int) site_i
                      - 1);
                  const unsigned int startJ = util::NumericalFunctions::max<int>(0, (int) site_j
                      - 1);
                  const unsigned int startK = util::NumericalFunctions::max<int>(0, (int) site_k
                      - 1);

                  const unsigned int endI =
                      util::NumericalFunctions::min<unsigned int>(iLatDat->GetXSiteCount() - 1,
                                                                  site_i + 1);
                  const unsigned int endJ =
                      util::NumericalFunctions::min<unsigned int>(iLatDat->GetYSiteCount() - 1,
                                                                  site_j + 1);
                  const unsigned int endK =
                      util::NumericalFunctions::min<unsigned int>(iLatDat->GetZSiteCount() - 1,
                                                                  site_k + 1);

                  for (unsigned int neigh_i = startI; neigh_i <= endI; neigh_i++)
                  {
                    for (unsigned int neigh_j = startJ; neigh_j <= endJ; neigh_j++)
                    {
                      for (unsigned int neigh_k = startK; neigh_k <= endK; neigh_k++)
                      {
                        const int *neigh_proc_id = iLatDat->GetProcIdFromGlobalCoords(neigh_i,
                                                                                      neigh_j,
                                                                                      neigh_k);

                        if (neigh_proc_id == NULL || *neigh_proc_id == BIG_NUMBER2)
                        {
                          continue;
                        }

                        initializeVelFieldBlock(iLatDat, neigh_i, neigh_j, neigh_k, *neigh_proc_id);

                        if ((int) mNetworkTopology->GetLocalRank() == *neigh_proc_id)
                        {
                          continue;
                        }

                        VelSiteData* vel_site_data_p = velSiteDataPointer(iLatDat, neigh_i,
                                                                          neigh_j, neigh_k);

                        if (vel_site_data_p->counter == counter)
                        {
                          continue;
                        }

                        vel_site_data_p->counter = counter;

                        bool seenSelf = false;
                        for (unsigned int mm = 0; mm < mNeighProcs.size() && !seenSelf; mm++)
                        {
                          if (*neigh_proc_id == (int) mNeighProcs[mm].id)
                          {
                            seenSelf = true;
                            ++mNeighProcs[mm].send_vs;
                          }
                        }
                        if (seenSelf)
                        {
                          continue;
                        }

                        NeighProc lNew;

                        lNew.id = *neigh_proc_id;
                        lNew.send_vs = 1;
                        mNeighProcs.push_back(lNew);
                      }
                    }
                  }

                  int lSiteIndex = lBlock->site_data[m];

                  // if the lattice site is not an inlet
                  if (iLatDat->GetSiteType(lSiteIndex) != geometry::LatticeData::INLET_TYPE)
                  {
                    continue;
                  }
                  ++inlet_sites;

                  if (inlet_sites % 50 != 0)
                  {
                    continue;
                  }

                  if (nParticleSeeds == particleSeedVec.capacity())
                  {
                    particleSeedVec.reserve(2 * particleSeedVec.capacity());
                  }
                  particleSeedVec[nParticleSeeds].x = (float) site_i;
                  particleSeedVec[nParticleSeeds].y = (float) site_j;
                  particleSeedVec[nParticleSeeds].z = (float) site_k;
                  particleSeedVec[nParticleSeeds].inlet_id = iLatDat->GetBoundaryId(lSiteIndex);
                  ++nParticleSeeds;
                }
              }
            }
          }
        }
      }

      shared_vs = 0;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        shared_vs += mNeighProcs[m].send_vs;
      }
      if (shared_vs > 0)
      {
        s_to_send = new unsigned int[3 * shared_vs];
        s_to_recv = new unsigned int[3 * shared_vs];

        v_to_send = new float[3 * shared_vs];
        v_to_recv = new float[3 * shared_vs];
      }
      shared_vs = 0;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m].s_to_send = &s_to_send[shared_vs * 3];
        mNeighProcs[m].s_to_recv = &s_to_recv[shared_vs * 3];

        mNeighProcs[m].v_to_send = &v_to_send[shared_vs * 3];
        mNeighProcs[m].v_to_recv = &v_to_recv[shared_vs * 3];

        shared_vs += mNeighProcs[m].send_vs;

        mNeighProcs[m].send_vs = 0;
      }

      particles_to_send_max = 1000;
      particles_to_recv_max = 1000;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m].p_to_send.reserve(5 * particles_to_send_max);
        mNeighProcs[m].p_to_recv.reserve(5 * particles_to_recv_max);
      }

      req = new MPI_Request[2 * iNetworkTopology->GetProcessorCount()];

      from_proc_id_to_neigh_proc_index = new unsigned int[iNetworkTopology->GetProcessorCount()];

      for (unsigned int m = 0; m < iNetworkTopology->GetProcessorCount(); m++)
      {
        from_proc_id_to_neigh_proc_index[m] = -1;
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        from_proc_id_to_neigh_proc_index[mNeighProcs[m].id] = m;
      }

      counter = 0;

      for (n = 0; n < iLatDat->GetBlockCount(); n++)
      {
        if (velocity_field[n].vel_site_data == NULL)
          continue;

        for (unsigned int m = 0; m < iLatDat->GetSitesPerBlockVolumeUnit(); m++)
        {
          velocity_field[n].vel_site_data[m].counter = counter;
        }
        if (iLatDat->GetBlock(n)->site_data == NULL)
          continue;

        for (unsigned int m = 0; m < iLatDat->GetSitesPerBlockVolumeUnit(); m++)
        {
          velocity_field[n].vel_site_data[m].site_id = iLatDat->GetBlock(n)->site_data[m];
        }
      }
      procs = iNetworkTopology->GetProcessorCount();
    }

    // Reset the streakline drawer.
    void StreaklineDrawer::Restart()
    {
      nParticles = 0;
    }

    // Communicate site ids to other processors.
    void StreaklineDrawer::communicateSiteIds()
    {
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Irecv(&mNeighProcs[m].recv_vs, 1, MPI_UNSIGNED, mNeighProcs[m].id, 30, MPI_COMM_WORLD,
                  &req[procs + mNeighProcs[m].id]);
        MPI_Isend(&mNeighProcs[m].send_vs, 1, MPI_UNSIGNED, mNeighProcs[m].id, 30, MPI_COMM_WORLD,
                  &req[mNeighProcs[m].id]);
        MPI_Wait(&req[mNeighProcs[m].id], status);
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Wait(&req[procs + mNeighProcs[m].id], status);

        if (mNeighProcs[m].recv_vs > 0)
        {
          MPI_Irecv(mNeighProcs[m].s_to_recv, mNeighProcs[m].recv_vs * 3, MPI_UNSIGNED,
                    mNeighProcs[m].id, 40, MPI_COMM_WORLD, &req[procs + mNeighProcs[m].id]);
        }
        if (mNeighProcs[m].send_vs > 0)
        {
          MPI_Isend(mNeighProcs[m].s_to_send, mNeighProcs[m].send_vs * 3, MPI_UNSIGNED,
                    mNeighProcs[m].id, 40, MPI_COMM_WORLD, &req[mNeighProcs[m].id]);

          MPI_Wait(&req[mNeighProcs[m].id], status);
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m].recv_vs > 0)
        {
          MPI_Wait(&req[procs + mNeighProcs[m].id], status);
        }
      }
    }

    // Communicate velocities to other processors.
    void StreaklineDrawer::communicateVelocities(geometry::LatticeData* iLatDat)
    {
      int site_i, site_j, site_k;
      int neigh_i, neigh_j, neigh_k;

      VelSiteData *vel_site_data_p;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m].send_vs > 0)
        {
          MPI_Irecv(mNeighProcs[m].v_to_recv, mNeighProcs[m].send_vs * 3, MPI_FLOAT,
                    mNeighProcs[m].id, 30, MPI_COMM_WORLD, &req[procs + mNeighProcs[m].id]);
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        for (unsigned int n = 0; n < mNeighProcs[m].recv_vs; n++)
        {
          site_i = mNeighProcs[m].s_to_recv[3 * n + 0];
          site_j = mNeighProcs[m].s_to_recv[3 * n + 1];
          site_k = mNeighProcs[m].s_to_recv[3 * n + 2];

          vel_site_data_p = velSiteDataPointer(iLatDat, site_i, site_j, site_k);

          if (vel_site_data_p != NULL)
          {
            mNeighProcs[m].v_to_send[3 * n + 0] = vel_site_data_p->vx;
            mNeighProcs[m].v_to_send[3 * n + 1] = vel_site_data_p->vy;
            mNeighProcs[m].v_to_send[3 * n + 2] = vel_site_data_p->vz;
          }
          else
          {
            mNeighProcs[m].v_to_send[3 * n + 0] = 0.;
            mNeighProcs[m].v_to_send[3 * n + 1] = 0.;
            mNeighProcs[m].v_to_send[3 * n + 2] = 0.;
          }
        }
        if (mNeighProcs[m].recv_vs > 0)
        {
          MPI_Isend(mNeighProcs[m].v_to_send, mNeighProcs[m].recv_vs * 3, MPI_FLOAT,
                    mNeighProcs[m].id, 30, MPI_COMM_WORLD, &req[mNeighProcs[m].id]);

          MPI_Wait(&req[mNeighProcs[m].id], status);
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m].send_vs <= 0)
          continue;

        MPI_Wait(&req[procs + mNeighProcs[m].id], status);

        for (unsigned int n = 0; n < mNeighProcs[m].send_vs; n++)
        {
          neigh_i = mNeighProcs[m].s_to_send[3 * n + 0];
          neigh_j = mNeighProcs[m].s_to_send[3 * n + 1];
          neigh_k = mNeighProcs[m].s_to_send[3 * n + 2];

          vel_site_data_p = velSiteDataPointer(iLatDat, neigh_i, neigh_j, neigh_k);

          if (vel_site_data_p != NULL)
          {
            vel_site_data_p->vx = mNeighProcs[m].v_to_recv[3 * n + 0];
            vel_site_data_p->vy = mNeighProcs[m].v_to_recv[3 * n + 1];
            vel_site_data_p->vz = mNeighProcs[m].v_to_recv[3 * n + 2];
          }
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m].send_vs = 0;
      }
    }

    // Update the velocity field.
    void StreaklineDrawer::updateVelField(int stage_id, geometry::LatticeData* iLatDat)
    {
      unsigned int particles_temp = nParticles;

      for (int n = (int) (particles_temp - 1); n >= 0; n--)
      {
        float v[2][2][2][3];
        int is_interior;
        localVelField(n, v, &is_interior, iLatDat);

        if (stage_id == 0 && !is_interior)
        {
          continue;
        }

        float interp_v[3];

        particleVelocity(&particleVec[n], v, interp_v);
        float vel = interp_v[0] * interp_v[0] + interp_v[1] * interp_v[1] + interp_v[2]
            * interp_v[2];

        if (vel > 1.0F)
        {
          particleVec[n].vel = 1.0F;
          particleVec[n].vx = interp_v[0] / sqrtf(vel);
          particleVec[n].vy = interp_v[1] / sqrtf(vel);
          particleVec[n].vz = interp_v[2] / sqrtf(vel);

        }
        else if (vel > 1.0e-8)
        {
          particleVec[n].vel = sqrtf(vel);
          particleVec[n].vx = interp_v[0];
          particleVec[n].vy = interp_v[1];
          particleVec[n].vz = interp_v[2];

        }
        else
        {
          deleteParticle(n);
        }
      }
    }

    // Update the particles.
    void StreaklineDrawer::updateParticles()
    {
      for (unsigned int n = 0; n < nParticles; n++)
      {
        // particle coords updating (dt = 1)
        particleVec[n].x += particleVec[n].vx;
        particleVec[n].y += particleVec[n].vy;
        particleVec[n].z += particleVec[n].vz;
      }
    }

    // Communicate that particles current state to other processors.
    void StreaklineDrawer::communicateParticles(geometry::LatticeData* iLatDat)
    {
      MPI_Status status;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Irecv(&mNeighProcs[m].recv_ps, 1, MPI_UNSIGNED, mNeighProcs[m].id, 30, MPI_COMM_WORLD,
                  &req[procs + mNeighProcs[m].id]);
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m].send_ps = 0;
      }

      unsigned int particles_temp = nParticles;

      for (int n = (int) (particles_temp - 1); n >= 0; n--)
      {
        unsigned int site_i = (unsigned int) particleVec[n].x;
        unsigned int site_j = (unsigned int) particleVec[n].y;
        unsigned int site_k = (unsigned int) particleVec[n].z;

        VelSiteData *vel_site_data_p = velSiteDataPointer(iLatDat, site_i, site_j, site_k);

        if (vel_site_data_p == NULL || (int) mNetworkTopology->GetLocalRank()
            == vel_site_data_p->proc_id || vel_site_data_p->proc_id == -1)
        {
          continue;
        }
        int m = from_proc_id_to_neigh_proc_index[vel_site_data_p->proc_id];

        if (mNeighProcs[m].send_ps == particles_to_send_max)
        {
          particles_to_send_max *= 2;
          mNeighProcs[m].p_to_send.reserve(5 * particles_to_send_max);
        }

        mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 0] = particleVec[n].x;
        mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 1] = particleVec[n].y;
        mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 2] = particleVec[n].z;
        mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 3] = particleVec[n].vel;
        mNeighProcs[m].p_to_send[5 * mNeighProcs[m].send_ps + 4] = particleVec[n].inlet_id + 0.1;
        ++mNeighProcs[m].send_ps;

        deleteParticle(n);
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Isend(&mNeighProcs[m].send_ps, 1, MPI_UNSIGNED, mNeighProcs[m].id, 30, MPI_COMM_WORLD,
                  &req[mNeighProcs[m].id]);

      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Wait(&req[procs + mNeighProcs[m].id], &status);
      }

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m].send_ps > 0)
        {
          MPI_Isend(&mNeighProcs[m].p_to_send[0], mNeighProcs[m].send_ps * 5, MPI_FLOAT,
                    mNeighProcs[m].id, 40, MPI_COMM_WORLD, &req[mNeighProcs[m].id]);

          MPI_Wait(&req[mNeighProcs[m].id], &status);
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m].recv_ps > 0)
        {
          if (mNeighProcs[m].recv_ps > particles_to_recv_max)
          {
            particles_to_recv_max *= 2;
            particles_to_recv_max = util::NumericalFunctions::max(particles_to_recv_max,
                                                                  mNeighProcs[m].recv_ps);
            mNeighProcs[m].p_to_recv.reserve(5 * particles_to_recv_max);
          }
          MPI_Irecv(&mNeighProcs[m].p_to_recv[0], mNeighProcs[m].recv_ps * 5, MPI_FLOAT,
                    mNeighProcs[m].id, 40, MPI_COMM_WORLD, &req[procs + mNeighProcs[m].id]);
          MPI_Wait(&req[procs + mNeighProcs[m].id], &status);

          for (unsigned int n = 0; n < mNeighProcs[m].recv_ps; n++)
          {
            createParticle(mNeighProcs[m].p_to_recv[5 * n + 0],
                           mNeighProcs[m].p_to_recv[5 * n + 1],
                           mNeighProcs[m].p_to_recv[5 * n + 2],
                           mNeighProcs[m].p_to_recv[5 * n + 3], (int) mNeighProcs[m].p_to_recv[5
                               * n + 4]);
          }
        }
      }
    }

    // Render the streaklines
    void StreaklineDrawer::render(geometry::LatticeData* iLatDat)
    {
      int pixels_x = mScreen->GetPixelsX();
      int pixels_y = mScreen->GetPixelsY();

      for (unsigned int n = 0; n < nParticles; n++)
      {
        float p1[3], p2[3];
        p1[0] = particleVec[n].x - float (iLatDat->GetXSiteCount() >> 1);
        p1[1] = particleVec[n].y - float (iLatDat->GetYSiteCount() >> 1);
        p1[2] = particleVec[n].z - float (iLatDat->GetZSiteCount() >> 1);

        int x[2];
        mViewpoint->Project(p1, p2);
        mScreen->Transform<int> (p2, x);

        if (! (x[0] < 0 || x[0] >= pixels_x || x[1] < 0 || x[1] >= pixels_y))
        {
          ColPixel col_pixel(x[0], x[1], particleVec[n].vel, p2[2], particleVec[n].inlet_id);
          mScreen->AddPixel(&col_pixel, mVisSettings);
        }
      }
    }

    // Draw streaklines
    void StreaklineDrawer::StreakLines(int time_steps,
                                       int time_steps_per_cycle,
                                       geometry::LatticeData* iLatDat)
    {
      unsigned int particle_creation_period =
          util::NumericalFunctions::max<unsigned int>(1, (unsigned int) (time_steps_per_cycle
              / 5000.0F));

      if (time_steps % (int) (0.5 + time_steps_per_cycle
          / mVisSettings->streaklines_per_pulsatile_period) <= (mVisSettings->streakline_length
          / 100.0F) * (time_steps_per_cycle / mVisSettings->streaklines_per_pulsatile_period)
          && time_steps % particle_creation_period == 0)
      {
        createSeedParticles();
      }

      ++counter;

      updateVelField(0, iLatDat);
      communicateSiteIds();
      communicateVelocities(iLatDat);
      updateVelField(1, iLatDat);
      updateParticles();
      communicateParticles(iLatDat);
    }

    // Destructor
    StreaklineDrawer::~StreaklineDrawer()
    {
      delete[] from_proc_id_to_neigh_proc_index;
      delete[] req;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m].p_to_recv.clear();
        mNeighProcs[m].p_to_send.clear();
      }

      if (shared_vs > 0)
      {
        delete[] v_to_recv;
        delete[] v_to_send;

        delete[] s_to_recv;
        delete[] s_to_send;
      }

      for (unsigned int m = 0; m < num_blocks; m++)
      {
        if (velocity_field[m].vel_site_data != NULL)
        {
          delete[] velocity_field[m].vel_site_data;
        }
      }

      delete[] velocity_field;
      particleSeedVec.clear();
      particleVec.clear();
    }

  }
}
