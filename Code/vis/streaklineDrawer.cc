//#include <stdlib.h>
#include <math.h>
#include <vector>

#include "utilityFunctions.h"

#include "vis/StreaklineDrawer.h"
#include "vis/Control.h"
#include "vis/ColPixel.h"

namespace hemelb
{
  namespace vis
  {

    // Function to initialise the velocity field at given coordinates.
    void StreaklineDrawer::initializeVelFieldBlock(lb::GlobalLatticeData &iGlobLatDat,
                                                   int site_i,
                                                   int site_j,
                                                   int site_k,
                                                   int proc_id)
    {
      if (site_i < 0 || site_i >= iGlobLatDat.SitesX || site_j < 0 || site_j
          >= iGlobLatDat.SitesY || site_k < 0 || site_k >= iGlobLatDat.SitesZ)
      {
        return;
      }
      int i = site_i >> iGlobLatDat.Log2BlockSize;
      int j = site_j >> iGlobLatDat.Log2BlockSize;
      int k = site_k >> iGlobLatDat.Log2BlockSize;

      int block_id = (i * iGlobLatDat.BlocksY + j) * iGlobLatDat.BlocksZ + k;
      int site_id;

      if (velocity_field[block_id].vel_site_data == NULL)
      {
        velocity_field[block_id].vel_site_data
            = new VelSiteData[iGlobLatDat.SitesPerBlockVolumeUnit];

        for (site_id = 0; site_id < iGlobLatDat.SitesPerBlockVolumeUnit; site_id++)
        {
          velocity_field[block_id].vel_site_data[site_id].proc_id = -1;
          velocity_field[block_id].vel_site_data[site_id].counter = 0;
        }
      }

      int ii = site_i - (i << iGlobLatDat.Log2BlockSize);
      int jj = site_j - (j << iGlobLatDat.Log2BlockSize);
      int kk = site_k - (k << iGlobLatDat.Log2BlockSize);

      site_id = ( ( (ii << iGlobLatDat.Log2BlockSize) + jj)
          << iGlobLatDat.Log2BlockSize) + kk;
      velocity_field[block_id].vel_site_data[site_id].proc_id = proc_id;
    }

    // Returns the velocity site data for a given index, or NULL if the index isn't valid / has
    // no data.
    StreaklineDrawer::VelSiteData *StreaklineDrawer::velSiteDataPointer(lb::GlobalLatticeData &iGlobLatDat,
                                                                        int site_i,
                                                                        int site_j,
                                                                        int site_k)
    {
      if (site_i < 0 || site_i >= iGlobLatDat.SitesX || site_j < 0 || site_j
          >= iGlobLatDat.SitesY || site_k < 0 || site_k >= iGlobLatDat.SitesZ)
      {
        return NULL;
      }
      int i = site_i >> iGlobLatDat.Log2BlockSize;
      int j = site_j >> iGlobLatDat.Log2BlockSize;
      int k = site_k >> iGlobLatDat.Log2BlockSize;

      int block_id = (i * iGlobLatDat.BlocksY + j) * iGlobLatDat.BlocksZ + k;

      if (velocity_field[block_id].vel_site_data == NULL)
      {
        return NULL;
      }
      int ii = site_i - (i << iGlobLatDat.Log2BlockSize);
      int jj = site_j - (j << iGlobLatDat.Log2BlockSize);
      int kk = site_k - (k << iGlobLatDat.Log2BlockSize);

      int site_id = ( ( (ii << iGlobLatDat.Log2BlockSize) + jj)
          << iGlobLatDat.Log2BlockSize) + kk;

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
    void StreaklineDrawer::createParticle(float x,
                                          float y,
                                          float z,
                                          float vel,
                                          int inlet_id)
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
        createParticle(particleSeedVec[n].x, particleSeedVec[n].y,
                       particleSeedVec[n].z, 0.0F, particleSeedVec[n].inlet_id);
      }
    }

    // Populate the matrix v with all the velocity field data at each index.
    void StreaklineDrawer::localVelField(int p_index,
                                         float v[2][2][2][3],
                                         int *is_interior,
                                         lb::GlobalLatticeData &iGlobLatDat,
                                         lb::LocalLatticeData &iLocalLatDat)
    {
      double vx, vy, vz;

      VelSiteData *vel_site_data_p;

      int site_i = int(particleVec[p_index].x);
      int site_j = int(particleVec[p_index].y);
      int site_k = int(particleVec[p_index].z);

      *is_interior = 1;

      int c1Plusc2 = 15;

      for (int i = 0; i < 2; i++)
      {
        int neigh_i = site_i + i;

        for (int j = 0; j < 2; j++)
        {
          int neigh_j = site_j + j;

          for (int k = 0; k < 2; k++)
          {
            int neigh_k = site_k + k;

            vel_site_data_p = velSiteDataPointer(iGlobLatDat, neigh_i, neigh_j,
                                                 neigh_k);

            if (vel_site_data_p == NULL || vel_site_data_p->proc_id == -1)
            {
              // it is a solid site and the velocity is
              // assumed to be zero
              v[i][j][k][0] = v[i][j][k][1] = v[i][j][k][2] = 0.0F;
              continue;
            }
            if (mNetworkTopology->LocalRank != vel_site_data_p->proc_id)
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
            else if (mNetworkTopology->LocalRank == vel_site_data_p->proc_id)
            {
              // the local counter is set equal to the global one
              // and the local velocity is calculated
              vel_site_data_p->counter = counter;
              double density;

              D3Q15::CalculateDensityAndVelocity(
                                                 &iLocalLatDat.FOld[vel_site_data_p->site_id
                                                     * c1Plusc2], density, vx,
                                                 vy, vz);

              v[i][j][k][0] = vel_site_data_p->vx = vx / density;
              v[i][j][k][1] = vel_site_data_p->vy = vy / density;
              v[i][j][k][2] = vel_site_data_p->vz = vz / density;
            }
            else
            {
              vel_site_data_p->counter = counter;

              int m =
                  from_proc_id_to_neigh_proc_index[vel_site_data_p->proc_id];

              mNeighProcs[m]->s_to_send[3 * mNeighProcs[m]->send_vs + 0]
                  = neigh_i;
              mNeighProcs[m]->s_to_send[3 * mNeighProcs[m]->send_vs + 1]
                  = neigh_j;
              mNeighProcs[m]->s_to_send[3 * mNeighProcs[m]->send_vs + 2]
                  = neigh_k;
              ++mNeighProcs[m]->send_vs;
            }
          }
        }
      }
    }

    // Constructor, populating fields from lattice data objects.
    StreaklineDrawer::StreaklineDrawer(const topology::NetworkTopology * iNetworkTopology,
                                       lb::LocalLatticeData &iLocalLatDat,
                                       lb::GlobalLatticeData &iGlobLatDat,
                                       bool & oSuccess)
    {
      oSuccess = true;
      mNetworkTopology = iNetworkTopology;

      int inlet_sites;
      const int *neigh_proc_id;

      lb::BlockData * lBlock;
      VelSiteData *vel_site_data_p;

      particleVec.reserve(10000);
      nParticles = 0;

      particleSeedVec.reserve(100);

      nParticleSeeds = 0;

      num_blocks = iGlobLatDat.BlockCount;
      velocity_field = new VelocityField[iGlobLatDat.BlockCount];

      for (int n = 0; n < iGlobLatDat.BlockCount; n++)
      {
        velocity_field[n].vel_site_data = NULL;
      }

      counter = 1;
      inlet_sites = 0;
      int n = -1;

      for (int i = 0; i < iGlobLatDat.SitesX; i += iGlobLatDat.BlockSize)
        for (int j = 0; j < iGlobLatDat.SitesY; j += iGlobLatDat.BlockSize)
          for (int k = 0; k < iGlobLatDat.SitesZ; k += iGlobLatDat.BlockSize)
          {

            lBlock = &iGlobLatDat.Blocks[++n];

            if (lBlock->site_data == NULL)
            {
              continue;
            }

            int m = -1;

            for (int site_i = i; site_i < i + iGlobLatDat.BlockSize; site_i++)
              for (int site_j = j; site_j < j + iGlobLatDat.BlockSize; site_j++)
                for (int site_k = k; site_k < k + iGlobLatDat.BlockSize; site_k++)
                {

                  m++;
                  if (mNetworkTopology->LocalRank
                      != lBlock->ProcessorRankForEachBlockSite[m])
                    continue;

                  for (int neigh_i = util::max(0, site_i - 1); neigh_i
                      <= util::min(iGlobLatDat.SitesX - 1, site_i + 1); neigh_i++)
                    for (int neigh_j = util::max(0, site_j - 1); neigh_j
                        <= util::min(iGlobLatDat.SitesY - 1, site_j + 1); neigh_j++)
                      for (int neigh_k = util::max(0, site_k - 1); neigh_k
                          <= util::min(iGlobLatDat.SitesZ - 1, site_k + 1); neigh_k++)
                      {

                        neigh_proc_id
                            = iGlobLatDat.GetProcIdFromGlobalCoords(neigh_i,
                                                                    neigh_j,
                                                                    neigh_k);

                        if (neigh_proc_id == NULL || *neigh_proc_id
                            == (1 << 30))
                        {
                          continue;
                        }

                        initializeVelFieldBlock(iGlobLatDat, neigh_i, neigh_j,
                                                neigh_k, *neigh_proc_id);

                        if (mNetworkTopology->LocalRank == *neigh_proc_id)
                          continue;

                        vel_site_data_p = velSiteDataPointer(iGlobLatDat,
                                                             neigh_i, neigh_j,
                                                             neigh_k);

                        if (vel_site_data_p->counter == counter)
                          continue;

                        vel_site_data_p->counter = counter;

                        bool seenSelf = false;
                        for (unsigned int mm = 0; mm < mNeighProcs.size() && !seenSelf; mm++)
                        {
                          if (*neigh_proc_id == mNeighProcs[mm]->id)
                          {
                            seenSelf = true;
                            ++mNeighProcs[mm]->send_vs;
                          }
                        }
                        if (seenSelf)
                          continue;

                        NeighProc * lNew = new NeighProc();

                        lNew->id = *neigh_proc_id;
                        lNew->send_vs = 1;
                        mNeighProcs.push_back(lNew);
                      }

                  int lSiteIndex = lBlock->site_data[m];

                  // if the lattice site is an not inlet one
                  if (iLocalLatDat.GetSiteType(lSiteIndex) != lb::INLET_TYPE)
                  {
                    continue;
                  }
                  ++inlet_sites;

                  if (inlet_sites % 50 != 0)
                    continue;

                  if (nParticleSeeds == particleSeedVec.capacity())
                  {
                    particleSeedVec.reserve(2 * particleSeedVec.capacity());
                  }
                  particleSeedVec[nParticleSeeds].x = (float) site_i;
                  particleSeedVec[nParticleSeeds].y = (float) site_j;
                  particleSeedVec[nParticleSeeds].z = (float) site_k;
                  particleSeedVec[nParticleSeeds].inlet_id
                      = iLocalLatDat.GetBoundaryId(lSiteIndex);
                  ++nParticleSeeds;
                }
          }
      shared_vs = 0;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        shared_vs += mNeighProcs[m]->send_vs;
      }
      if (shared_vs > 0)
      {
        s_to_send = new short int[3 * shared_vs];
        s_to_recv = new short int[3 * shared_vs];

        v_to_send = new float[3 * shared_vs];
        v_to_recv = new float[3 * shared_vs];
      }
      shared_vs = 0;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m]->s_to_send = &s_to_send[shared_vs * 3];
        mNeighProcs[m]->s_to_recv = &s_to_recv[shared_vs * 3];

        mNeighProcs[m]->v_to_send = &v_to_send[shared_vs * 3];
        mNeighProcs[m]->v_to_recv = &v_to_recv[shared_vs * 3];

        shared_vs += mNeighProcs[m]->send_vs;

        mNeighProcs[m]->send_vs = 0;
      }

      particles_to_send_max = 1000;
      particles_to_recv_max = 1000;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m]->p_to_send.reserve(5 * particles_to_send_max);
        mNeighProcs[m]->p_to_recv.reserve(5 * particles_to_recv_max);
      }

      req = new MPI_Request[2 * iNetworkTopology->ProcessorCount];

      from_proc_id_to_neigh_proc_index
          = new short int[iNetworkTopology->ProcessorCount];

      for (int m = 0; m < iNetworkTopology->ProcessorCount; m++)
      {
        from_proc_id_to_neigh_proc_index[m] = -1;
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        from_proc_id_to_neigh_proc_index[mNeighProcs[m]->id] = m;
      }

      counter = 0;

      for (n = 0; n < iGlobLatDat.BlockCount; n++)
      {
        if (velocity_field[n].vel_site_data == NULL)
          continue;

        for (int m = 0; m < iGlobLatDat.SitesPerBlockVolumeUnit; m++)
        {
          velocity_field[n].vel_site_data[m].counter = counter;
        }
        if (iGlobLatDat.Blocks[n].site_data == NULL)
          continue;

        for (int m = 0; m < iGlobLatDat.SitesPerBlockVolumeUnit; m++)
        {
          velocity_field[n].vel_site_data[m].site_id
              = iGlobLatDat.Blocks[n].site_data[m];
        }
      }
      procs = iNetworkTopology->ProcessorCount;
    }

    // Reset the streakline drawer.
    void StreaklineDrawer::Restart()
    {
      nParticles = 0;
    }

    // Communicate site ids to other processors.
    void StreaklineDrawer::communicateSiteIds()
    {
#ifndef NOMPI
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Irecv(&mNeighProcs[m]->recv_vs, 1, MPI_INT, mNeighProcs[m]->id, 30,
                  MPI_COMM_WORLD, &req[procs + mNeighProcs[m]->id]);
        MPI_Isend(&mNeighProcs[m]->send_vs, 1, MPI_INT, mNeighProcs[m]->id, 30,
                  MPI_COMM_WORLD, &req[mNeighProcs[m]->id]);
        MPI_Wait(&req[mNeighProcs[m]->id], status);
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Wait(&req[procs + mNeighProcs[m]->id], status);

        if (mNeighProcs[m]->recv_vs > 0)
        {
          MPI_Irecv(mNeighProcs[m]->s_to_recv, mNeighProcs[m]->recv_vs * 3,
                    MPI_SHORT, mNeighProcs[m]->id, 40, MPI_COMM_WORLD,
                    &req[procs + mNeighProcs[m]->id]);
        }
        if (mNeighProcs[m]->send_vs > 0)
        {
          MPI_Isend(mNeighProcs[m]->s_to_send, mNeighProcs[m]->send_vs * 3,
                    MPI_SHORT, mNeighProcs[m]->id, 40, MPI_COMM_WORLD,
                    &req[mNeighProcs[m]->id]);
          MPI_Wait(&req[mNeighProcs[m]->id], status);
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m]->recv_vs > 0)
        {
          MPI_Wait(&req[procs + mNeighProcs[m]->id], status);
        }
      }
#endif // NOMPI
    }

    // Communicate velocities to other processors.
    void StreaklineDrawer::communicateVelocities(lb::GlobalLatticeData &iGlobLatDat)
    {
#ifndef NOMPI
      int site_i, site_j, site_k;
      int neigh_i, neigh_j, neigh_k;

      VelSiteData *vel_site_data_p;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m]->send_vs > 0)
        {
          MPI_Irecv(mNeighProcs[m]->v_to_recv, mNeighProcs[m]->send_vs * 3,
                    MPI_FLOAT, mNeighProcs[m]->id, 30, MPI_COMM_WORLD,
                    &req[procs + mNeighProcs[m]->id]);
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        for (int n = 0; n < mNeighProcs[m]->recv_vs; n++)
        {
          site_i = mNeighProcs[m]->s_to_recv[3 * n + 0];
          site_j = mNeighProcs[m]->s_to_recv[3 * n + 1];
          site_k = mNeighProcs[m]->s_to_recv[3 * n + 2];

          vel_site_data_p = velSiteDataPointer(iGlobLatDat, site_i, site_j,
                                               site_k);

          if (vel_site_data_p != NULL)
          {
            mNeighProcs[m]->v_to_send[3 * n + 0] = vel_site_data_p->vx;
            mNeighProcs[m]->v_to_send[3 * n + 1] = vel_site_data_p->vy;
            mNeighProcs[m]->v_to_send[3 * n + 2] = vel_site_data_p->vz;
          }
          else
          {
            mNeighProcs[m]->v_to_send[3 * n + 0] = 0.;
            mNeighProcs[m]->v_to_send[3 * n + 1] = 0.;
            mNeighProcs[m]->v_to_send[3 * n + 2] = 0.;
          }
        }
        if (mNeighProcs[m]->recv_vs > 0)
        {
          MPI_Isend(mNeighProcs[m]->v_to_send, mNeighProcs[m]->recv_vs * 3,
                    MPI_FLOAT, mNeighProcs[m]->id, 30, MPI_COMM_WORLD,
                    &req[mNeighProcs[m]->id]);
          MPI_Wait(&req[mNeighProcs[m]->id], status);
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m]->send_vs <= 0)
          continue;

        MPI_Wait(&req[procs + mNeighProcs[m]->id], status);

        for (int n = 0; n < mNeighProcs[m]->send_vs; n++)
        {
          neigh_i = mNeighProcs[m]->s_to_send[3 * n + 0];
          neigh_j = mNeighProcs[m]->s_to_send[3 * n + 1];
          neigh_k = mNeighProcs[m]->s_to_send[3 * n + 2];

          vel_site_data_p = velSiteDataPointer(iGlobLatDat, neigh_i, neigh_j,
                                               neigh_k);

          if (vel_site_data_p != NULL)
          {
            vel_site_data_p->vx = mNeighProcs[m]->v_to_recv[3 * n + 0];
            vel_site_data_p->vy = mNeighProcs[m]->v_to_recv[3 * n + 1];
            vel_site_data_p->vz = mNeighProcs[m]->v_to_recv[3 * n + 2];
          }
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m]->send_vs = 0;
      }
#endif // NOMPI
    }

    // Update the velocity field.
    void StreaklineDrawer::updateVelField(int stage_id,
                                          lb::GlobalLatticeData &iGlobLatDat,
                                          lb::LocalLatticeData &iLocalLatDat)
    {
      float v[2][2][2][3], interp_v[3];
      float vel;

      int particles_temp;
      int is_interior;
      int n;

      particles_temp = nParticles;

      for (n = particles_temp - 1; n >= 0; n--)
      {
        localVelField(n, v, &is_interior, iGlobLatDat, iLocalLatDat);

        if (stage_id == 0 && !is_interior)
          continue;

        particleVelocity(&particleVec[n], v, interp_v);
        vel = interp_v[0] * interp_v[0] + interp_v[1] * interp_v[1]
            + interp_v[2] * interp_v[2];

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
    void StreaklineDrawer::communicateParticles(lb::GlobalLatticeData &iGlobLatDat)
    {
#ifndef NOMPI
      int site_i, site_j, site_k;
      int particles_temp;
      MPI_Status status;

      VelSiteData *vel_site_data_p;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Irecv(&mNeighProcs[m]->recv_ps, 1, MPI_INT, mNeighProcs[m]->id, 30,
                  MPI_COMM_WORLD, &req[procs + mNeighProcs[m]->id]);
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m]->send_ps = 0;
      }

      particles_temp = nParticles;

      for (int n = particles_temp - 1; n >= 0; n--)
      {
        site_i = (int) particleVec[n].x;
        site_j = (int) particleVec[n].y;
        site_k = (int) particleVec[n].z;

        vel_site_data_p = velSiteDataPointer(iGlobLatDat, site_i, site_j,
                                             site_k);

        if (vel_site_data_p == NULL || mNetworkTopology->LocalRank
            == vel_site_data_p->proc_id || vel_site_data_p->proc_id == -1)
        {
          continue;
        }
        int m = from_proc_id_to_neigh_proc_index[vel_site_data_p->proc_id];

        if (mNeighProcs[m]->send_ps == particles_to_send_max)
        {
          particles_to_send_max *= 2;
          mNeighProcs[m]->p_to_send.reserve(5 * particles_to_send_max);

        }

        mNeighProcs[m]->p_to_send[5 * mNeighProcs[m]->send_ps + 0]
            = particleVec[n].x;
        mNeighProcs[m]->p_to_send[5 * mNeighProcs[m]->send_ps + 1]
            = particleVec[n].y;
        mNeighProcs[m]->p_to_send[5 * mNeighProcs[m]->send_ps + 2]
            = particleVec[n].z;
        mNeighProcs[m]->p_to_send[5 * mNeighProcs[m]->send_ps + 3]
            = particleVec[n].vel;
        mNeighProcs[m]->p_to_send[5 * mNeighProcs[m]->send_ps + 4]
            = particleVec[n].inlet_id + 0.1;
        ++mNeighProcs[m]->send_ps;

        deleteParticle(n);
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Isend(&mNeighProcs[m]->send_ps, 1, MPI_INT, mNeighProcs[m]->id, 30,
                  MPI_COMM_WORLD, &req[mNeighProcs[m]->id]);
        MPI_Wait(&req[mNeighProcs[m]->id], &status);
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        MPI_Wait(&req[procs + mNeighProcs[m]->id], &status);
      }

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m]->send_ps > 0)
        {
          MPI_Isend(&mNeighProcs[m]->p_to_send[0], mNeighProcs[m]->send_ps * 5,
                    MPI_FLOAT, mNeighProcs[m]->id, 40, MPI_COMM_WORLD,
                    &req[mNeighProcs[m]->id]);
          MPI_Wait(&req[mNeighProcs[m]->id], &status);
        }
      }
      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        if (mNeighProcs[m]->recv_ps > 0)
        {
          if (mNeighProcs[m]->recv_ps > particles_to_recv_max)
          {
            particles_to_recv_max *= 2;
            particles_to_recv_max = util::max(particles_to_recv_max,
                                              mNeighProcs[m]->recv_ps);
            mNeighProcs[m]->p_to_recv.reserve(5 * particles_to_recv_max);
          }
          MPI_Irecv(&mNeighProcs[m]->p_to_recv[0], mNeighProcs[m]->recv_ps * 5,
                    MPI_FLOAT, mNeighProcs[m]->id, 40, MPI_COMM_WORLD,
                    &req[procs + mNeighProcs[m]->id]);
          MPI_Wait(&req[procs + mNeighProcs[m]->id], &status);

          for (int n = 0; n < mNeighProcs[m]->recv_ps; n++)
          {
            createParticle(mNeighProcs[m]->p_to_recv[5 * n + 0],
                           mNeighProcs[m]->p_to_recv[5 * n + 1],
                           mNeighProcs[m]->p_to_recv[5 * n + 2],
                           mNeighProcs[m]->p_to_recv[5 * n + 3],
                           (int) mNeighProcs[m]->p_to_recv[5 * n + 4]);
          }
        }
      }
#endif // NOMPI
    }

    // Render the streaklines
    void StreaklineDrawer::render(lb::GlobalLatticeData &iGlobLatDat)
    {
      float screen_max[2];
      float scale[2];
      float p1[3], p2[3];

      ColPixel col_pixel;

      int pixels_x = vis::controller->mScreen.PixelsX;
      int pixels_y = vis::controller->mScreen.PixelsY;

      screen_max[0] = vis::controller->mScreen.MaxXValue;
      screen_max[1] = vis::controller->mScreen.MaxYValue;

      scale[0] = vis::controller->mScreen.ScaleX;
      scale[1] = vis::controller->mScreen.ScaleY;

      for (unsigned int n = 0; n < nParticles; n++)
      {
        p1[0] = particleVec[n].x - float(iGlobLatDat.SitesX >> 1);
        p1[1] = particleVec[n].y - float(iGlobLatDat.SitesY >> 1);
        p1[2] = particleVec[n].z - float(iGlobLatDat.SitesZ >> 1);

        vis::controller->project(p1, p2);

        p2[0] = int(scale[0] * (p2[0] + screen_max[0]));
        p2[1] = int(scale[1] * (p2[1] + screen_max[1]));

        int i = int(p2[0]);
        int j = int(p2[1]);

        if (! (i < 0 || i >= pixels_x || j < 0 || j >= pixels_y))
        {
          col_pixel.particle_vel = particleVec[n].vel;
          col_pixel.particle_z = p2[2];
          col_pixel.particle_inlet_id = particleVec[n].inlet_id;
          col_pixel.i = PixelId(i, j);
          col_pixel.i.isStreakline = true;

          vis::controller->writePixel(&col_pixel);
        }
      }
    }

    // Draw streaklines
    void StreaklineDrawer::StreakLines(int time_steps,
                                       int time_steps_per_cycle,
                                       lb::GlobalLatticeData &iGlobLatDat,
                                       lb::LocalLatticeData &iLocalLatDat)
    {
      int particle_creation_period = util::max(1, (int) (time_steps_per_cycle
          / 5000.0F));

      if (time_steps % (int) (time_steps_per_cycle
          / vis::controller->streaklines_per_pulsatile_period)
          <= (vis::controller->streakline_length / 100.0F)
              * (time_steps_per_cycle
                  / vis::controller->streaklines_per_pulsatile_period)
          && time_steps % particle_creation_period == 0)
      {
        createSeedParticles();
      }

      ++counter;

      updateVelField(0, iGlobLatDat, iLocalLatDat);
      communicateSiteIds();
      communicateVelocities(iGlobLatDat);
      updateVelField(1, iGlobLatDat, iLocalLatDat);
      updateParticles();
      communicateParticles(iGlobLatDat);
    }

    // Destructor
    StreaklineDrawer::~StreaklineDrawer()
    {
      delete from_proc_id_to_neigh_proc_index;
      delete req;

      for (unsigned int m = 0; m < mNeighProcs.size(); m++)
      {
        mNeighProcs[m]->p_to_recv.clear();
        mNeighProcs[m]->p_to_send.clear();
      }

      if (shared_vs > 0)
      {
        delete v_to_recv;
        delete v_to_send;

        delete s_to_recv;
        delete s_to_send;
      }

      for (int m = 0; m < num_blocks; m++)
      {
        if (velocity_field[m].vel_site_data != NULL)
        {
          delete velocity_field[m].vel_site_data;
        }
      }

      delete velocity_field;
      particleSeedVec.clear();
      particleVec.clear();
    }

  }
}
