#include <math.h>
#include <stdlib.h>
#include <vector>
#include <limits>

#include "debug/Debugger.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "vis/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    // TODO RENAME THIS FUNCTION
    void RayTracer::AABBvsRay(const AABB* aabb,
                              const float inverseDirection[3],
                              const bool xyzComponentIsPositive[3],
                              float* t_near,
                              float* t_far)
    {
      float tx0 = (xyzComponentIsPositive[0]
        ? aabb->acc_2
        : aabb->acc_1) * inverseDirection[0];
      float tx1 = (xyzComponentIsPositive[0]
        ? aabb->acc_1
        : aabb->acc_2) * inverseDirection[0];
      float ty0 = (xyzComponentIsPositive[1]
        ? aabb->acc_4
        : aabb->acc_3) * inverseDirection[1];
      float ty1 = (xyzComponentIsPositive[1]
        ? aabb->acc_3
        : aabb->acc_4) * inverseDirection[1];
      float tz0 = (xyzComponentIsPositive[2]
        ? aabb->acc_6
        : aabb->acc_5) * inverseDirection[2];
      float tz1 = (xyzComponentIsPositive[2]
        ? aabb->acc_5
        : aabb->acc_6) * inverseDirection[2];

      *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
      *t_far = fminf(tx1, fminf(ty1, tz1));
    }

    /**
     * Update a colour vector to include a section of known length through a
     * solid of known colour.
     *
     * @param dt
     * @param palette
     * @param col
     */
    void RayTracer::UpdateColour(float dt, const float palette[3], float col[3])
    {
      for (int ii = 0; ii < 3; ++ii)
      {
        col[ii] += dt * palette[ii];
      }
    }

    void RayTracer::UpdateRayData(const float flow_field[3],
                                  float ray_t,
                                  float ray_segment,
                                  Ray* bCurrentRay)
    {
      if (flow_field[0] < 0.0F)
      {
        return; // solid voxel
      }

      float palette[3];

      // update the volume rendering of the velocity flow field
      ColPixel::PickColour(flow_field[1] * mDomainStats->velocity_threshold_max_inv, palette);

      UpdateColour(ray_segment, palette, bCurrentRay->VelocityColour);

      if (mVisSettings->mStressType != lb::ShearStress)
      {
        // update the volume rendering of the von Mises stress flow field
        float scaled_stress = flow_field[2] * mDomainStats->stress_threshold_max_inv;

        ColPixel::PickColour(scaled_stress, palette);

        UpdateColour(ray_segment, palette, bCurrentRay->StressColour);
      }

      bCurrentRay->Length += ray_segment;

      if (bCurrentRay->Density >= 0.0F)
      {
        return;
      }

      bCurrentRay->MinT = ray_t;

      // keep track of the density nearest to the view point
      bCurrentRay->Density = flow_field[0];

      // keep track of the stress nearest to the view point
      bCurrentRay->Stress = flow_field[2];
    }

    void RayTracer::TraverseVoxels(const float block_min[3],
                                   const float block_x[3],
                                   const float voxel_flow_field[],
                                   float t,
                                   Ray* bCurrentRay,
                                   const bool xyz_is_1[3])
    {
      unsigned int i_vec[3];

      for (int i = 0; i < 3; i++)
      {
        i_vec[i] = util::NumericalFunctions::enforceBounds<unsigned int>(block_x[i],
                                                                         0,
                                                                         block_size_1);
      }

      float t_max[3];
      for (int i = 0; i < 3; i++)
      {
        t_max[i] = (block_min[i] + (float) (xyz_is_1[i]
          ? i_vec[i] + 1
          : i_vec[i])) * bCurrentRay->InverseDirection[i];
      }

      unsigned int i = i_vec[0] * block_size2;
      unsigned int j = i_vec[1] * mLatDat->GetBlockSize();
      unsigned int k = i_vec[2];

      while (true)
      {
        if (t_max[0] < t_max[1])
        {
          if (t_max[0] < t_max[2])
          {
            UpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS],
                          t,
                          t_max[0] - t,
                          bCurrentRay);

            if (xyz_is_1[0])
            {
              if ( (i += block_size2) >= block_size3)
              {
                return;
              }
              t = t_max[0];
              t_max[0] += bCurrentRay->InverseDirection[0];
            }
            else
            {
              if (i < block_size2)
              {
                return;
              }
              else
              {
                i -= block_size2;
              }
              t = t_max[0];
              t_max[0] -= bCurrentRay->InverseDirection[0];
            }
          }
          else
          {
            UpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS],
                          t,
                          t_max[2] - t,
                          bCurrentRay);

            if (xyz_is_1[2])
            {
              if (++k >= mLatDat->GetBlockSize())
              {
                return;
              }
              t = t_max[2];
              t_max[2] += bCurrentRay->InverseDirection[2];
            }
            else
            {
              if (k == 0)
              {
                return;
              }
              else
              {
                --k;
              }
              t = t_max[2];
              t_max[2] -= bCurrentRay->InverseDirection[2];
            }
          }
        }
        else
        {
          if (t_max[1] < t_max[2])
          {
            UpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS],
                          t,
                          t_max[1] - t,
                          bCurrentRay);

            if (xyz_is_1[1])
            {
              if ( (j += mLatDat->GetBlockSize()) >= block_size2)
              {
                return;
              }
              t = t_max[1];
              t_max[1] += bCurrentRay->InverseDirection[1];
            }
            else
            {
              if (j < mLatDat->GetBlockSize())
              {
                return;
              }
              else
              {
                j -= mLatDat->GetBlockSize();
              }
              t = t_max[1];
              t_max[1] -= bCurrentRay->InverseDirection[1];
            }
          }
          else
          {
            UpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS],
                          t,
                          t_max[2] - t,
                          bCurrentRay);

            if (xyz_is_1[2])
            {
              if (++k >= mLatDat->GetBlockSize())
              {
                return;
              }
              t = t_max[2];
              t_max[2] += bCurrentRay->InverseDirection[2];
            }
            else
            {
              if (k == 0)
              {
                return;
              }
              else
              {
                --k;
              }
              t = t_max[2];
              t_max[2] -= bCurrentRay->InverseDirection[2];
            }
          }
        }
      }
    }

    void RayTracer::TraverseBlocks(const Cluster* cluster,
                                   const bool xyz_Is_1[3],
                                   const float ray_dx[3],
                                   float **block_flow_field,
                                   Ray *bCurrentRay)
    {

      int cluster_blocks_vec[3];
      cluster_blocks_vec[0] = cluster->blocks_x - 1;
      cluster_blocks_vec[1] = cluster->blocks_y - 1;
      cluster_blocks_vec[2] = cluster->blocks_z - 1;
      int cluster_blocks_z = cluster->blocks_z;
      int cluster_blocks_yz = (int) cluster->blocks_y * (int) cluster->blocks_z;
      int cluster_blocks = (int) cluster->blocks_x * cluster_blocks_yz;

      int i_vec[3];
      float block_min[3];

      for (int l = 0; l < 3; l++)
      {
        i_vec[l] = util::NumericalFunctions::enforceBounds(cluster_blocks_vec[l],
                                                           0,
                                                           (int) (mBlockSizeInverse * ray_dx[l]));
        block_min[l] = (float) i_vec[l] * mBlockSizeFloat - ray_dx[l];
      }

      int i = i_vec[0] * cluster_blocks_yz;
      int j = i_vec[1] * cluster_blocks_z;
      int k = i_vec[2];

      float block_x[3];
      if (block_flow_field[i + j + k] != NULL)
      {
        block_x[0] = -block_min[0];
        block_x[1] = -block_min[1];
        block_x[2] = -block_min[2];

        TraverseVoxels(block_min, block_x, block_flow_field[i + j + k], 0.0F, bCurrentRay, xyz_Is_1);
      }

      float t_max[3];
      float t_delta[3];
      for (int l = 0; l < 3; l++)
      {
        t_max[l] = (xyz_Is_1[l]
          ? block_min[l] + mBlockSizeFloat
          : block_min[l]) * bCurrentRay->InverseDirection[l];
        t_delta[l] = mBlockSizeFloat * bCurrentRay->InverseDirection[l];
      }

      while (true)
      {
        if (t_max[0] < t_max[1])
        {
          if (t_max[0] < t_max[2])
          {
            if (xyz_Is_1[0])
            {
              if ( (i += cluster_blocks_yz) >= cluster_blocks)
                return;
              block_min[0] += mBlockSizeFloat;
            }
            else
            {
              if ( (i -= cluster_blocks_yz) < 0)
                return;
              block_min[0] -= mBlockSizeFloat;
            }

            if (block_flow_field[i + j + k] != NULL)
            {
              block_x[0] = t_max[0] * bCurrentRay->Direction[0] - block_min[0];
              block_x[1] = t_max[0] * bCurrentRay->Direction[1] - block_min[1];
              block_x[2] = t_max[0] * bCurrentRay->Direction[2] - block_min[2];

              TraverseVoxels(block_min,
                             block_x,
                             block_flow_field[i + j + k],
                             t_max[0],
                             bCurrentRay,
                             xyz_Is_1);
            }

            t_max[0] = xyz_Is_1[0]
              ? t_max[0] + t_delta[0]
              : t_max[0] - t_delta[0];
          }
          else
          {
            if (xyz_Is_1[2])
            {
              if (++k >= cluster_blocks_z)
                return;
              block_min[2] += mBlockSizeFloat;
            }
            else
            {
              if (--k < 0)
                return;
              block_min[2] -= mBlockSizeFloat;
            }

            if (block_flow_field[i + j + k] != NULL)
            {
              block_x[0] = t_max[2] * bCurrentRay->Direction[0] - block_min[0];
              block_x[1] = t_max[2] * bCurrentRay->Direction[1] - block_min[1];
              block_x[2] = t_max[2] * bCurrentRay->Direction[2] - block_min[2];

              TraverseVoxels(block_min,
                             block_x,
                             block_flow_field[i + j + k],
                             t_max[2],
                             bCurrentRay,
                             xyz_Is_1);
            }

            t_max[2] = xyz_Is_1[2]
              ? t_max[2] + t_delta[2]
              : t_max[2] - t_delta[2];
          }
        }
        else
        {
          if (t_max[1] < t_max[2])
          {
            if (xyz_Is_1[1])
            {
              if ( (j += cluster_blocks_z) >= cluster_blocks_yz)
                return;
              block_min[1] += mBlockSizeFloat;
            }
            else
            {
              if ( (j -= cluster_blocks_z) < 0)
                return;
              block_min[1] -= mBlockSizeFloat;
            }

            if (block_flow_field[i + j + k] != NULL)
            {
              block_x[0] = t_max[1] * bCurrentRay->Direction[0] - block_min[0];
              block_x[1] = t_max[1] * bCurrentRay->Direction[1] - block_min[1];
              block_x[2] = t_max[1] * bCurrentRay->Direction[2] - block_min[2];

              TraverseVoxels(block_min,
                             block_x,
                             block_flow_field[i + j + k],
                             t_max[1],
                             bCurrentRay,
                             xyz_Is_1);
            }

            t_max[1] = xyz_Is_1[1]
              ? t_max[1] + t_delta[1]
              : t_max[1] - t_delta[1];
          }
          else
          {
            if (xyz_Is_1[2])
            {
              if (++k >= cluster_blocks_z)
                return;
              block_min[2] += mBlockSizeFloat;
            }
            else
            {
              if (--k < 0)
                return;
              block_min[2] -= mBlockSizeFloat;
            }

            if (block_flow_field[i + j + k] != NULL)
            {
              block_x[0] = t_max[2] * bCurrentRay->Direction[0] - block_min[0];
              block_x[1] = t_max[2] * bCurrentRay->Direction[1] - block_min[1];
              block_x[2] = t_max[2] * bCurrentRay->Direction[2] - block_min[2];

              TraverseVoxels(block_min,
                             block_x,
                             block_flow_field[i + j + k],
                             t_max[2],
                             bCurrentRay,
                             xyz_Is_1);
            }

            t_max[2] = xyz_Is_1[2]
              ? t_max[2] + t_delta[2]
              : t_max[2] - t_delta[2];
          }
        }
      }
    }

    void RayTracer::BuildClusters()
    {
      std::vector<MinLocation> clusterBlockMins;

      const int n_x[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, +0, +0, +0, +0, +0, +0, +0, +0, +1,
                          +1, +1, +1, +1, +1, +1, +1, +1 };
      const int n_y[] = { -1, -1, -1, +0, +0, +0, +1, +1, +1, -1, -1, -1, +0, +0, +1, +1, +1, -1,
                          -1, -1, +0, +0, +0, +1, +1, +1 };
      const int n_z[] = { -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +1, -1, +0, +1, -1,
                          +0, +1, -1, +0, +1, -1, +0, +1 };

      short int *cluster_id = new short int[mLatDat->GetBlockCount()];
      for (unsigned int n = 0; n < mLatDat->GetBlockCount(); n++)
      {
        cluster_id[n] = -1;
      }

      int blocks_buffer_size = 10000;

      std::vector<BlockLocation> *block_location_a =
          new std::vector<BlockLocation>(blocks_buffer_size);
      std::vector<BlockLocation> *block_location_b =
          new std::vector<BlockLocation>(blocks_buffer_size);

      bool *is_block_visited = new bool[mLatDat->GetBlockCount()];
      for (unsigned int n = 0; n < mLatDat->GetBlockCount(); n++)
      {
        is_block_visited[n] = false;
      }

      int n = -1;

      for (unsigned int i = 0; i < mLatDat->GetXBlockCount(); i++)
      {
        for (unsigned int j = 0; j < mLatDat->GetYBlockCount(); j++)
        {
          for (unsigned int k = 0; k < mLatDat->GetZBlockCount(); k++)
          {
            n++;

            geometry::LatticeData::BlockData * lBlock = mLatDat->GetBlock((unsigned int) n);
            if (lBlock->ProcessorRankForEachBlockSite == NULL || is_block_visited[n])
            {
              continue;
            }

            is_block_visited[n] = 1;

            int blocks_a = 0;

            for (unsigned int m = 0; m < mLatDat->GetSitesPerBlockVolumeUnit(); m++)
            {
              if ((int) mNetworkTopology->GetLocalRank()
                  == lBlock->ProcessorRankForEachBlockSite[m])
              {
                BlockLocation& tempBlockLoc = block_location_a->at(0);
                tempBlockLoc.i = i;
                tempBlockLoc.j = j;
                tempBlockLoc.k = k;
                blocks_a = 1;
                break;
              }
            }

            if (blocks_a == 0)
            {
              continue;
            }

            Cluster lNewCluster;
            MinLocation blockMinimum;

            blockMinimum.i = i;
            blockMinimum.j = j;
            blockMinimum.k = k;

            unsigned short int cluster_block_max_i = i;
            unsigned short int cluster_block_max_j = j;
            unsigned short int cluster_block_max_k = k;

            cluster_id[n] = mClusters.size();
            bool are_blocks_incrementing = true;

            while (are_blocks_incrementing)
            {
              int blocks_b = 0;
              are_blocks_incrementing = false;

              for (int index_a = 0; index_a < blocks_a; index_a++)
              {
                const BlockLocation& tempBlockLoc = block_location_a->at(index_a);

                for (int l = 0; l < 26; l++)
                {
                  int neigh_i = tempBlockLoc.i + n_x[l];
                  int neigh_j = tempBlockLoc.j + n_y[l];
                  int neigh_k = tempBlockLoc.k + n_z[l];

                  if (neigh_i == -1 || neigh_i == (int) mLatDat->GetXBlockCount())
                    continue;
                  if (neigh_j == -1 || neigh_j == (int) mLatDat->GetYBlockCount())
                    continue;
                  if (neigh_k == -1 || neigh_k == (int) mLatDat->GetZBlockCount())
                    continue;

                  int block_id = (neigh_i * mLatDat->GetYBlockCount() + neigh_j)
                      * mLatDat->GetZBlockCount() + neigh_k;

                  if (is_block_visited[block_id]
                      || (lBlock = mLatDat->GetBlock(block_id))->ProcessorRankForEachBlockSite
                          == NULL)
                  {
                    continue;
                  }

                  bool is_site_found = false;
                  for (unsigned int m = 0; m < mLatDat->GetSitesPerBlockVolumeUnit(); m++)
                  {
                    if ((int) mNetworkTopology->GetLocalRank()
                        == lBlock->ProcessorRankForEachBlockSite[m])
                    {
                      is_site_found = true;
                      break;
                    }
                  }

                  if (!is_site_found)
                  {
                    continue;
                  }

                  is_block_visited[block_id] = 1;

                  if (blocks_b == blocks_buffer_size)
                  {
                    blocks_buffer_size *= 2;
                    block_location_a->resize(blocks_buffer_size);
                    block_location_b->resize(blocks_buffer_size);
                  }
                  are_blocks_incrementing = true;

                  BlockLocation& tempBlockLoc = block_location_b->at(blocks_b);
                  tempBlockLoc.i = neigh_i;
                  tempBlockLoc.j = neigh_j;
                  tempBlockLoc.k = neigh_k;
                  ++blocks_b;

                  blockMinimum.i = util::NumericalFunctions::min((int) blockMinimum.i, neigh_i);
                  blockMinimum.j = util::NumericalFunctions::min((int) blockMinimum.j, neigh_j);
                  blockMinimum.k = util::NumericalFunctions::min((int) blockMinimum.k, neigh_k);

                  cluster_block_max_i = util::NumericalFunctions::max((int) cluster_block_max_i,
                                                                      neigh_i);
                  cluster_block_max_j = util::NumericalFunctions::max((int) cluster_block_max_j,
                                                                      neigh_j);
                  cluster_block_max_k = util::NumericalFunctions::max((int) cluster_block_max_k,
                                                                      neigh_k);

                  cluster_id[block_id] = mClusters.size();
                }
              }

              // swap pointers in block_location_a/_b
              std::vector<BlockLocation>* tempBlockLocation = block_location_a;
              block_location_a = block_location_b;
              block_location_b = tempBlockLocation;

              blocks_a = blocks_b;
            }

            lNewCluster.x[0] = blockMinimum.i * mLatDat->GetBlockSize() - 0.5F
                * mLatDat->GetXSiteCount();
            lNewCluster.x[1] = blockMinimum.j * mLatDat->GetBlockSize() - 0.5F
                * mLatDat->GetYSiteCount();
            lNewCluster.x[2] = blockMinimum.k * mLatDat->GetBlockSize() - 0.5F
                * mLatDat->GetZSiteCount();

            lNewCluster.blocks_x = 1 + cluster_block_max_i - blockMinimum.i;
            lNewCluster.blocks_y = 1 + cluster_block_max_j - blockMinimum.j;
            lNewCluster.blocks_z = 1 + cluster_block_max_k - blockMinimum.k;

            mClusters.push_back(lNewCluster);
            clusterBlockMins.push_back(blockMinimum);
          }
        }
      }

      delete block_location_b;
      delete block_location_a;

      delete[] is_block_visited;

      cluster_voxel = new float *[mLatDat->GetLocalFluidSiteCount() * VIS_FIELDS];

      cluster_flow_field = new float **[mClusters.size()];

      for (unsigned int lThisClusterId = 0; lThisClusterId < mClusters.size(); lThisClusterId++)
      {
        Cluster* cluster_p = &mClusters[lThisClusterId];

        cluster_flow_field[lThisClusterId] = new float *[cluster_p->blocks_x * cluster_p->blocks_y
            * cluster_p->blocks_z];

        int voxel_min[3], voxel_max[3];
        for (int l = 0; l < 3; l++)
        {
          voxel_min[l] = std::numeric_limits<int>::max();
          voxel_max[l] = std::numeric_limits<int>::min();
        }

        n = -1;

        for (int i = 0; i < cluster_p->blocks_x; i++)
        {
          for (int j = 0; j < cluster_p->blocks_y; j++)
          {
            for (int k = 0; k < cluster_p->blocks_z; k++)
            {
              ++n;

              int block_coord[3];
              MinLocation* mins = &clusterBlockMins[lThisClusterId];
              block_coord[0] = (i + mins->i) * mLatDat->GetBlockSize();
              block_coord[1] = (j + mins->j) * mLatDat->GetBlockSize();
              block_coord[2] = (k + mins->k) * mLatDat->GetBlockSize();

              int block_id = ( (i + mins->i) * mLatDat->GetYBlockCount() + (j + mins->j))
                  * mLatDat->GetZBlockCount() + (k + mins->k);

              cluster_flow_field[lThisClusterId][n] = NULL;

              if (cluster_id[block_id] != (short int) lThisClusterId)
              {
                continue;
              }

              geometry::LatticeData::BlockData * lBlock = mLatDat->GetBlock(block_id);

              cluster_flow_field[lThisClusterId][n]
                  = new float[mLatDat->GetSitesPerBlockVolumeUnit() * VIS_FIELDS];

              int m = -1;

              float siteLocOnBlock[3];
              for (siteLocOnBlock[0] = 0; siteLocOnBlock[0] < mLatDat->GetBlockSize(); siteLocOnBlock[0]++)
                for (siteLocOnBlock[1] = 0; siteLocOnBlock[1] < mLatDat->GetBlockSize(); siteLocOnBlock[1]++)
                  for (siteLocOnBlock[2] = 0; siteLocOnBlock[2] < mLatDat->GetBlockSize(); siteLocOnBlock[2]++)
                  {
                    unsigned int my_site_id;
                    my_site_id = lBlock->site_data[++m];

                    if (my_site_id & BIG_NUMBER3)
                    {
                      for (int l = 0; l < VIS_FIELDS; l++)
                        cluster_flow_field[lThisClusterId][n][m * VIS_FIELDS + l] = -1.0F;

                      continue;
                    }

                    for (int l = 0; l < VIS_FIELDS; l++)
                    {
                      cluster_flow_field[lThisClusterId][n][m * VIS_FIELDS + l] = 1.0F;
                    }

                    for (int l = 0; l < VIS_FIELDS; l++)
                    {
                      cluster_voxel[my_site_id * VIS_FIELDS + l]
                          = &cluster_flow_field[lThisClusterId][n][m * VIS_FIELDS + l];
                    }

                    for (int l = 0; l < 3; l++)
                    {
                      voxel_min[l] = util::NumericalFunctions::min<int>(voxel_min[l],
                                                                        siteLocOnBlock[l]
                                                                            + block_coord[l]);
                      voxel_max[l] = util::NumericalFunctions::max<int>(voxel_max[l],
                                                                        siteLocOnBlock[l]
                                                                            + block_coord[l]);
                    }

                  } // for ii[0..2]

            } // for k
          } // for j
        } // for i

        cluster_p->minmax_x[0] = (float) voxel_min[0] - 0.5F * (float) mLatDat->GetXSiteCount();
        cluster_p->minmax_y[0] = (float) voxel_min[1] - 0.5F * (float) mLatDat->GetYSiteCount();
        cluster_p->minmax_z[0] = (float) voxel_min[2] - 0.5F * (float) mLatDat->GetZSiteCount();

        cluster_p->minmax_x[1] = (float) (voxel_max[0] + 1) - 0.5F
            * (float) mLatDat->GetXSiteCount();
        cluster_p->minmax_y[1] = (float) (voxel_max[1] + 1) - 0.5F
            * (float) mLatDat->GetYSiteCount();
        cluster_p ->minmax_z[1] = (float) (voxel_max[2] + 1) - 0.5F
            * (float) mLatDat->GetZSiteCount();
      }
      delete[] cluster_id;
    }

    RayTracer::RayTracer(const topology::NetworkTopology * iNetworkTopology,
                         const geometry::LatticeData* iLatDat,
                         const DomainStats* iDomainStats,
                         Screen* iScreen,
                         Viewpoint* iViewpoint,
                         VisSettings* iVisSettings) :
      mNetworkTopology(iNetworkTopology), mLatDat(iLatDat), mDomainStats(iDomainStats),
          mScreen(iScreen), mViewpoint(iViewpoint), mVisSettings(iVisSettings),
          mBlockSizeFloat((float) mLatDat->GetBlockSize()),
          mBlockSizeInverse(1.F / mBlockSizeFloat), block_size2(mLatDat->GetBlockSize()
              * mLatDat->GetBlockSize()), block_size3(mLatDat->GetBlockSize() * block_size2),
          block_size_1(mLatDat->GetBlockSize() - 1), blocks_yz(mLatDat->GetYBlockCount()
              * mLatDat->GetZBlockCount())
    {
      BuildClusters();
    }

    void RayTracer::Render()
    {
      const float* const projectedUnitX = mScreen->GetUnitVectorProjectionX();
      const float* const projectedUnitY = mScreen->GetUnitVectorProjectionY();

      const float* viewpointCentre = mViewpoint->GetViewpointCentre();

      for (unsigned int clusterId = 0; clusterId < mClusters.size(); clusterId++)
      {
        const Cluster* thisCluster = &mClusters[clusterId];

        // the image-based projection of the mClusters bounding box is
        // calculated here
        float cluster_x[3];
        for (int l = 0; l < 3; l++)
        {
          cluster_x[l] = thisCluster->x[l] - viewpointCentre[l];
        }

        float **block_flow_field = cluster_flow_field[clusterId];

        float subimageMins[2], subimageMaxes[2];

        subimageMins[0] = subimageMins[1] = std::numeric_limits<float>::max();
        subimageMaxes[0] = subimageMaxes[1] = std::numeric_limits<float>::min();

        float p1[3];
        for (int i = 0; i < 2; i++)
        {
          p1[0] = thisCluster->minmax_x[i];

          for (int j = 0; j < 2; j++)
          {
            p1[1] = thisCluster->minmax_y[j];

            for (int k = 0; k < 2; k++)
            {
              p1[2] = thisCluster->minmax_z[k];

              float p2[3];
              mViewpoint->Project(p1, p2);

              subimageMins[0] = fminf(subimageMins[0], p2[0]);
              subimageMaxes[0] = fmaxf(subimageMaxes[0], p2[0]);

              subimageMins[1] = fminf(subimageMins[1], p2[1]);
              subimageMaxes[1] = fmaxf(subimageMaxes[1], p2[1]);
            }
          }
        }

        int subimageMinXY[2], subimageMaxXY[2];

        mScreen->Transform<int> (subimageMins, subimageMinXY);
        mScreen->Transform<int> (subimageMaxes, subimageMaxXY);

        // If the entire sub-image is off the screen, continue to the next cluster.
        if (subimageMinXY[0] >= mScreen->GetPixelsX() || subimageMaxXY[0] < 0 || subimageMinXY[1]
            >= mScreen->GetPixelsY() || subimageMaxXY[1] < 0)
        {
          continue;
        }

        // Crop the sub-image to the screen.
        subimageMinXY[0] = util::NumericalFunctions::max(subimageMinXY[0], 0);
        subimageMaxXY[0] = util::NumericalFunctions::min(subimageMaxXY[0], mScreen->GetPixelsX()
            - 1);
        subimageMinXY[1] = util::NumericalFunctions::max(subimageMinXY[1], 0);
        subimageMaxXY[1] = util::NumericalFunctions::min(subimageMaxXY[1], mScreen->GetPixelsY()
            - 1);

        AABB aabb;
        aabb.acc_1 = thisCluster->minmax_x[1] - viewpointCentre[0];
        aabb.acc_2 = thisCluster->minmax_x[0] - viewpointCentre[0];
        aabb.acc_3 = thisCluster->minmax_y[1] - viewpointCentre[1];
        aabb.acc_4 = thisCluster->minmax_y[0] - viewpointCentre[1];
        aabb.acc_5 = thisCluster->minmax_z[1] - viewpointCentre[2];
        aabb.acc_6 = thisCluster->minmax_z[0] - viewpointCentre[2];

        float par3[3];
        const float* vtx = mScreen->GetVtx();
        for (int l = 0; l < 3; l++)
        {
          par3[l] = vtx[l] + subimageMinXY[0] * projectedUnitX[l] + subimageMinXY[1]
              * projectedUnitY[l];
        }

        for (int subImageX = subimageMinXY[0]; subImageX <= subimageMaxXY[0]; ++subImageX)
        {
          float lRayDirection[3];
          for (int l = 0; l < 3; l++)
          {
            lRayDirection[l] = par3[l];
          }

          for (int subImageY = subimageMinXY[1]; subImageY <= subimageMaxXY[1]; ++subImageY)
          {
            Ray lRay;

            lRay.Direction[0] = lRayDirection[0];
            lRay.Direction[1] = lRayDirection[1];
            lRay.Direction[2] = lRayDirection[2];

            float lInverseDirectionMagnitude = 1.0F / sqrtf(lRayDirection[0] * lRayDirection[0]
                + lRayDirection[1] * lRayDirection[1] + lRayDirection[2] * lRayDirection[2]);

            lRay.Direction[0] *= lInverseDirectionMagnitude;
            lRay.Direction[1] *= lInverseDirectionMagnitude;
            lRay.Direction[2] *= lInverseDirectionMagnitude;

            lRay.InverseDirection[0] = 1.0F / lRay.Direction[0];
            lRay.InverseDirection[1] = 1.0F / lRay.Direction[1];
            lRay.InverseDirection[2] = 1.0F / lRay.Direction[2];

            bool lRayInPositiveDirection[3];
            lRayInPositiveDirection[0] = lRay.Direction[0] > 0.0F;
            lRayInPositiveDirection[1] = lRay.Direction[1] > 0.0F;
            lRayInPositiveDirection[2] = lRay.Direction[2] > 0.0F;

            lRayDirection[0] += projectedUnitY[0];
            lRayDirection[1] += projectedUnitY[1];
            lRayDirection[2] += projectedUnitY[2];

            float t_near, t_far;
            AABBvsRay(&aabb, lRay.InverseDirection, lRayInPositiveDirection, &t_near, &t_far);

            if (t_near > t_far)
            {
              continue;
            }

            float ray_dx[3];
            ray_dx[0] = t_near * lRay.Direction[0] - cluster_x[0];
            ray_dx[1] = t_near * lRay.Direction[1] - cluster_x[1];
            ray_dx[2] = t_near * lRay.Direction[2] - cluster_x[2];

            lRay.VelocityColour[0] = 0.0F;
            lRay.VelocityColour[1] = 0.0F;
            lRay.VelocityColour[2] = 0.0F;

            lRay.StressColour[0] = 0.0F;
            lRay.StressColour[1] = 0.0F;
            lRay.StressColour[2] = 0.0F;

            lRay.Length = 0.0F;
            lRay.MinT = std::numeric_limits<float>::max();
            lRay.Density = -1.0F;

            TraverseBlocks(thisCluster, lRayInPositiveDirection, ray_dx, block_flow_field, &lRay);

            if (lRay.MinT == std::numeric_limits<float>::max())
            {
              continue;
            }

            ColPixel col_pixel(subImageX,
                               subImageY,
                               lRay.MinT + t_near,
                               lRay.Length,
                                (lRay.Density - mDomainStats->density_threshold_min)
                                   * mDomainStats->density_threshold_minmax_inv,
                               lRay.Stress != std::numeric_limits<float>::max()
                                 ? lRay.Stress * mDomainStats->stress_threshold_max_inv
                                 : std::numeric_limits<float>::max(),
                               lRay.VelocityColour,
                               lRay.StressColour);

            mScreen->AddPixel(&col_pixel, mVisSettings);
          }
          par3[0] += projectedUnitX[0];
          par3[1] += projectedUnitX[1];
          par3[2] += projectedUnitX[2];
        }
      }
    }

    void RayTracer::UpdateClusterVoxel(int i, float density, float velocity, float stress)
    {
      cluster_voxel[3 * i][0] = density;
      cluster_voxel[3 * i][1] = velocity;
      cluster_voxel[3 * i][2] = stress;
    }

    RayTracer::~RayTracer()
    {
      for (unsigned int n = 0; n < mClusters.size(); n++)
      {
        for (int m = 0; m < (mClusters[n].blocks_x * mClusters[n].blocks_y * mClusters[n].blocks_z); m++)
        {
          if (cluster_flow_field[n][m] != NULL)
          {
            delete[] cluster_flow_field[n][m];
          }
        }
        delete[] cluster_flow_field[n];
      }

      delete[] cluster_flow_field;
      delete[] cluster_voxel;
    }

  }
}
