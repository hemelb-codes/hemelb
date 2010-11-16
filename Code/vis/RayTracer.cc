#include <math.h>
#include <stdlib.h>
#include <vector>
#include <limits>

#include "lb/LbmParameters.h"
#include "utilityFunctions.h"
#include "vis/RayTracer.h"
// TODO ideally this wouldn't be here.
#include "vis/Control.h"

namespace hemelb
{
  namespace vis
  {
    // TODO RENAME THIS FUNCTION AND MAKE IT MORE EFFICIENT.
    void RayTracer::rtAABBvsRayFn(const AABB &aabb,
                                  const float &inv_x,
                                  const float &inv_y,
                                  const float &inv_z,
                                  const bool xyz_sign_is_1[],
                                  float &t_near,
                                  float &t_far)
    {
      float tx0 = (xyz_sign_is_1[0]
        ? aabb.acc_2
        : aabb.acc_1) * inv_x;
      float tx1 = (xyz_sign_is_1[0]
        ? aabb.acc_1
        : aabb.acc_2) * inv_x;
      float ty0 = (xyz_sign_is_1[1]
        ? aabb.acc_4
        : aabb.acc_3) * inv_y;
      float ty1 = (xyz_sign_is_1[1]
        ? aabb.acc_3
        : aabb.acc_4) * inv_y;
      float tz0 = (xyz_sign_is_1[2]
        ? aabb.acc_6
        : aabb.acc_5) * inv_z;
      float tz1 = (xyz_sign_is_1[2]
        ? aabb.acc_5
        : aabb.acc_6) * inv_z;

      t_near = fmaxf(tx0, fmaxf(ty0, tz0));
      t_far = fminf(tx1, fminf(ty1, tz1));
    }

    void RayTracer::rtUpdateColour(float dt, float palette[], float col[])
    {
      col[0] += dt * palette[0];
      col[1] += dt * palette[1];
      col[2] += dt * palette[2];
    }

    void RayTracer::rtUpdateRayData(float *flow_field,
                                    float ray_t,
                                    float ray_segment,
                                    Ray *bCurrentRay,
                                    void(*ColourPalette)(float value,
                                                         float col[]),
                                    const lb::StressTypes iLbmStressType)
    {
      if (*flow_field < 0.0F)
        return; // solid voxel

      float palette[3];

      // update the volume rendering of the velocity flow field
      float scaled_velocity = * (flow_field + 1)
          * vis::controller->velocity_threshold_max_inv;

      ColourPalette(scaled_velocity, palette);

      bCurrentRay->VelocityColour[0] += ray_segment * palette[0];
      bCurrentRay->VelocityColour[1] += ray_segment * palette[1];
      bCurrentRay->VelocityColour[2] += ray_segment * palette[2];

      if (iLbmStressType != lb::ShearStress)
      {
        // update the volume rendering of the von Mises stress flow field
        float scaled_stress = * (flow_field + 2)
            * vis::controller->stress_threshold_max_inv;

        ColourPalette(scaled_stress, palette);

        bCurrentRay->StressColour[0] += ray_segment * palette[0];
        bCurrentRay->StressColour[1] += ray_segment * palette[1];
        bCurrentRay->StressColour[2] += ray_segment * palette[2];
      }
      bCurrentRay->Length += ray_segment;

      if (bCurrentRay->Density >= 0.0F)
        return;

      bCurrentRay->MinT = ray_t;

      // keep track of the density nearest to the view point
      bCurrentRay->Density = *flow_field;

      // keep track of the stress nearest to the view point
      bCurrentRay->Stress = * (flow_field + 2);
    }

    void RayTracer::rtTraverseVoxels(float block_min[],
                                     float block_x[],
                                     float voxel_flow_field[],
                                     float t,
                                     Ray *bCurrentRay,
                                     void(*ColourPalette)(float value,
                                                          float col[]),
                                     bool xyz_is_1[],
                                     const lb::StressTypes iLbmStressType)
    {
      float t_max[3];
      int i_vec[3];

      for (int i = 0; i < 3; i++)
      {
        i_vec[i] = (int) block_x[i];
      }
      for (int i = 0; i < 3; i++)
      {
        i_vec[i] = (i_vec[i] < 0)
          ? 0
          : i_vec[i];
        i_vec[i] = (i_vec[i] > block_size_1)
          ? block_size_1
          : i_vec[i];
      }

      for (int i = 0; i < 3; i++)
      {
        t_max[i] = (block_min[i] + (float) (xyz_is_1[i]
          ? i_vec[i] + 1
          : i_vec[i])) * bCurrentRay->InverseDirection[i];
      }

      int i = i_vec[0] * block_size2;
      int j = i_vec[1] * mGlobLatDat->BlockSize;
      int k = i_vec[2];

      while (true)
      {
        if (t_max[0] < t_max[1])
        {
          if (t_max[0] < t_max[2])
          {
            rtUpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS], t,
                            t_max[0] - t, bCurrentRay, ColourPalette,
                            iLbmStressType);

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
              if ( (i -= block_size2) < 0)
              {
                return;
              }
              t = t_max[0];
              t_max[0] -= bCurrentRay->InverseDirection[0];
            }
          }
          else
          {
            rtUpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS], t,
                            t_max[2] - t, bCurrentRay, ColourPalette,
                            iLbmStressType);

            if (xyz_is_1[2])
            {
              if (++k >= mGlobLatDat->BlockSize)
              {
                return;
              }
              t = t_max[2];
              t_max[2] += bCurrentRay->InverseDirection[2];
            }
            else
            {
              if (--k < 0)
              {
                return;
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
            rtUpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS], t,
                            t_max[1] - t, bCurrentRay, ColourPalette,
                            iLbmStressType);

            if (xyz_is_1[1])
            {
              if ( (j += mGlobLatDat->BlockSize) >= block_size2)
              {
                return;
              }
              t = t_max[1];
              t_max[1] += bCurrentRay->InverseDirection[1];
            }
            else
            {
              if ( (j -= mGlobLatDat->BlockSize) < 0)
              {
                return;
              }
              t = t_max[1];
              t_max[1] -= bCurrentRay->InverseDirection[1];
            }
          }
          else
          {
            rtUpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS], t,
                            t_max[2] - t, bCurrentRay, ColourPalette,
                            iLbmStressType);

            if (xyz_is_1[2])
            {
              if (++k >= mGlobLatDat->BlockSize)
              {
                return;
              }
              t = t_max[2];
              t_max[2] += bCurrentRay->InverseDirection[2];
            }
            else
            {
              if (--k < 0)
              {
                return;
              }
              t = t_max[2];
              t_max[2] -= bCurrentRay->InverseDirection[2];
            }
          }
        }
      }
    }

    void RayTracer::rtTraverseBlocksFn(float ray_dx[],
                                       float **block_flow_field,
                                       Ray *bCurrentRay,
                                       void(*ColourPalette)(float value,
                                                            float col[]),
                                       bool xyz_Is_1[],
                                       const lb::StressTypes iLbmStressType)
    {
      float block_min[3];
      float t_max[3];
      float block_x[3];
      float t_delta[3];
      float dx[3];

      int i_vec[3];

      for (int i = 0; i < 3; i++)
      {
        dx[i] = ray_dx[i];
      }
      for (int l = 0; l < 3; l++)
      {
        i_vec[l] = util::enforceBounds(cluster_blocks_vec[l], 0,
                                       (int) (mBlockSizeInverse * dx[l]));
        block_min[l] = (float) i_vec[l] * mBlockSizeFloat - dx[l];
      }
      int i = i_vec[0] * cluster_blocks_yz;
      int j = i_vec[1] * cluster_blocks_z;
      int k = i_vec[2];

      if (block_flow_field[i + j + k] != NULL)
      {
        block_x[0] = -block_min[0];
        block_x[1] = -block_min[1];
        block_x[2] = -block_min[2];

        rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k], 0.0F,
                         bCurrentRay, ColourPalette, xyz_Is_1, iLbmStressType);
      }
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

              rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k],
                               t_max[0], bCurrentRay, ColourPalette, xyz_Is_1,
                               iLbmStressType);
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

              rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k],
                               t_max[2], bCurrentRay, ColourPalette, xyz_Is_1,
                               iLbmStressType);
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

              rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k],
                               t_max[1], bCurrentRay, ColourPalette, xyz_Is_1,
                               iLbmStressType);
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

              rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k],
                               t_max[2], bCurrentRay, ColourPalette, xyz_Is_1,
                               iLbmStressType);
            }

            t_max[2] = xyz_Is_1[2]
              ? t_max[2] + t_delta[2]
              : t_max[2] - t_delta[2];
          }
        }
      }
    }

    void RayTracer::rtBuildClusters()
    {
      const int n_x[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, +0, +0, +0, +0,
                          +0, +0, +0, +0, +1, +1, +1, +1, +1, +1, +1, +1, +1 };
      const int n_y[] = { -1, -1, -1, +0, +0, +0, +1, +1, +1, -1, -1, -1, +0,
                          +0, +1, +1, +1, -1, -1, -1, +0, +0, +0, +1, +1, +1 };
      const int n_z[] = { -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1,
                          +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1 };

      unsigned short int *cluster_id =
          new unsigned short int[mGlobLatDat->BlockCount];
      for (int n = 0; n < mGlobLatDat->BlockCount; n++)
      {
        cluster_id[n] = -1;
      }

      mClusters = std::vector<Cluster*>();

      int blocks_buffer_size = 10000;

      std::vector<BlockLocation> *block_location_a = new std::vector<
          BlockLocation>(blocks_buffer_size);
      std::vector<BlockLocation> *block_location_b = new std::vector<
          BlockLocation>(blocks_buffer_size);

      bool *is_block_visited = new bool[mGlobLatDat->BlockCount];
      for (int n = 0; n < mGlobLatDat->BlockCount; n++)
      {
        is_block_visited[n] = false;
      }

      int block_min_x = std::numeric_limits<int>::max();
      int block_min_y = std::numeric_limits<int>::max();
      int block_min_z = std::numeric_limits<int>::max();
      int block_max_x = std::numeric_limits<int>::min();
      int block_max_y = std::numeric_limits<int>::min();
      int block_max_z = std::numeric_limits<int>::min();

      int n = -1;

      for (int i = 0; i < mGlobLatDat->BlocksX; i++)
      {
        for (int j = 0; j < mGlobLatDat->BlocksY; j++)
        {
          for (int k = 0; k < mGlobLatDat->BlocksZ; k++)
          {
            n++;

            lb::BlockData * lBlock = &mGlobLatDat->Blocks[n];
            if (lBlock->ProcessorRankForEachBlockSite == NULL)
            {
              continue;
            }

            block_min_x = util::min(block_min_x, i);
            block_min_y = util::min(block_min_y, j);
            block_min_z = util::min(block_min_z, k);
            block_max_x = util::max(block_max_x, i);
            block_max_y = util::max(block_max_y, j);
            block_max_z = util::max(block_max_z, k);

            if (is_block_visited[n])
            {
              continue;
            }

            is_block_visited[n] = 1;

            int blocks_a = 0;

            for (int m = 0; m < mGlobLatDat->SitesPerBlockVolumeUnit; m++)
            {
              if (mNetworkTopology->LocalRank
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

            Cluster *lNewCluster = new Cluster();

            lNewCluster->block_min[0] = i;
            lNewCluster->block_min[1] = j;
            lNewCluster->block_min[2] = k;

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
                const BlockLocation& tempBlockLoc =
                    block_location_a->at(index_a);

                for (int l = 0; l < 26; l++)
                {
                  int neigh_i = tempBlockLoc.i + n_x[l];
                  int neigh_j = tempBlockLoc.j + n_y[l];
                  int neigh_k = tempBlockLoc.k + n_z[l];

                  if (neigh_i == -1 || neigh_i == mGlobLatDat->BlocksX)
                    continue;
                  if (neigh_j == -1 || neigh_j == mGlobLatDat->BlocksY)
                    continue;
                  if (neigh_k == -1 || neigh_k == mGlobLatDat->BlocksZ)
                    continue;

                  int block_id = (neigh_i * mGlobLatDat->BlocksY + neigh_j)
                      * mGlobLatDat->BlocksZ + neigh_k;

                  if (is_block_visited[block_id]
                      || (lBlock = &mGlobLatDat->Blocks[block_id])->ProcessorRankForEachBlockSite
                          == NULL)
                  {
                    continue;
                  }

                  bool is_site_found = false;
                  for (int m = 0; m < mGlobLatDat->SitesPerBlockVolumeUnit; m++)
                  {
                    if (mNetworkTopology->LocalRank
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

                  lNewCluster->block_min[0]
                      = util::min((int) lNewCluster->block_min[0], neigh_i);
                  lNewCluster->block_min[1]
                      = util::min((int) lNewCluster->block_min[1], neigh_j);
                  lNewCluster->block_min[2]
                      = util::min((int) lNewCluster->block_min[2], neigh_k);

                  cluster_block_max_i = util::max((int) cluster_block_max_i,
                                                  neigh_i);
                  cluster_block_max_j = util::max((int) cluster_block_max_j,
                                                  neigh_j);
                  cluster_block_max_k = util::max((int) cluster_block_max_k,
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

            lNewCluster->x[0] = lNewCluster->block_min[0]
                * mGlobLatDat->BlockSize - 0.5F * mGlobLatDat->SitesX;
            lNewCluster->x[1] = lNewCluster->block_min[1]
                * mGlobLatDat->BlockSize - 0.5F * mGlobLatDat->SitesY;
            lNewCluster->x[2] = lNewCluster->block_min[2]
                * mGlobLatDat->BlockSize - 0.5F * mGlobLatDat->SitesZ;

            lNewCluster->blocks_x = 1 + cluster_block_max_i
                - lNewCluster->block_min[0];
            lNewCluster->blocks_y = 1 + cluster_block_max_j
                - lNewCluster->block_min[1];
            lNewCluster->blocks_z = 1 + cluster_block_max_k
                - lNewCluster->block_min[2];

            mClusters.push_back(lNewCluster);
          }
        }
      }

      delete block_location_b;
      delete block_location_a;

      delete[] is_block_visited;

      vis::controller->ctr_x = 0.5F * mGlobLatDat->BlockSize * (block_min_x
          + block_max_x);
      vis::controller->ctr_y = 0.5F * mGlobLatDat->BlockSize * (block_min_y
          + block_max_y);
      vis::controller->ctr_z = 0.5F * mGlobLatDat->BlockSize * (block_min_z
          + block_max_z);

      cluster_voxel = new float *[mLocalLatDat->LocalFluidSites * VIS_FIELDS];

      cluster_flow_field = new float **[mClusters.size()];

      for (unsigned int lThisClusterId = 0; lThisClusterId < mClusters.size(); lThisClusterId++)
      {
        Cluster *cluster_p = mClusters[lThisClusterId];

        cluster_flow_field[lThisClusterId] = new float *[cluster_p->blocks_x
            * cluster_p->blocks_y * cluster_p->blocks_z];

        int voxel_min[3], voxel_max[3];
        for (int l = 0; l < 3; l++)
        {
          voxel_min[l] = std::numeric_limits<int>::max();
          voxel_max[l] = std::numeric_limits<int>::min();
        }

        n = -1;

        int block_coord[3];
        for (int i = 0; i < cluster_p->blocks_x; i++)
        {
          block_coord[0] = (i + cluster_p->block_min[0])
              * mGlobLatDat->BlockSize;

          for (int j = 0; j < cluster_p->blocks_y; j++)
          {
            block_coord[1] = (j + cluster_p->block_min[1])
                * mGlobLatDat->BlockSize;

            for (int k = 0; k < cluster_p->blocks_z; k++)
            {
              block_coord[2] = (k + cluster_p->block_min[2])
                  * mGlobLatDat->BlockSize;

              int block_id = ( (i + cluster_p->block_min[0])
                  * mGlobLatDat->BlocksY + (j + cluster_p->block_min[1]))
                  * mGlobLatDat->BlocksZ + (k + cluster_p->block_min[2]);

              cluster_flow_field[lThisClusterId][++n] = NULL;

              if (cluster_id[block_id] != lThisClusterId)
              {
                continue;
              }

              lb::BlockData * lBlock = &mGlobLatDat->Blocks[block_id];

              cluster_flow_field[lThisClusterId][n]
                  = new float[mGlobLatDat->SitesPerBlockVolumeUnit * VIS_FIELDS];

              int m = -1;

              float ii[3];
              for (ii[0] = 0; ii[0] < mGlobLatDat->BlockSize; ii[0]++)
                for (ii[1] = 0; ii[1] < mGlobLatDat->BlockSize; ii[1]++)
                  for (ii[2] = 0; ii[2] < mGlobLatDat->BlockSize; ii[2]++)
                  {
                    unsigned int my_site_id;
                    my_site_id = lBlock->site_data[++m];

                    if (my_site_id & (1U << 31U))
                    {
                      for (int l = 0; l < VIS_FIELDS; l++)
                        cluster_flow_field[lThisClusterId][n][m * VIS_FIELDS
                            + l] = -1.0F;

                      continue;
                    }

                    for (int l = 0; l < VIS_FIELDS; l++)
                    {
                      cluster_flow_field[lThisClusterId][n][m * VIS_FIELDS + l]
                          = 1.0F;
                    }

                    for (int l = 0; l < VIS_FIELDS; l++)
                    {
                      cluster_voxel[my_site_id * VIS_FIELDS + l]
                          = &cluster_flow_field[lThisClusterId][n][m
                              * VIS_FIELDS + l];
                    }

                    for (int l = 0; l < 3; l++)
                    {
                      voxel_min[l] = util::min(voxel_min[l], ii[l]
                          + block_coord[l]);
                      voxel_max[l] = util::max(voxel_max[l], ii[l]
                          + block_coord[l]);
                    }

                  } // for ii[0..2]

            } // for k
          } // for j
        } // for i

        cluster_p->minmax_x[0] = (float) voxel_min[0] - 0.5F
            * (float) mGlobLatDat->SitesX;
        cluster_p->minmax_y[0] = (float) voxel_min[1] - 0.5F
            * (float) mGlobLatDat->SitesY;
        cluster_p->minmax_z[0] = (float) voxel_min[2] - 0.5F
            * (float) mGlobLatDat->SitesZ;

        cluster_p->minmax_x[1] = (float) (voxel_max[0] + 1) - 0.5F
            * (float) mGlobLatDat->SitesX;
        cluster_p->minmax_y[1] = (float) (voxel_max[1] + 1) - 0.5F
            * (float) mGlobLatDat->SitesY;
        cluster_p->minmax_z[1] = (float) (voxel_max[2] + 1) - 0.5F
            * (float) mGlobLatDat->SitesZ;
      }
      delete[] cluster_id;
    }

    RayTracer::RayTracer(const topology::NetworkTopology * iNetworkTopology,
                         const lb::LocalLatticeData* iLocalLatDat,
                         const lb::GlobalLatticeData* iGlobLatDat)
    {
      mNetworkTopology = iNetworkTopology;
      mLocalLatDat = iLocalLatDat;
      mGlobLatDat = iGlobLatDat;

      // Init globals
      blocks_yz = iGlobLatDat->BlocksY * iGlobLatDat->BlocksZ;
      mBlockSizeFloat = float(iGlobLatDat->BlockSize);
      block_size2 = iGlobLatDat ->BlockSize * iGlobLatDat->BlockSize;
      block_size3 = iGlobLatDat->BlockSize * block_size2;
      block_size_1 = iGlobLatDat->BlockSize - 1;

      mBlockSizeInverse = 1.F / mBlockSizeFloat;

      rtBuildClusters();
    }

    void RayTracer::Render(const lb::StressTypes iLbmStressType)
    {
      float lScreenMaxX = vis::controller->mScreen.MaxXValue;
      float lScreenMaxY = vis::controller->mScreen.MaxYValue;

      int lPixelsX = vis::controller->mScreen.PixelsX;
      int lPixelsY = vis::controller->mScreen.PixelsY;

      float p0[3];
      for (int l = 0; l < 3; l++)
      {
        p0[l] = vis::controller->mViewpoint.x[l];
      }

      float par1[3], par2[3];
      float screen_vtx[3];
      for (int l = 0; l < 3; l++)
      {
        par1[l] = vis::controller->mScreen.UnitVectorProjectionX[l];
        par2[l] = vis::controller->mScreen.UnitVectorProjectionY[l];
        screen_vtx[l] = vis::controller->mScreen.vtx[l];
      }

      float lScaleX = vis::controller->mScreen.ScaleX;
      float lScaleY = vis::controller->mScreen.ScaleY;

      for (unsigned int cluster_id = 0; cluster_id < mClusters.size(); cluster_id++)
      {
        Cluster *cluster_p = mClusters[cluster_id];

        // the image-based projection of the mClusters bounding box is
        // calculated here
        float cluster_x[3];
        for (int l = 0; l < 3; l++)
        {
          cluster_x[l] = cluster_p->x[l] - p0[l];
        }
        cluster_blocks_vec[0] = cluster_p->blocks_x - 1;
        cluster_blocks_vec[1] = cluster_p->blocks_y - 1;
        cluster_blocks_vec[2] = cluster_p->blocks_z - 1;
        cluster_blocks_z = cluster_p->blocks_z;
        cluster_blocks_yz = (int) cluster_p->blocks_y
            * (int) cluster_p->blocks_z;
        cluster_blocks = (int) cluster_p->blocks_x * cluster_blocks_yz;

        float **block_flow_field = cluster_flow_field[cluster_id];

        float subimage_vtx[4];
        subimage_vtx[0] = 1.0e+30F;
        subimage_vtx[1] = -1.0e+30F;
        subimage_vtx[2] = 1.0e+30F;
        subimage_vtx[3] = -1.0e+30F;

        float p1[3];
        for (int i = 0; i < 2; i++)
        {
          p1[0] = cluster_p->minmax_x[i];

          for (int j = 0; j < 2; j++)
          {
            p1[1] = cluster_p->minmax_y[j];

            for (int k = 0; k < 2; k++)
            {
              p1[2] = cluster_p->minmax_z[k];

              float p2[3];
              vis::controller->project(p1, p2);

              subimage_vtx[0] = fminf(subimage_vtx[0], p2[0]);
              subimage_vtx[1] = fmaxf(subimage_vtx[1], p2[0]);
              subimage_vtx[2] = fminf(subimage_vtx[2], p2[1]);
              subimage_vtx[3] = fmaxf(subimage_vtx[3], p2[1]);
            }
          }
        }

        int subimage_pix[4];
        subimage_pix[0] = (int) (lScaleX * (subimage_vtx[0] + lScreenMaxX)
            + 0.5F);
        subimage_pix[1] = (int) (lScaleX * (subimage_vtx[1] + lScreenMaxX)
            + 0.5F);
        subimage_pix[2] = (int) (lScaleY * (subimage_vtx[2] + lScreenMaxY)
            + 0.5F);
        subimage_pix[3] = (int) (lScaleY * (subimage_vtx[3] + lScreenMaxY)
            + 0.5F);

        if (subimage_pix[0] >= lPixelsX || subimage_pix[1] < 0
            || subimage_pix[2] >= lPixelsY || subimage_pix[3] < 0)
        {
          continue;
        }
        subimage_pix[0] = util::max(subimage_pix[0], 0);
        subimage_pix[1] = util::min(subimage_pix[1], lPixelsX - 1);
        subimage_pix[2] = util::max(subimage_pix[2], 0);
        subimage_pix[3] = util::min(subimage_pix[3], lPixelsY - 1);

        AABB aabb;
        aabb.acc_1 = cluster_p->minmax_x[1] - p0[0];
        aabb.acc_2 = cluster_p->minmax_x[0] - p0[0];
        aabb.acc_3 = cluster_p->minmax_y[1] - p0[1];
        aabb.acc_4 = cluster_p->minmax_y[0] - p0[1];
        aabb.acc_5 = cluster_p->minmax_z[1] - p0[2];
        aabb.acc_6 = cluster_p->minmax_z[0] - p0[2];

        float par3[3];
        for (int l = 0; l < 3; l++)
        {
          par3[l] = screen_vtx[l] + subimage_pix[0] * par1[l] + subimage_pix[2]
              * par2[l];
        }

        for (int i = subimage_pix[0]; i <= subimage_pix[1]; i++)
        {
          float lRayDirection[3];
          for (int l = 0; l < 3; l++)
          {
            lRayDirection[l] = par3[l];
          }
          for (int j = subimage_pix[2]; j <= subimage_pix[3]; j++)
          {
            Ray lRay;

            lRay.Direction[0] = lRayDirection[0];
            lRay.Direction[1] = lRayDirection[1];
            lRay.Direction[2] = lRayDirection[2];

            float lInverseDirectionMagnitude = 1.0F / sqrtf(lRayDirection[0]
                * lRayDirection[0] + lRayDirection[1] * lRayDirection[1]
                + lRayDirection[2] * lRayDirection[2]);

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

            lRayDirection[0] += par2[0];
            lRayDirection[1] += par2[1];
            lRayDirection[2] += par2[2];

            float t_near, t_far;
            rtAABBvsRayFn(aabb, lRay.InverseDirection[0],
                          lRay.InverseDirection[1], lRay.InverseDirection[2],
                          lRayInPositiveDirection, t_near, t_far);

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

            if (iLbmStressType != lb::ShearStress)
            {
              lRay.StressColour[0] = 0.0F;
              lRay.StressColour[1] = 0.0F;
              lRay.StressColour[2] = 0.0F;
            }
            lRay.Length = 0.0F;
            lRay.MinT = std::numeric_limits<float>::max();
            lRay.Density = -1.0F;

            rtTraverseBlocksFn(ray_dx, block_flow_field, &lRay,
                               ColourPalette::pickColour,
                               lRayInPositiveDirection, iLbmStressType);

            if (lRay.MinT == std::numeric_limits<float>::max())
            {
              continue;
            }

            ColPixel col_pixel;
            col_pixel.vel_r = lRay.VelocityColour[0] * 255.0F;
            col_pixel.vel_g = lRay.VelocityColour[1] * 255.0F;
            col_pixel.vel_b = lRay.VelocityColour[2] * 255.0F;

            if (iLbmStressType != lb::ShearStress)
            {
              col_pixel.stress_r = lRay.StressColour[0] * 255.0F;
              col_pixel.stress_g = lRay.StressColour[1] * 255.0F;
              col_pixel.stress_b = lRay.StressColour[2] * 255.0F;
            }
            col_pixel.dt = lRay.Length;
            col_pixel.t = lRay.MinT + t_near;
            col_pixel.density = (lRay.Density
                - vis::controller->density_threshold_min)
                * vis::controller->density_threshold_minmax_inv;

            if (lRay.Stress != std::numeric_limits<float>::max())
            {
              col_pixel.stress = lRay.Stress
                  * vis::controller->stress_threshold_max_inv;
            }
            else
            {
              col_pixel.stress = std::numeric_limits<float>::max();
            }
            col_pixel.i = PixelId(i, j);
            col_pixel.i.isRt = true;

            vis::controller->writePixel(&col_pixel);
          }
          par3[0] += par1[0];
          par3[1] += par1[1];
          par3[2] += par1[2];
        }
      }
    }

    void RayTracer::UpdateClusterVoxel(const int &i,
                                       const float &density,
                                       const float &velocity,
                                       const float &stress)
    {
      cluster_voxel[3 * i][0] = density;
      cluster_voxel[3 * i][1] = velocity;
      cluster_voxel[3 * i][2] = stress;
    }

    RayTracer::~RayTracer()
    {
      for (unsigned int n = 0; n < mClusters.size(); n++)
      {
        for (int m = 0; m < (mClusters[n]->blocks_x * mClusters[n]->blocks_y
            * mClusters[n]->blocks_z); m++)
        {
          if (cluster_flow_field[n][m] != NULL)
          {
            delete[] cluster_flow_field[n][m];
          }
        }
        delete[] cluster_flow_field[n];

        delete mClusters[n];
      }

      delete[] cluster_flow_field;
      delete[] cluster_voxel;

      mClusters.clear();
    }

  }
}
