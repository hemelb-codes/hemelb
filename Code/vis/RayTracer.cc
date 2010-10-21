#include <math.h>
#include <stdlib.h>
#include <vector>

#include "utilityFunctions.h"

#include "vis/RayTracer.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace vis
  {

    //TODO put these somewhere more sensible.
    std::vector<Cluster> cluster;
    float **cluster_voxel = NULL;
    float ***cluster_flow_field = NULL;

    int clusters;

    float ray_dir[3];
    float ray_inv[3];
    float ray_vel_col[3];
    float ray_stress_col[3];
    float ray_length;
    float ray_t_min;
    float ray_density;
    float ray_stress;

    int cluster_blocks_vec[3];
    int cluster_blocks_z, cluster_blocks_yz, cluster_blocks;

    float block_size_f;
    float block_size_inv;
    int block_size2, block_size3, block_size_1;
    int blocks_yz;

    // TODO RENAME THIS FUNCTION AND MAKE IT MORE EFFICIENT.
    void rtAABBvsRayFn(AABB *aabb,
                       float inv_x,
                       float inv_y,
                       float inv_z,
                       float *t_near,
                       float *t_far,
                       bool xyz_sign_is_1[])
    {
      float tx0, ty0, tz0;
      float tx1, ty1, tz1;

      tx0 = (xyz_sign_is_1[0]
        ? aabb->acc_2
        : aabb->acc_1) * inv_x;
      tx1 = (xyz_sign_is_1[0]
        ? aabb->acc_1
        : aabb->acc_2) * inv_x;
      ty0 = (xyz_sign_is_1[1]
        ? aabb->acc_4
        : aabb->acc_3) * inv_y;
      ty1 = (xyz_sign_is_1[1]
        ? aabb->acc_3
        : aabb->acc_4) * inv_y;
      tz0 = (xyz_sign_is_1[2]
        ? aabb->acc_6
        : aabb->acc_5) * inv_z;
      tz1 = (xyz_sign_is_1[2]
        ? aabb->acc_5
        : aabb->acc_6) * inv_z;

      *t_near = fmaxf(tx0, fmaxf(ty0, tz0));
      *t_far = fminf(tx1, fminf(ty1, tz1));
    }

    void rtUpdateColour(float dt, float palette[], float col[])
    {
      col[0] += dt * palette[0];
      col[1] += dt * palette[1];
      col[2] += dt * palette[2];
    }

    void rayTracer::rtUpdateRayData(float *flow_field,
                                    float ray_t,
                                    float ray_segment,
                                    void(*ColourPalette)(float value,
                                                         float col[]),
                                    const float iLbmStressType)
    {
      if (*flow_field < 0.0F)
        return; // solid voxel

      float palette[3];

      // update the volume rendering of the velocity flow field
      float scaled_velocity = * (flow_field + 1)
          * vis::controller->velocity_threshold_max_inv;

      ColourPalette(scaled_velocity, palette);

      ray_vel_col[0] += ray_segment * palette[0];
      ray_vel_col[1] += ray_segment * palette[1];
      ray_vel_col[2] += ray_segment * palette[2];

      if (iLbmStressType != SHEAR_STRESS)
      {
        // update the volume rendering of the von Mises stress flow field
        float scaled_stress = * (flow_field + 2)
            * vis::controller->stress_threshold_max_inv;

        ColourPalette(scaled_stress, palette);

        ray_stress_col[0] += ray_segment * palette[0];
        ray_stress_col[1] += ray_segment * palette[1];
        ray_stress_col[2] += ray_segment * palette[2];
      }
      ray_length += ray_segment;

      if (ray_density >= 0.0F)
        return;

      ray_t_min = ray_t;

      // keep track of the density nearest to the view point
      ray_density = *flow_field;

      // keep track of the stress nearest to the view point
      ray_stress = * (flow_field + 2);
    }

    void rayTracer::rtTraverseVoxels(float block_min[],
                                     float block_x[],
                                     float voxel_flow_field[],
                                     float t,
                                     void(*ColourPalette)(float value,
                                                          float col[]),
                                     bool xyz_is_1[],
                                     const float iLbmStressType)
    {
      float t_max[3];
      int i_vec[3];
      int i, j, k;

      for (i = 0; i < 3; i++)
      {
        i_vec[i] = (int) block_x[i];
      }
      for (i = 0; i < 3; i++)
      {
        i_vec[i] = (i_vec[i] < 0)
          ? 0
          : i_vec[i];
        i_vec[i] = (i_vec[i] > block_size_1)
          ? block_size_1
          : i_vec[i];
      }

      for (i = 0; i < 3; i++)
      {
        t_max[i] = (block_min[i] + (float) (xyz_is_1[i]
          ? i_vec[i] + 1
          : i_vec[i])) * ray_inv[i];
      }

      i = i_vec[0] * block_size2;
      j = i_vec[1] * block_size;
      k = i_vec[2];

      while (true)
      {
        if (t_max[0] < t_max[1])
        {
          if (t_max[0] < t_max[2])
          {
            rtUpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS], t,
                            t_max[0] - t, ColourPalette, iLbmStressType);

            if (xyz_is_1[0])
            {
              if ( (i += block_size2) >= block_size3)
                return;
              t = t_max[0];
              t_max[0] += ray_inv[0];
            }
            else
            {
              if ( (i -= block_size2) < 0)
                return;
              t = t_max[0];
              t_max[0] -= ray_inv[0];
            }
          }
          else
          {
            rtUpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS], t,
                            t_max[2] - t, ColourPalette, iLbmStressType);

            if (xyz_is_1[2])
            {
              if (++k >= block_size)
                return;
              t = t_max[2];
              t_max[2] += ray_inv[2];
            }
            else
            {
              if (--k < 0)
                return;
              t = t_max[2];
              t_max[2] -= ray_inv[2];
            }
          }
        }
        else
        {
          if (t_max[1] < t_max[2])
          {
            rtUpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS], t,
                            t_max[1] - t, ColourPalette, iLbmStressType);

            if (xyz_is_1[1])
            {
              if ( (j += block_size) >= block_size2)
                return;
              t = t_max[1];
              t_max[1] += ray_inv[1];
            }
            else
            {
              if ( (j -= block_size) < 0)
                return;
              t = t_max[1];
              t_max[1] -= ray_inv[1];
            }
          }
          else
          {
            rtUpdateRayData(&voxel_flow_field[ (i + j + k) * VIS_FIELDS], t,
                            t_max[2] - t, ColourPalette, iLbmStressType);

            if (xyz_is_1[2])
            {
              if (++k >= block_size)
                return;
              t = t_max[2];
              t_max[2] += ray_inv[2];
            }
            else
            {
              if (--k < 0)
                return;
              t = t_max[2];
              t_max[2] -= ray_inv[2];
            }
          }
        }
      }
    }

    void rayTracer::rtTraverseBlocksFn(float ray_dx[],
                                       float **block_flow_field,
                                       void(*ColourPalette)(float value,
                                                            float col[]),
                                       bool xyz_Is_1[],
                                       const double iLbmStressType)
    {
      float block_min[3];
      float t_max[3];
      float block_x[3];
      float t_delta[3];
      float dx[3];

      int i_vec[3];
      int i, j, k, l;

      for (i = 0; i < 3; i++)
      {
        dx[i] = ray_dx[i];
      }
      for (l = 0; l < 3; l++)
      {
        i_vec[l] = util::enforceBounds(cluster_blocks_vec[l], 0,
                                       (int) (block_size_inv * dx[l]));
        block_min[l] = (float) i_vec[l] * block_size_f - dx[l];
      }
      i = i_vec[0] * cluster_blocks_yz;
      j = i_vec[1] * cluster_blocks_z;
      k = i_vec[2];

      if (block_flow_field[i + j + k] != NULL)
      {
        block_x[0] = -block_min[0];
        block_x[1] = -block_min[1];
        block_x[2] = -block_min[2];

        rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k], 0.0F,
                         ColourPalette, xyz_Is_1, iLbmStressType);
      }
      for (l = 0; l < 3; l++)
      {
        t_max[l] = (xyz_Is_1[l]
          ? block_min[l] + block_size_f
          : block_min[l]) * ray_inv[l];
        t_delta[l] = block_size_f * ray_inv[l];
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
              block_min[0] += block_size_f;
            }
            else
            {
              if ( (i -= cluster_blocks_yz) < 0)
                return;
              block_min[0] -= block_size_f;
            }

            if (block_flow_field[i + j + k] != NULL)
            {
              block_x[0] = t_max[0] * ray_dir[0] - block_min[0];
              block_x[1] = t_max[0] * ray_dir[1] - block_min[1];
              block_x[2] = t_max[0] * ray_dir[2] - block_min[2];

              rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k],
                               t_max[0], ColourPalette, xyz_Is_1,
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
              block_min[2] += block_size_f;
            }
            else
            {
              if (--k < 0)
                return;
              block_min[2] -= block_size_f;
            }

            if (block_flow_field[i + j + k] != NULL)
            {
              block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
              block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
              block_x[2] = t_max[2] * ray_dir[2] - block_min[2];

              rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k],
                               t_max[2], ColourPalette, xyz_Is_1,
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
              block_min[1] += block_size_f;
            }
            else
            {
              if ( (j -= cluster_blocks_z) < 0)
                return;
              block_min[1] -= block_size_f;
            }

            if (block_flow_field[i + j + k] != NULL)
            {
              block_x[0] = t_max[1] * ray_dir[0] - block_min[0];
              block_x[1] = t_max[1] * ray_dir[1] - block_min[1];
              block_x[2] = t_max[1] * ray_dir[2] - block_min[2];

              rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k],
                               t_max[1], ColourPalette, xyz_Is_1,
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
              block_min[2] += block_size_f;
            }
            else
            {
              if (--k < 0)
                return;
              block_min[2] -= block_size_f;
            }

            if (block_flow_field[i + j + k] != NULL)
            {
              block_x[0] = t_max[2] * ray_dir[0] - block_min[0];
              block_x[1] = t_max[2] * ray_dir[1] - block_min[1];
              block_x[2] = t_max[2] * ray_dir[2] - block_min[2];

              rtTraverseVoxels(block_min, block_x, block_flow_field[i + j + k],
                               t_max[2], ColourPalette, xyz_Is_1,
                               iLbmStressType);
            }

            t_max[2] = xyz_Is_1[2]
              ? t_max[2] + t_delta[2]
              : t_max[2] - t_delta[2];
          }
        }
      }
    }

    void rtBuildClusters(Net *net)
    {
      int n_x[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, +0, +0, +0, +0, +0, +0,
                    +0, +0, +1, +1, +1, +1, +1, +1, +1, +1, +1 };
      int n_y[] = { -1, -1, -1, +0, +0, +0, +1, +1, +1, -1, -1, -1, +0, +0, +1,
                    +1, +1, -1, -1, -1, +0, +0, +0, +1, +1, +1 };
      int n_z[] = { -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1, -1, +1, -1,
                    +0, +1, -1, +0, +1, -1, +0, +1, -1, +0, +1 };

      int neigh_i, neigh_j, neigh_k;
      int voxel_min[3], voxel_max[3];
      int block_coord[3], ii[3];
      int cluster_id;
      int block_id;
      int blocks_a, blocks_b;
      int index_a;
      int blocks_buffer_size;
      int are_blocks_incrementing;
      int is_site_found;
      int i, j, k;
      int l, m, n;
      int clusters_max;
      int block_min_x, block_min_y, block_min_z;
      int block_max_x, block_max_y, block_max_z;
      int dummy = 0;

      unsigned int my_site_id;

      unsigned short int cluster_block_max_i, cluster_block_max_j,
          cluster_block_max_k;

      bool *is_block_visited;

      std::vector<BlockLocation> *block_location_a, *block_location_b;

      DataBlock *map_block_p;
      ProcBlock *proc_block_p;

      Cluster *cluster_p;

      net->cluster_id = new short int[blocks];

      clusters_max = 20;
      clusters = 0;
      cluster = std::vector<Cluster>(clusters_max);

      for (n = 0; n < blocks; n++)
      {
        net->cluster_id[n] = -1;
      }

      cluster_p = NULL;
      cluster_block_max_i = dummy;
      cluster_block_max_j = dummy;
      cluster_block_max_k = dummy;

      is_block_visited = new bool[blocks];

      blocks_buffer_size = 10000;
      block_location_a = new std::vector<BlockLocation>(blocks_buffer_size);
      block_location_b = new std::vector<BlockLocation>(blocks_buffer_size);

      for (n = 0; n < blocks; n++)
      {
        is_block_visited[n] = 0;
      }

      block_min_x = +1000000000;
      block_min_y = +1000000000;
      block_min_z = +1000000000;
      block_max_x = -1000000000;
      block_max_y = -1000000000;
      block_max_z = -1000000000;

      n = -1;

      for (i = 0; i < blocks_x; i++)
      {
        for (j = 0; j < blocks_y; j++)
        {
          for (k = 0; k < blocks_z; k++)
          {
            if ( (proc_block_p = &net->proc_block[++n])->proc_id == NULL)
              continue;

            block_min_x = util::min(block_min_x, i);
            block_min_y = util::min(block_min_y, j);
            block_min_z = util::min(block_min_z, k);
            block_max_x = util::max(block_max_x, i);
            block_max_y = util::max(block_max_y, j);
            block_max_z = util::max(block_max_z, k);

            if (is_block_visited[n])
              continue;

            is_block_visited[n] = 1;

            blocks_a = 0;

            for (m = 0; m < sites_in_a_block; m++)
            {
              if (net->IsCurrentProcRank(proc_block_p->proc_id[m]))
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
              continue;

            if (clusters == clusters_max)
            {
              clusters_max <<= 1;
              cluster.resize(clusters_max);
            }

            cluster_p = &cluster[clusters];
            ++clusters;

            cluster_p->block_min[0] = i;
            cluster_p->block_min[1] = j;
            cluster_p->block_min[2] = k;

            cluster_block_max_i = i;
            cluster_block_max_j = j;
            cluster_block_max_k = k;

            net->cluster_id[n] = clusters - 1;

            are_blocks_incrementing = 1;

            while (are_blocks_incrementing)
            {
              blocks_b = 0;
              are_blocks_incrementing = 0;

              for (index_a = 0; index_a < blocks_a; index_a++)
              {
                const BlockLocation& tempBlockLoc =
                    block_location_a->at(index_a);

                for (l = 0; l < 26; l++)
                {
                  neigh_i = tempBlockLoc.i + n_x[l];
                  neigh_j = tempBlockLoc.j + n_y[l];
                  neigh_k = tempBlockLoc.k + n_z[l];

                  if (neigh_i == -1 || neigh_i == blocks_x)
                    continue;
                  if (neigh_j == -1 || neigh_j == blocks_y)
                    continue;
                  if (neigh_k == -1 || neigh_k == blocks_z)
                    continue;

                  block_id = (neigh_i * blocks_y + neigh_j) * blocks_z
                      + neigh_k;

                  if (is_block_visited[block_id] || (proc_block_p
                      = &net->proc_block[block_id])->proc_id == NULL)
                  {
                    continue;
                  }

                  is_site_found = 0;

                  for (m = 0; m < sites_in_a_block; m++)
                  {
                    if (net->IsCurrentProcRank(proc_block_p->proc_id[m]))
                    {
                      is_site_found = 1;
                      break;
                    }
                  }

                  if (!is_site_found)
                    continue;

                  is_block_visited[block_id] = 1;

                  are_blocks_incrementing = 1;

                  if (blocks_b == blocks_buffer_size)
                  {
                    blocks_buffer_size *= 2;

                    block_location_a->resize(blocks_buffer_size);
                    block_location_b->resize(blocks_buffer_size);
                  }

                  BlockLocation& tempBlockLoc = block_location_b->at(blocks_b);
                  tempBlockLoc.i = neigh_i;
                  tempBlockLoc.j = neigh_j;
                  tempBlockLoc.k = neigh_k;
                  ++blocks_b;

                  cluster_p->block_min[0]
                      = util::min((int) cluster_p->block_min[0], neigh_i);
                  cluster_p->block_min[1]
                      = util::min((int) cluster_p->block_min[1], neigh_j);
                  cluster_p->block_min[2]
                      = util::min((int) cluster_p->block_min[2], neigh_k);

                  cluster_block_max_i = util::max((int) cluster_block_max_i,
                                                  neigh_i);
                  cluster_block_max_j = util::max((int) cluster_block_max_j,
                                                  neigh_j);
                  cluster_block_max_k = util::max((int) cluster_block_max_k,
                                                  neigh_k);

                  net->cluster_id[block_id] = clusters - 1;
                }
              }

              // swap pointers in block_location_a/_b
              std::vector<BlockLocation>* tempBlockLocation = block_location_a;
              block_location_a = block_location_b;
              block_location_b = tempBlockLocation;

              blocks_a = blocks_b;
            }

            cluster_p->x[0] = cluster_p->block_min[0] * block_size - 0.5F
                * sites_x;
            cluster_p->x[1] = cluster_p->block_min[1] * block_size - 0.5F
                * sites_y;
            cluster_p->x[2] = cluster_p->block_min[2] * block_size - 0.5F
                * sites_z;

            cluster_p->blocks_x = 1 + cluster_block_max_i
                - cluster_p->block_min[0];
            cluster_p->blocks_y = 1 + cluster_block_max_j
                - cluster_p->block_min[1];
            cluster_p->blocks_z = 1 + cluster_block_max_k
                - cluster_p->block_min[2];
          }
        }
      }

      delete block_location_b;
      delete block_location_a;

      delete[] is_block_visited;

      vis::controller->ctr_x = 0.5F * block_size * (block_min_x + block_max_x);
      vis::controller->ctr_y = 0.5F * block_size * (block_min_y + block_max_y);
      vis::controller->ctr_z = 0.5F * block_size * (block_min_z + block_max_z);

      cluster_voxel = new float *[net->my_sites * VIS_FIELDS];

      cluster_flow_field = new float **[clusters];

      for (cluster_id = 0; cluster_id < clusters; cluster_id++)
      {
        cluster_p = &cluster[cluster_id];

        cluster_flow_field[cluster_id] = new float *[cluster_p->blocks_x
            * cluster_p->blocks_y * cluster_p->blocks_z];

        for (l = 0; l < 3; l++)
        {
          voxel_min[l] = +1000000000;
          voxel_max[l] = -1000000000;
        }

        n = -1;

        for (i = 0; i < cluster_p->blocks_x; i++)
        {
          block_coord[0] = (i + cluster_p->block_min[0]) * block_size;

          for (j = 0; j < cluster_p->blocks_y; j++)
          {
            block_coord[1] = (j + cluster_p->block_min[1]) * block_size;

            for (k = 0; k < cluster_p->blocks_z; k++)
            {
              block_coord[2] = (k + cluster_p->block_min[2]) * block_size;

              block_id = ( (i + cluster_p->block_min[0]) * blocks_y + (j
                  + cluster_p->block_min[1])) * blocks_z + (k
                  + cluster_p->block_min[2]);

              cluster_flow_field[cluster_id][++n] = NULL;

              if (net->cluster_id[block_id] != cluster_id)
              {
                continue;
              }

              map_block_p = &net->map_block[block_id];

              cluster_flow_field[cluster_id][n] = new float[sites_in_a_block
                  * VIS_FIELDS];

              m = -1;

              for (ii[0] = 0; ii[0] < block_size; ii[0]++)
                for (ii[1] = 0; ii[1] < block_size; ii[1]++)
                  for (ii[2] = 0; ii[2] < block_size; ii[2]++)
                  {

                    my_site_id = map_block_p->site_data[++m];

                    if (my_site_id & (1U << 31U))
                    {
                      for (l = 0; l < VIS_FIELDS; l++)
                        cluster_flow_field[cluster_id][n][m * VIS_FIELDS + l]
                            = -1.0F;

                      continue;
                    }

                    for (l = 0; l < VIS_FIELDS; l++)
                      cluster_flow_field[cluster_id][n][m * VIS_FIELDS + l]
                          = 1.0F;

                    for (l = 0; l < VIS_FIELDS; l++)
                      cluster_voxel[my_site_id * VIS_FIELDS + l]
                          = &cluster_flow_field[cluster_id][n][m * VIS_FIELDS
                              + l];

                    for (l = 0; l < 3; l++)
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

        cluster_p->minmax_x[0] = (float) voxel_min[0] - 0.5F * (float) sites_x;
        cluster_p->minmax_y[0] = (float) voxel_min[1] - 0.5F * (float) sites_y;
        cluster_p->minmax_z[0] = (float) voxel_min[2] - 0.5F * (float) sites_z;

        cluster_p->minmax_x[1] = (float) (voxel_max[0] + 1) - 0.5F
            * (float) sites_x;
        cluster_p->minmax_y[1] = (float) (voxel_max[1] + 1) - 0.5F
            * (float) sites_y;
        cluster_p->minmax_z[1] = (float) (voxel_max[2] + 1) - 0.5F
            * (float) sites_z;
      }
      delete[] net->cluster_id;
    }

    rayTracer::rayTracer(Net *net)
    {
      // Init globals
      blocks_yz = blocks_y * blocks_z;
      block_size_f = float(block_size);
      block_size2 = block_size * block_size;
      block_size3 = block_size * block_size2;
      block_size_1 = block_size - 1;

      block_size_inv = 1.F / (float) block_size;

      rtBuildClusters(net);
    }

    void rayTracer::render(const float iLbmStressType)
    {
      // the volume rendering is performed here

      float screen_max[4];
      float screen_vtx[4];
      float p0[3], p1[3], p2[3];
      float ray_dx[3];
      float cluster_x[3];
      float dir[3];
      float par1[3], par2[3], par3[3];
      float subimage_vtx[4];
      float scale_vec[4];
      float t_near, t_far;
      float **block_flow_field;
      float temp1;

      int subimage_pix[4];
      bool ray_sign[3];

      int pixels_x, pixels_y;
      int cluster_id;
      // int viewpoint_flag;
      int i, j, k, l;

      AABB aabb;

      Cluster *cluster_p;

      ColPixel col_pixel;

      pixels_x = vis::controller->screen.pixels_x;
      pixels_y = vis::controller->screen.pixels_y;

      screen_max[0] = vis::controller->screen.max_x;
      screen_max[1] = vis::controller->screen.max_x;
      screen_max[2] = vis::controller->screen.max_y;
      screen_max[3] = vis::controller->screen.max_y;

      for (l = 0; l < 3; l++)
      {
        p0[l] = vis::controller->viewpoint.x[l];
      }
      for (l = 0; l < 3; l++)
      {
        par1[l] = vis::controller->screen.dir1[l];
        par2[l] = vis::controller->screen.dir2[l];
        screen_vtx[l] = vis::controller->screen.vtx[l];
      }
      scale_vec[0] = scale_vec[1] = vis::controller->screen.scale_x;
      scale_vec[2] = scale_vec[3] = vis::controller->screen.scale_y;

      for (cluster_id = 0; cluster_id < clusters; cluster_id++)
      {
        cluster_p = &cluster[cluster_id];

        // the image-based projection of the cluster bounding box is
        // calculated here

        for (l = 0; l < 3; l++)
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

        block_flow_field = cluster_flow_field[cluster_id];

        subimage_vtx[0] = 1.0e+30F;
        subimage_vtx[1] = -1.0e+30F;
        subimage_vtx[2] = 1.0e+30F;
        subimage_vtx[3] = -1.0e+30F;

        for (i = 0; i < 2; i++)
        {
          p1[0] = cluster_p->minmax_x[i];

          for (j = 0; j < 2; j++)
          {
            p1[1] = cluster_p->minmax_y[j];

            for (k = 0; k < 2; k++)
            {
              p1[2] = cluster_p->minmax_z[k];

              vis::controller->project(p1, p2);

              subimage_vtx[0] = fminf(subimage_vtx[0], p2[0]);
              subimage_vtx[1] = fmaxf(subimage_vtx[1], p2[0]);
              subimage_vtx[2] = fminf(subimage_vtx[2], p2[1]);
              subimage_vtx[3] = fmaxf(subimage_vtx[3], p2[1]);
            }
          }
        }
        subimage_pix[0] = (int) (scale_vec[0] * (subimage_vtx[0]
            + screen_max[0]) + 0.5F);
        subimage_pix[1] = (int) (scale_vec[1] * (subimage_vtx[1]
            + screen_max[1]) + 0.5F);
        subimage_pix[2] = (int) (scale_vec[2] * (subimage_vtx[2]
            + screen_max[2]) + 0.5F);
        subimage_pix[3] = (int) (scale_vec[3] * (subimage_vtx[3]
            + screen_max[3]) + 0.5F);

        if (subimage_pix[0] >= pixels_x || subimage_pix[1] < 0
            || subimage_pix[2] >= pixels_y || subimage_pix[3] < 0)
        {
          continue;
        }
        subimage_pix[0] = util::max(subimage_pix[0], 0);
        subimage_pix[1] = util::min(subimage_pix[1], pixels_x - 1);
        subimage_pix[2] = util::max(subimage_pix[2], 0);
        subimage_pix[3] = util::min(subimage_pix[3], pixels_y - 1);

        // if (p0[0] >= v[0][0] && p0[1] >= v[0][1] && p0[2] >= v[0][2] &&
        // 	  p0[0] <= v[1][0] && p0[1] <= v[1][1] && p0[2] <= v[1][2])
        // 	{
        // 	  viewpoint_flag = 1;
        // 	}
        // else
        // 	{
        // 	  viewpoint_flag = 0;
        // 	}
        aabb.acc_1 = cluster_p->minmax_x[1] - p0[0];
        aabb.acc_2 = cluster_p->minmax_x[0] - p0[0];
        aabb.acc_3 = cluster_p->minmax_y[1] - p0[1];
        aabb.acc_4 = cluster_p->minmax_y[0] - p0[1];
        aabb.acc_5 = cluster_p->minmax_z[1] - p0[2];
        aabb.acc_6 = cluster_p->minmax_z[0] - p0[2];

        for (l = 0; l < 3; l++)
        {
          par3[l] = screen_vtx[l] + subimage_pix[0] * par1[l] + subimage_pix[2]
              * par2[l];
        }
        for (i = subimage_pix[0]; i <= subimage_pix[1]; i++)
        {
          for (l = 0; l < 3; l++)
          {
            dir[l] = par3[l];
          }
          for (j = subimage_pix[2]; j <= subimage_pix[3]; j++)
          {
            ray_dir[0] = dir[0];
            ray_dir[1] = dir[1];
            ray_dir[2] = dir[2];

            temp1 = 1.0F / sqrtf(dir[0] * dir[0] + dir[1] * dir[1] + dir[2]
                * dir[2]);

            ray_dir[0] *= temp1;
            ray_dir[1] *= temp1;
            ray_dir[2] *= temp1;

            ray_inv[0] = 1.0F / ray_dir[0];
            ray_inv[1] = 1.0F / ray_dir[1];
            ray_inv[2] = 1.0F / ray_dir[2];

            ray_sign[0] = ray_dir[0] > 0.0F;
            ray_sign[1] = ray_dir[1] > 0.0F;
            ray_sign[2] = ray_dir[2] > 0.0F;

            dir[0] += par2[0];
            dir[1] += par2[1];
            dir[2] += par2[2];

            // t_near = 0.0F;

            // if (!viewpoint_flag)
            {
              rtAABBvsRayFn(&aabb, ray_inv[0], ray_inv[1], ray_inv[2], &t_near,
                            &t_far, ray_sign);

              if (t_near > t_far)
                continue;

              ray_dx[0] = t_near * ray_dir[0] - cluster_x[0];
              ray_dx[1] = t_near * ray_dir[1] - cluster_x[1];
              ray_dx[2] = t_near * ray_dir[2] - cluster_x[2];
            }
            // else
            // 	{
            // 	  ray_dx[0] = cluster_x[0];
            // 	  ray_dx[1] = cluster_x[1];
            // 	  ray_dx[2] = cluster_x[2];
            // 	}
            ray_vel_col[0] = 0.0F;
            ray_vel_col[1] = 0.0F;
            ray_vel_col[2] = 0.0F;

            if (iLbmStressType != SHEAR_STRESS)
            {
              ray_stress_col[0] = 0.0F;
              ray_stress_col[1] = 0.0F;
              ray_stress_col[2] = 0.0F;
            }
            ray_length = 0.0F;
            ray_t_min = 1.0e+30F;
            ray_density = -1.0F;

            rtTraverseBlocksFn(ray_dx, block_flow_field,
                               ColourPalette::pickColour, ray_sign,
                               iLbmStressType);

            if (ray_t_min >= 1.e+30F)
              continue;

            col_pixel.vel_r = ray_vel_col[0] * 255.0F;
            col_pixel.vel_g = ray_vel_col[1] * 255.0F;
            col_pixel.vel_b = ray_vel_col[2] * 255.0F;

            if (iLbmStressType != SHEAR_STRESS)
            {
              col_pixel.stress_r = ray_stress_col[0] * 255.0F;
              col_pixel.stress_g = ray_stress_col[1] * 255.0F;
              col_pixel.stress_b = ray_stress_col[2] * 255.0F;
            }
            col_pixel.dt = ray_length;
            col_pixel.t = ray_t_min + t_near;
            col_pixel.density = (ray_density
                - vis::controller->density_threshold_min)
                * vis::controller->density_threshold_minmax_inv;

            if (ray_stress < 1.0e+30F)
            {
              col_pixel.stress = ray_stress
                  * vis::controller->stress_threshold_max_inv;
            }
            else
            {
              col_pixel.stress = 1.0e+30F;
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

    void rtUpdateClusterVoxel(int i,
                              float density,
                              float velocity,
                              float stress)
    {
      *cluster_voxel[3 * i] = density;
      *cluster_voxel[3 * i + 1] = velocity;
      *cluster_voxel[3 * i + 2] = stress;
    }

    rayTracer::~rayTracer()
    {
      int m, n;

      for (n = 0; n < clusters; n++)
      {
        for (m = 0; m < (cluster[n].blocks_x * cluster[n].blocks_y
            * cluster[n].blocks_z); m++)
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

      cluster.clear();
    }

  }
}
