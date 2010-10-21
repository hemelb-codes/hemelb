#ifndef HEMELB_VIS_RAYTRACER_H
#define HEMELB_VIS_RAYTRACER_H

#include "constants.h"
#include "net.h"

// TODO this could probably be reduced to the net class and some visualisation class.

namespace hemelb
{
  namespace vis
  {
    // Horrifying globals
    extern float block_size_f;
    extern float block_size_inv;
    extern int block_size2, block_size3, block_size_1;

    class rayTracer
    {
      public:
        rayTracer(Net *net);
        ~rayTracer();

        void rtUpdateClusterVoxel(int i,
                                  float density,
                                  float velocity,
                                  float stress);

        void render(const float iLbmStressType);

      private:

        struct AABB
        {
            float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
        };

        struct Cluster
        {
            float minmax_x[2], minmax_y[2], minmax_z[2];

            float x[3];

            unsigned short int blocks_x, blocks_y, blocks_z;
            unsigned short int block_min[3];
        };

        // Some sort of coordinates.
        struct BlockLocation
        {
            short int i, j, k;
        };

        void rtUpdateRayData(float *flow_field,
                             float ray_t,
                             float ray_segment,
                             void(*ColourPalette)(float value, float col[]),
                             const float iLbmStressType);

        void rtTraverseVoxels(float block_min[],
                              float block_x[],
                              float voxel_flow_field[],
                              float t,
                              void(*ColourPalette)(float value, float col[]),
                              bool xyz_is_1[],
                              const float iLbmStressType);

        void rtTraverseBlocksFn(float ray_dx[],
                                float **block_flow_field,
                                void(*ColourPalette)(float value, float col[]),
                                bool xyz_Is_1[],
                                const double iLbmStressType);

        void rtAABBvsRayFn(AABB *aabb,
                           float inv_x,
                           float inv_y,
                           float inv_z,
                           float *t_near,
                           float *t_far,
                           bool xyz_sign_is_1[]);

        void rtUpdateColour(float dt, float palette[], float col[]);

        void rtBuildClusters(Net *net);

        std::vector<Cluster> cluster;
        float **cluster_voxel;
        float ***cluster_flow_field;

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
    };

  }
}

#endif // HEMELB_VIS_RAYTRACER_H
