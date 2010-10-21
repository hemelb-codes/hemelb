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

    // Some sort of coordinates.
    struct BlockLocation
    {
        short int i, j, k;
    };

    struct Cluster
    {
        float minmax_x[2], minmax_y[2], minmax_z[2];

        float x[3];

        unsigned short int blocks_x, blocks_y, blocks_z;
        unsigned short int block_min[3];
    };

    struct AABB
    {
        float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
    };

    void rtUpdateClusterVoxel(int i,
                              float density,
                              float velocity,
                              float stress);

    class rayTracer
    {
      public:
        rayTracer(Net *net);
        ~rayTracer();

        void rtUpdateRayData(float *flow_field,
                             float ray_t,
                             float ray_segment,
                             void(*ColourPalette)(float value, float col[]),
                             const float iLbmStressType);
        void render(const float iLbmStressType);

      private:
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
    };

  }
}

#endif // HEMELB_VIS_RAYTRACER_H
