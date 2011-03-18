#ifndef HEMELB_VIS_RAYTRACER_H
#define HEMELB_VIS_RAYTRACER_H

#include "constants.h"
#include "geometry/LocalLatticeData.h"
#include "geometry/GlobalLatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/DomainStats.h"
#include "vis/Screen.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    class RayTracer
    {
      public:
        // Constructor and destructor do all the usual stuff.
        RayTracer(const topology::NetworkTopology * iNetworkTopology,
                  const geometry::LocalLatticeData* iLocalLatDat,
                  const geometry::GlobalLatticeData* iGlobLatDat,
                  DomainStats* iDomainStats,
                  Screen* iScreen,
                  Viewpoint* iViewpoint,
                  VisSettings* iVisSettings);
        ~RayTracer();

        // Method to update the voxel corresponding to site i with its
        // newly calculated density, velocity and stress.
        void UpdateClusterVoxel(const int &i,
                                const float &density,
                                const float &velocity,
                                const float &stress);

        // Render the current state into an image.
        void Render(const lb::StressTypes iLbmStressType);

      private:

        struct Ray
        {
            float Direction[3];
            float InverseDirection[3];
            float Length;

            float VelocityColour[3];
            float StressColour[3];
            float Stress;
            float Density;
            float MinT;
        };

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
                             Ray *bCurrentRay,
                             void(*ColourPalette)(float value, float col[]),
                             const lb::StressTypes iLbmStressType);

        void rtTraverseVoxels(float block_min[],
                              float block_x[],
                              float voxel_flow_field[],
                              float t,
                              Ray *bCurrentRay,
                              void(*ColourPalette)(float value, float col[]),
                              bool xyz_is_1[],
                              const lb::StressTypes iLbmStressType);

        void rtTraverseBlocksFn(float ray_dx[],
                                float **block_flow_field,
                                Ray *bCurrentRay,
                                void(*ColourPalette)(float value, float col[]),
                                bool xyz_Is_1[],
                                const lb::StressTypes iLbmStressType);

        void rtAABBvsRayFn(const AABB &aabb,
                           const float &inv_x,
                           const float &inv_y,
                           const float &inv_z,
                           const bool xyz_sign_is_1[],
                           float &t_near,
                           float &t_far);

        void rtUpdateColour(float dt, float palette[], float col[]);

        void rtBuildClusters();

        const topology::NetworkTopology * mNetworkTopology;
        const geometry::LocalLatticeData* mLocalLatDat;
        const geometry::GlobalLatticeData* mGlobLatDat;

        DomainStats* mDomainStats;
        Screen* mScreen;
        Viewpoint* mViewpoint;
        VisSettings* mVisSettings;

        std::vector<Cluster*> mClusters;
        float **cluster_voxel;
        float ***cluster_flow_field;

        int cluster_blocks_vec[3];
        int cluster_blocks_z, cluster_blocks_yz, cluster_blocks;

        float mBlockSizeFloat;
        float mBlockSizeInverse;
        unsigned int block_size2, block_size3, block_size_1;
        unsigned int blocks_yz;
    };

  }
}

#endif // HEMELB_VIS_RAYTRACER_H
