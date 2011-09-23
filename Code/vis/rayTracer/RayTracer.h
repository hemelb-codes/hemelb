#ifndef HEMELB_VIS_RAYTRACER_H
#define HEMELB_VIS_RAYTRACER_H

//#define NDEBUG;
#include <assert.h>

#include <map>
#include <stack>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <limits>

#include "constants.h"
#include "debug/Debugger.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"
#include "log/Logger.h"
#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h" 
#include "util/Vector3D.h"
#include "vis/DomainStats.h"
#include "vis/Screen.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"
#include "vis/XYCoordinates.h"
#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/ClusterBuilder.h"
#include "vis/rayTracer/ClusterRayTracer.h"
#include "vis/rayTracer/Ray.h"
#include "vis/rayTracer/RayTracer.h"
#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      template<typename ClusterType, typename RayDataType>
      class RayTracer
      {
        public:
          // Constructor and destructor do all the usual stuff.
          RayTracer(const geometry::LatticeData* iLatDat,
                    const DomainStats* iDomainStats,
                    Screen* iScreen,
                    Viewpoint* iViewpoint,
                    VisSettings* iVisSettings) :
              mClusterBuilder(iLatDat), mLatDat(iLatDat), mDomainStats(iDomainStats), mScreen(iScreen), mViewpoint(iViewpoint), mVisSettings(iVisSettings)
          {
          }

          ~RayTracer()
          {
          }

          //Calls the cluster builder to build the clusters
          void BuildClusters()
          {
            mClusterBuilder.BuildClusters();
          }

          // Method to update the voxel corresponding to site i with its
          // newly calculated density, velocity and stress.
          void UpdateClusterVoxel(site_t i,
                                  distribn_t density,
                                  distribn_t velocity,
                                  distribn_t stress)
          {
            assert(static_cast<site_t>(static_cast<unsigned int>(i)) == i);

            mClusterBuilder.GetClusterVoxelDataPointer(i)->SetDensity(static_cast<float>(density));
            mClusterBuilder.GetClusterVoxelDataPointer(i)->SetVelocity(static_cast<float>(velocity));
            mClusterBuilder.GetClusterVoxelDataPointer(i)->SetStress(static_cast<float>(stress));
          }

          // Render the current state into an image.
          void Render()
          {
            ClusterRayTracer<ClusterType, RayDataType> lClusterRayTracer(*mViewpoint,
                                                                         *mScreen,
                                                                         *mDomainStats,
                                                                         *mVisSettings,
                                                                         *mLatDat);

            for (unsigned int clusterId = 0; clusterId < mClusterBuilder.GetClusters().size();
                clusterId++)
                {
              lClusterRayTracer.RenderCluster(mClusterBuilder.GetClusters()[clusterId]);
            }
          }

        private:
          ClusterBuilder<ClusterType> mClusterBuilder;
          const geometry::LatticeData* mLatDat;

          const DomainStats* mDomainStats;
          Screen* mScreen;
          Viewpoint* mViewpoint;
          VisSettings* mVisSettings;
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_H
