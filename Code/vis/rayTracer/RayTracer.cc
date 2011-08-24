//#define NDEBUG;
#include <assert.h>


#include <math.h>
#include <stdlib.h>
#include <vector>
#include <limits>

#include "debug/Debugger.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h" 
#include "vis/Vector3D.h"
#include "vis/XYCoordinates.h"
#include "vis/rayTracer/ClusterRayTracer.h"
#include "vis/rayTracer/RayTracer.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      // TODO RENAME THIS FUNCTION
            void RayTracer::BuildClusters()
      {
	mClusterBuilder = new ClusterBuilder(mLatDat);
	mClusterBuilder->BuildClusters();
      }

 

            
      RayTracer::RayTracer(const geometry::LatticeData* iLatDat,
			   const DomainStats* iDomainStats,
			   Screen* iScreen,
			   Viewpoint* iViewpoint,
			   VisSettings* iVisSettings) :
  	mLatDat(iLatDat), mDomainStats(iDomainStats), mScreen(iScreen), mViewpoint(iViewpoint),
	mVisSettings(iVisSettings), mBlockSizeFloat((float) mLatDat->GetBlockSize()),
	mBlockSizeInverse(1.F / mBlockSizeFloat), block_size2(mLatDat->GetBlockSize()
							      * mLatDat->GetBlockSize()), block_size3(mLatDat->GetBlockSize() * block_size2),
	block_size_1(mLatDat->GetBlockSize() - 1), blocksYz(mLatDat->GetYBlockCount()
							    * mLatDat->GetZBlockCount())
      {
	mClusterBuilder = NULL;
      }

      void RayTracer::Render()
      {
	ClusterRayTracer lClusterRayTracer(*mViewpoint, *mScreen, *mDomainStats, *mVisSettings, *mLatDat);

	for (unsigned int clusterId = 0; clusterId < mClusterBuilder->GetClusters().size(); clusterId++)
	{
	  lClusterRayTracer.RenderCluster(*mClusterBuilder->GetClusters()[clusterId]);
	}
      }      

      void RayTracer::UpdateClusterVoxel(site_t i,
					 distribn_t density,
					 distribn_t velocity,
					 distribn_t stress)
      {
	assert(static_cast<site_t>(static_cast<unsigned int>(i)) == i);
	
	
	mClusterBuilder->GetClusterVoxelDataPointer(i)->Density =
	  (float) density;
	mClusterBuilder->GetClusterVoxelDataPointer(i)->Velocity =
	  (float) velocity;
	mClusterBuilder->GetClusterVoxelDataPointer(i)->Stress =
	  (float) stress;
      }

      RayTracer::~RayTracer()
      {
	delete mClusterBuilder;
      }
    }
  }
}
