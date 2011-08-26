//#define NDEBUG;
#include <assert.h>

#include "vis/Vector3D.h"
#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/ClusterTraverser.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterTraverser::ClusterTraverser
      (const Cluster& iCluster)
	: VolumeTraverser(),
	  mCluster(iCluster)
      {
      }

      ClusterTraverser::~ClusterTraverser()
      {
      }
     	  
      site_t ClusterTraverser::GetXCount() 
      {
	return mCluster.blocksX;
      }

      site_t ClusterTraverser::GetYCount()
      {
	return mCluster.blocksY;
      }

      site_t ClusterTraverser::GetZCount() 
      {
	return mCluster.blocksZ;
      }

    }
  }
}
