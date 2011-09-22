#ifndef HEMELB_VIS_CLUSTERTRAVERSER_H
#define HEMELB_VIS_CLUSTERTRAVERSER_H

#include "geometry/VolumeTraverser.h"
#include "vis/rayTracer/Cluster.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      //ClusterTraverser is used to traverse the cluster
      template <typename ClusterType>
	class ClusterTraverser : public geometry::VolumeTraverser
      {
      public:
        ClusterTraverser(const ClusterType& iCluster) :
	mCluster(iCluster)
	{
	}

	virtual ~ClusterTraverser()
	{
	}

	virtual site_t GetXCount()
	{
	  return mCluster.blocksX;
	}

	virtual site_t GetYCount()
	{
	  return mCluster.blocksY;
	}

	virtual site_t GetZCount()
	{
	  return mCluster.blocksZ;
	}
	
      private:
	const ClusterType& mCluster;
      }; 
    }
  }
}


#endif // HEMELB_VIS_CLUSTERTRAVERSER_H
