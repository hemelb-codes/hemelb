#ifndef HEMELB_VIS_CLUSTERTRAVERSER_H
#define HEMELB_VIS_CLUSTERTRAVERSER_H

#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/VolumeTraverser.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      //ClusterTraverser is used to traverse the cluster
      class ClusterTraverser : public VolumeTraverser
      {
      public:
	ClusterTraverser(const Cluster& iCluster);
	~ClusterTraverser();

	site_t CurrentClusterNumber();

	Vector3D<site_t> GetSiteCoordinatesOfLowestSiteInCurrentCluster();

	virtual site_t GetXCount();

	virtual site_t GetYCount();

	virtual site_t GetZCount();

	const Cluster& mCluster;
      }; 
    }
  }
}


#endif // HEMELB_VIS_CLUSTERTRAVERSER_H
