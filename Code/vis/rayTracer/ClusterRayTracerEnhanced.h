#ifndef HEMELB_VIS_CLUSTERRAYTRACERENHANCED_H
#define HEMELB_VIS_CLUSTERRAYTRACERENHANCED_H

#include "ClusterRayTracer.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class ClusterRayTracerEnhanced : public ClusterRayTracer
      {
      public:
	ClusterRayTracerEnhanced(const Viewpoint& iViewpoint,
				 Screen& iScreen,
				 const DomainStats& iDomainStats,
				 const VisSettings& iVisSettings,
				 const hemelb::geometry::LatticeData& iLatticeData);
   
      private:
	virtual void CastRayForPixel(const Cluster& iCluster,
			     const XYCoordinates<int>& iPixel,
			     const Vector3D<float>& iRayDirection);

	virtual void UpdateRayData
	  (const Cluster& iCluster,
	   site_t iBlockNumber,
	   site_t iSiteNumber,
	   float iLengthFromClusterFirstIntersectionToVoxel,
	   float iRayLengthInVoxel,
	   Ray& ioRay);

      };
    }
  }
}

#endif // HEMELB_VIS_CLUSTERRAYTRACERENHANCED_H
