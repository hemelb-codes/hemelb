#ifndef HEMELB_VIS_RAYTRACERENHANCED_H
#define HEMELB_VIS_RAYTRACERENHANCED_H

#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class RayTracerEnhanced : public RayTracer
      {
      public:
	RayTracerEnhanced(const geometry::LatticeData* iLatDat,
			  const DomainStats* iDomainStats,
			  Screen* iScreen,
			  Viewpoint* iViewpoint,
			  VisSettings* iVisSettings);

	virtual void BuildClusters();

	virtual void Render();
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACERENHANCED_H
