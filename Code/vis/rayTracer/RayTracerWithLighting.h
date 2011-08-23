#ifndef HEMELB_VIS_RAYTRACERWITHLIGHTING_H
#define HEMELB_VIS_RAYTRACERWITHLIGHTING_H

#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class RayTracerWithLighting : public RayTracer
      {
      public:
	virtual void BuildClusters();
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACERWITHLIGHTING_H
