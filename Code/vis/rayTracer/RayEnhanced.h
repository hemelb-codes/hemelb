#ifndef HEMELB_VIS_RAYENHANCED_H
#define HEMELB_VIS_RAYENHANCED_H

#include "vis/Vector3D.h"
#include "vis/rayTracer/Ray.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class RayEnhanced : public Ray
      {
      public:
	RayEnhanced();
	RayEnhanced(Vector3D<float> iDirection);

	float GetIntensity();

	void HandleWallIntersection(Vector3D<float>& iWallNormal);

      private:
 	 float mIntensity;

	 static const float mMinLogIntensityMultiple;
      };	
    }
  }
}

#endif // HEMELB_VIS_RAYENHANCED_H
