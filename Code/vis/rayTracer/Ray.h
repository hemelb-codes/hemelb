#ifndef HEMELB_VIS_RAY_H
#define HEMELB_VIS_RAY_H

#include "vis/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class Ray
      {
      public:
	Ray();
	Ray(Vector3D<float> iDirection);

	Vector3D<float> GetDirection() const;
	Vector3D<float> GetInverseDirection() const;

	bool XIncreasing() const;
	bool YIncreasing() const;
	bool ZIncreasing() const;
	
	float Length;

	float VelocityColour[3];
	float StressColour[3];
	float Stress;
	float Density;
	float LengthToFirstRayIntersection;
	
      private:
	Vector3D<float> mDirection;
	Vector3D<float> mInverseDirection;

      };	
    }
  }
}

#endif // HEMELB_VIS_RAY_H
