#include <limits>

#include "vis/rayTracer/Ray.h"
#include "vis/Vector3D.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      Ray::Ray(Vector3D<float> iDirection) :
	mDirection(iDirection)
	//Caution: ensure normalised 
      {
	mInverseDirection = Vector3D<float>
	  ( 1.0F/iDirection.x,
	    1.0F/iDirection.y,
	    1.0F/iDirection.z);

	VelocityColour[0] = 0.0F;
	VelocityColour[1] = 0.0F;
	VelocityColour[2] = 0.0F;

	StressColour[0] = 0.0F;
	StressColour[1] = 0.0F;
	StressColour[2] = 0.0F;

	Length = 0.0F;
	LengthToFirstRayIntersection = std::numeric_limits<float>::max();
	Density = -1.0F;
      }
      
      Vector3D<float> Ray::GetDirection() const
      {
	return mDirection;
      }

      Vector3D<float> Ray::GetInverseDirection() const
      {
	return mInverseDirection;
      }

      bool Ray::XIncreasing() const
      {
	return GetDirection().x > 0.0F;
      }

      bool Ray::YIncreasing() const
      {
	return GetDirection().y > 0.0F;
      }

      bool Ray::ZIncreasing() const
      {
	return GetDirection().z > 0.0F;
      }

    }
  }
}
