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
      }
      
      Vector3D<float> Ray::GetDirection() const
      {
	return mDirection;
      }

      Vector3D<float> Ray::GetInverseDirection() const
      {
	return mInverseDirection;
      }

    }
  }
}
