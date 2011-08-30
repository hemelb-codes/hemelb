#include <limits>

#include "vis/rayTracer/RayEnhanced.h"
#include "vis/Vector3D.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {

      const float RayEnhanced::mMinLogIntensityMultiple = 0.5F;

      RayEnhanced::RayEnhanced(Vector3D<float> iDirection) :
	Ray(iDirection)
	//Caution: ensure normalised 
      {
	mIntensity = 0.0F;
      }

      float RayEnhanced::GetIntensity()
      {
	return mIntensity;
      }
      
      void RayEnhanced::HandleWallIntersection(Vector3D<float>& iWallNormal)
      {
	//Calculate the absolute dot product of the wall
	//vector normal and the ray direction
	float lDotProduct = 
	  GetDirection().DotProduct(iWallNormal);

	//We want the attentuation to be between 
	//mMinLogIntensityMultiple and 1
	float lIntensityMultiple = mMinLogIntensityMultiple + 
	  (1.0F-mMinLogIntensityMultiple)*lDotProduct;

	//Scale the current intensity of the ray
	mIntensity *= lIntensityMultiple;
      }
    }
  }
}
