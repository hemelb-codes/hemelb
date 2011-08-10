#include <math.h>


#include "vis/Location.h"
#include "vis/Viewpoint.h"

namespace hemelb
{
  namespace vis
  {
    Viewpoint::Viewpoint()
      //mViewpointCentre(0.0F) 
    {}
    
    Location<float> Viewpoint::RotateToViewpoint(const Location<float>& iVector) const
    {
      // A rotation of iThetaX clockwise about the x-axis
      // Followed by a rotation of iThetaY anticlockwise about the y-axis.

      return Rotate(SinXRotation, CosXRotation, SinYRotation, CosYRotation, iVector);
    }

    Location <float> Viewpoint::Rotate
    (float iSinX,
     float iCosX,
     float iSinY,
     float iCosY,
     const Location<float>& iVector)
    {
      // A rotation of iThetaX clockwise about the x-axis
      // Followed by a rotation of iThetaY anticlockwise about the y-axis.
      // In matrices:
      //       (cos(iThetaY)  0 sin(iThetaY)) (1 0            0              )
      // Out = (0             1 0           ) (0 cos(-iThetaX) -sin(-iThetaX)) In
      //       (-sin(iThetaY) 0 cos(iThetaY)) (0 sin(-iThetaX) cos(-iThetaX) )
      //
      //       (Xcos(iThetaY) + Zsin(iThetaY)cos(iThetaX) - Ysin(iThetaY)sin(iThetaX))
      // Out = (Ycos(iThetaX) + Zsin(iThetaX)                                        )
      //       (Zcos(iThetaX)cos(iThetaY) - Ysin(iThetaX)cos(iThetaY) - Xsin(iThetaY))

      const float lTemp = iVector.z * iCosX - iVector.y * iSinX;
      
      return Location <float>(
	lTemp*iSinY + iVector.x*iCosY,
	iVector.z * iSinX + iVector.y * iCosX,
	lTemp * iCosY - iVector.x * iSinY);
    }

    Location<float> Viewpoint::Project(const Location<float>& p1) const
    {
      Location<float> x1;
      

      x1 = p1 - mViewpointCentre;
      float temp1 = x1.x;
      x1.x = x1.y;
      x1.y = temp1;
      
      Location<float> x2 = Rotate( -SinYRotation, CosYRotation, 
				  -SinXRotation, CosXRotation, x1);
      
      

      float temp2 = mDistance / (-x2.z);
     
      return Location <float> (temp2 * x2.y,
			       temp2 * x2.x,
			       -x2.z);

    }

    /**
     * Set the position of the viewpoint.
     *
     * @param longitude in radians.
     * @param latitude in radians.
     * @param localCentre
     * @param distance
     */
    void Viewpoint::SetViewpointPosition(
      float longitude,
      float latitude,
      const Location<float>& iLocalCentre,
      float rad,
      float iDistance)
    {
      SinYRotation = sinf(longitude);
      CosYRotation = cosf(longitude);

      SinXRotation = sinf(latitude);
      CosXRotation = cosf(latitude);

      mViewpointCentre = RotateToViewpoint(Location<float>(0., 0., rad));

      mViewpointCentre += iLocalCentre;

      mDistance = iDistance;
    }

    const Location<float>& Viewpoint::GetViewpointCentre() const
    {
      return mViewpointCentre;
    }

  }
}
