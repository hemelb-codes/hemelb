#include <math.h>

#include "log/Logger.h"
#include "vis/Vector3D.h"
#include "vis/Viewpoint.h"
#include "vis/XYCoordinates.h"

namespace hemelb
{
  namespace vis
  {
    Viewpoint::Viewpoint() :
      mViewpointLocation(0.0F) 
    {}
    
    Vector3D<float> Viewpoint::RotateCameraCoordinatesToWorldCoordinates(const Vector3D<float>& iVector) const
    {
      // A rotation of iThetaX clockwise looking up the x-axis (increasing)
      // Followed by a rotation of iThetaY anticlockwise looking up the y-axis.
      return Rotate(mSinLatitude, mCosLatitude, mSinLongitude, mCosLongitude, iVector);
    }

    Vector3D<float> Viewpoint::RotateWorldToCameraCoordinates(const Vector3D<float>& iVector) const
    {
      // The reverse of the above
      return UnRotate(mSinLatitude, mCosLatitude, mSinLongitude, mCosLongitude, iVector);
    }

    Vector3D <float> Viewpoint::Rotate
    (float iSinThetaX,
     float iCosThetaX,
     float iSinThetaY,
     float iCosThetaY,
     const Vector3D<float>& iVector)
    {
      // A rotation of iThetaX clockwise looking down the x-axis
      // Followed by a rotation of iThetaY anticlockwise looking down the y-axis.
      // In matrices:
      //       (cos(iThetaY)  0 sin(iThetaY)) (1 0            0              )
      // Out = (0             1 0           ) (0 cos(-iThetaX) -sin(-iThetaX)) In
      //       (-sin(iThetaY) 0 cos(iThetaY)) (0 sin(-iThetaX) cos(-iThetaX) )
      //
      //       (Xcos(iThetaY) + Zsin(iThetaY)cos(iThetaX) - Ysin(iThetaY)sin(iThetaX))
      // Out = (Ycos(iThetaX) + Zsin(iThetaX)                                        )
      //       (Zcos(iThetaX)cos(iThetaY) - Ysin(iThetaX)cos(iThetaY) - Xsin(iThetaY))

      const float lTemp = iVector.z * iCosThetaX - iVector.y * iSinThetaX;
      
      return Vector3D <float>(
	lTemp*iSinThetaY + iVector.x*iCosThetaY,
	iVector.z * iSinThetaX + iVector.y * iCosThetaX,
	lTemp * iCosThetaY - iVector.x * iSinThetaY);
    }

    Vector3D <float> Viewpoint::UnRotate
    (float iSinThetaX,
     float iCosThetaX,
     float iSinThetaY,
     float iCosThetaY,
     const Vector3D<float>& iVector)
    {
      // A rotation of iThetaY aniclockwise looking down the y-axis
      // Followed by a rotation of iThetaX clockwise looking down the x-axis.
      // In matrices:
      //       (1 0             0            )(cos(-iThetaY) 0 sin(-iThetaY)) 
      // Out = (0 cos(iThetaX)  -sin(iThetaX))(0             1 0           ) In
      //       (0 sin(iThetaX) cos(iThetaX) )(-sin(-iThetaY) 0 cos(-iThetaY)) 
      //
      // This is the Rotation matrix inversed / transposted
      // ie (AB)^-1 = (AB)^t = B^t A^T

      const float lTemp = iVector.x * iSinThetaY + iVector.z * iCosThetaY;
      
      return Vector3D <float>(
	iVector.x*iCosThetaY - iVector.z*iSinThetaY,
	-lTemp*iSinThetaX + iVector.y*iCosThetaX,
	lTemp*iCosThetaX + iVector.y*iSinThetaX);
    }

    Vector3D<float> Viewpoint::Project(const Vector3D<float>& iWorldLocation) const
    {
      Vector3D<float> lLocationCamCoordinates = GetViewPointLocationInCameraCoordinates(iWorldLocation);

      // NB - the oringinal code was not doing this but a reflection rotation
      // and then back reflectiond which produced similar images.
      // It was believed that this was in error. 
      
      //Carry out a perspective projection on an infinite spanning screen 
      //between the eye and the subject.
      //Reverse the sign such that depth is positive (I believe).  
      return Vector3D <float> ( mDistanceFromEyeToScreen / (-lLocationCamCoordinates.z)
				* lLocationCamCoordinates.x,
				 mDistanceFromEyeToScreen / (-lLocationCamCoordinates.z) 
				* lLocationCamCoordinates.y,
				-lLocationCamCoordinates.z);
    }

    XYCoordinates<float> Viewpoint::FlatProject(const Vector3D<float>& iWorldLocation) const
    {
      Vector3D<float> lLocationCamCoordinates = GetViewPointLocationInCameraCoordinates(iWorldLocation);
      
      return XYCoordinates<float> ( mDistanceFromEyeToScreen / (-lLocationCamCoordinates.z)
				* lLocationCamCoordinates.x,
				 mDistanceFromEyeToScreen / (-lLocationCamCoordinates.z) 
				* lLocationCamCoordinates.y);
    }


    void Viewpoint::SetViewpointPosition(
      float iLongitude,
      float iLatitude,
      const Vector3D<float>& iLocalCentre,
      float iRadius,
      float iDistanceFromEyeToScreen)
    {
      hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>
       	("Latitude: %f / Longitude: %f", iLatitude, iLongitude);
      mSinLongitude = sinf(iLongitude);
      mCosLongitude = cosf(iLongitude);
    
      mSinLatitude = sinf(iLatitude);
      mCosLatitude = cosf(iLatitude);

      // Z in the camera co-ordinates indicates the distance from the centre
      // the viewpoint into the camera prior rotation.
      mViewpointLocation = RotateCameraCoordinatesToWorldCoordinates(Vector3D<float>(0., 0., iRadius));

      //Translate the camera location to allow it to point at 
      //a local centre rather than the world centre
      mViewpointLocation += iLocalCentre;

      mDistanceFromEyeToScreen = iDistanceFromEyeToScreen;
    }

    const Vector3D<float>& Viewpoint::GetViewpointCentreLocation() const
    {
      return mViewpointLocation;
    }

    Vector3D<float> Viewpoint::GetViewPointLocationInCameraCoordinates(const Vector3D<float>& iWorldLocation) const
    {
      //Calculate the location of the point relative to the viewpoint 
      Vector3D<float> lLocationRelativeToViewPoint = iWorldLocation - mViewpointLocation;
      
      //Rotate the location vector in the opposite manner to that of the camera
      return RotateWorldToCameraCoordinates(lLocationRelativeToViewPoint);
    }

  }
}
