// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <cmath>

#include "log/Logger.h"
#include "util/Vector3D.h"
#include "vis/Viewpoint.h"
#include "vis/XYCoordinates.h"

namespace hemelb
{
  namespace vis
  {
    Viewpoint::Viewpoint() :
      mViewpointLocationInWorldCoordinates(0.0F)
    {
    }

    util::Vector3D<float> Viewpoint::RotateCameraCoordinatesToWorldCoordinates(const util::Vector3D<
        float>& iVector) const
    {
      // A rotation of iThetaX clockwise looking up the x-axis (increasing)
      // Followed by a rotation of iThetaY anticlockwise looking up the y-axis.
      return Rotate(mSinLatitude, mCosLatitude, mSinLongitude, mCosLongitude, iVector);
    }

    util::Vector3D<float> Viewpoint::RotateWorldToCameraCoordinates(const util::Vector3D<float>& iVector) const
    {
      // The reverse of the above
      return UnRotate(mSinLatitude, mCosLatitude, mSinLongitude, mCosLongitude, iVector);
    }

    util::Vector3D<float> Viewpoint::Rotate(float iSinThetaX,
                                            float iCosThetaX,
                                            float iSinThetaY,
                                            float iCosThetaY,
                                            const util::Vector3D<float>& iVector)
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

      return util::Vector3D<float>(lTemp * iSinThetaY + iVector.x * iCosThetaY,
                                   iVector.z * iSinThetaX + iVector.y * iCosThetaX,
                                   lTemp * iCosThetaY - iVector.x * iSinThetaY);
    }

    util::Vector3D<float> Viewpoint::UnRotate(float iSinThetaX,
                                              float iCosThetaX,
                                              float iSinThetaY,
                                              float iCosThetaY,
                                              const util::Vector3D<float>& iVector)
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

      return util::Vector3D<distribn_t>(iVector.x * iCosThetaY - iVector.z * iSinThetaY,
                                        -lTemp * iSinThetaX + iVector.y * iCosThetaX,
                                        lTemp * iCosThetaX + iVector.y * iSinThetaX);
    }

    util::Vector3D<float> Viewpoint::Project(const util::Vector3D<float>& iWorldLocation) const
    {
      util::Vector3D<float> lLocationCamCoordinates =
          GetLocationInCameraCoordinates(iWorldLocation);

      //Carry out a perspective projection on an infinite spanning screen 
      //between the camera and the subject.
      //Reverse the sign such that depth is positive (I believe).  
      return util::Vector3D<float>(mDistanceFromCameraToScreen / (-lLocationCamCoordinates.z)
                                       * lLocationCamCoordinates.x,
                                   mDistanceFromCameraToScreen / (-lLocationCamCoordinates.z)
                                       * lLocationCamCoordinates.y,
                                   -lLocationCamCoordinates.z);
    }

    XYCoordinates<float> Viewpoint::FlatProject(const util::Vector3D<float>& iWorldLocation) const
    {
      util::Vector3D<float> lLocationCamCoordinates =
          GetLocationInCameraCoordinates(iWorldLocation);

      return XYCoordinates<float> (mDistanceFromCameraToScreen / (-lLocationCamCoordinates.z)
                                       * lLocationCamCoordinates.x,
                                   mDistanceFromCameraToScreen / (-lLocationCamCoordinates.z)
                                       * lLocationCamCoordinates.y);
    }

    void Viewpoint::SetViewpointPosition(float iLongitude,
                                         float iLatitude,
                                         const util::Vector3D<float>& iLocalCentre,
                                         float iRadius,
                                         float iDistanceFromCameraToScreen)
    {

      mSinLongitude = sinf(iLongitude);
      mCosLongitude = cosf(iLongitude);

      mSinLatitude = sinf(iLatitude);
      mCosLatitude = cosf(iLatitude);

      //The camera is located at 0,0,radius from the world centre in camera co-ordinates 
      mViewpointLocationInWorldCoordinates
          = RotateCameraCoordinatesToWorldCoordinates(util::Vector3D<float>(0., 0., iRadius));

      //Translate the camera location to allow it to point at 
      //a local centre rather than the world centre
      mViewpointLocationInWorldCoordinates += iLocalCentre;

      mDistanceFromCameraToScreen = iDistanceFromCameraToScreen;
    }

    const util::Vector3D<float>& Viewpoint::GetViewpointLocation() const
    {
      return mViewpointLocationInWorldCoordinates;
    }

    float Viewpoint::GetDistanceFromCameraToScreen() const
    {
      return mDistanceFromCameraToScreen;
    }

    util::Vector3D<float> Viewpoint::GetLocationInCameraCoordinates(const util::Vector3D<float>& iWorldLocation) const
    {
      //Calculate the location of the point relative to the viewpoint 
      util::Vector3D<float> lLocationRelativeToViewPoint = iWorldLocation
          - mViewpointLocationInWorldCoordinates;

      //Rotate the location vector in the opposite manner to that of the camera
      return RotateWorldToCameraCoordinates(lLocationRelativeToViewPoint);
    }

  }
}
