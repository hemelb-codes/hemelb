// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_VIEWPOINT_H
#define HEMELB_VIS_VIEWPOINT_H

#include "util/Vector3D.h"
#include "vis/XYCoordinates.h"

namespace hemelb
{
  namespace vis
  {
    /**
     * Holds the position of the viewpoint or camera location and performs 
     * projection of points in world coordinates into 3D screen co-ordinates
     */
    class Viewpoint
    {
        //See http://en.wikipedia.org/wiki/3D_projection#Perspective_projection
      public:
        Viewpoint();

        /**
         * Changes a location in camera coordinates into world coordinates
         * Note that both co-ordinate systems share the same (0,0,0) point
         * and the camera is at (0,0,radius) in camera coordinates
         * 
         */
        util::Vector3D<float> RotateCameraCoordinatesToWorldCoordinates(const util::Vector3D<float>& iVector) const;

        /**
         * Does the reverse of the above
         *
         */
        util::Vector3D<float> RotateWorldToCameraCoordinates(const util::Vector3D<float>& iVector) const;

        /**
         * Projects a location in world coordinates onto the infinite screen,
         * by translating and rotating such as to give coordinates relative
         * to the camera, and then applying a perspective projection
         */
        util::Vector3D<float> Project(const util::Vector3D<float>& p1) const;

        /**
         * Same as project but doesn't return a z component
         *
         */
        XYCoordinates<float> FlatProject(const util::Vector3D<float>& iWorldLocation) const;

        /**
         * Sets the position of the Camera or Viewpoint
         *
         * In world coordinates, the camera is located at radius from
         * the local centre, and pointing at an angle indicated by the
         * latitude and longitude.
         *
         * @param iLongitude - longitude in radians.
         * @param iLatitude - latitude in radians.
         * @param iLocalCentre - where the camera should point in world-coordinates
         * @param iRadius - the distance of the camera from this local centre
         * @param iDistanceFromCameraToScreen - the distance of the infinite screen
         * from the viewer. Allows for zoom
         */
        void SetViewpointPosition(float iLongitude,
                                  float iLatitude,
                                  const util::Vector3D<float>& iLocalCentre,
                                  float iRadius,
                                  float iDistanceFromCameraToScreen);

        const util::Vector3D<float>& GetViewpointLocation() const;

        float GetDistanceFromCameraToScreen() const;

      private:
        util::Vector3D<float> GetLocationInCameraCoordinates(const util::Vector3D<float>& iWorldLocation) const;

        //Performs a vector rotation using stored
        //Sin and Cosine Values
        static util::Vector3D<float> Rotate(float iSinThetaX,
                                            float iCosThetaX,
                                            float iSinThetaY,
                                            float iCosThetaY,
                                            const util::Vector3D<float>& iVector);

        //Reverses a vector rotation of the above
        static util::Vector3D<float> UnRotate(float iSinThetaX,
                                              float iCosThetaX,
                                              float iSinThetaY,
                                              float iCosThetaY,
                                              const util::Vector3D<float>& iVector);

        float mSinLongitude;
        float mCosLongitude;
        float mSinLatitude;
        float mCosLatitude;

        //Stores the viewpoint Location in world co-ordinate
        util::Vector3D<float> mViewpointLocationInWorldCoordinates;

        float mDistanceFromCameraToScreen;
    };
  }
}

#endif /* HEMELB_VIS_VIEWPOINT_H */
