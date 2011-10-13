#include <stdlib.h>

#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include "vis/Screen.h"

namespace hemelb
{
  namespace vis
  {

    Screen::Screen()
    {
    }

    Screen::~Screen()
    {
    }

    void Screen::Set(float maxX,
                     float maxY,
                     int pixelsX,
                     int pixelsY,
                     float rad,
                     const Viewpoint* viewpoint)
    {
      MaxXValue = maxX;
      MaxYValue = maxY;

      UnitVectorProjectionX
          = viewpoint-> RotateCameraCoordinatesToWorldCoordinates(Vector3D<float> (MaxXValue,
                                                                                   0.0F,
                                                                                   0.0F));
      UnitVectorProjectionY
          = viewpoint-> RotateCameraCoordinatesToWorldCoordinates(Vector3D<float> (0.0F,
                                                                                   MaxYValue,
                                                                                   0.0F));

      Resize(pixelsX, pixelsY);

      ScaleX = (float) xPixels / (2.F * MaxXValue);
      ScaleY = (float) yPixels / (2.F * MaxYValue);

      Vector3D<float> radVector = viewpoint-> RotateCameraCoordinatesToWorldCoordinates(Vector3D<
          float> (0.F, 0.F, -rad));

      mVtx = radVector * 0.5F - UnitVectorProjectionX - UnitVectorProjectionY;

      UnitVectorProjectionX.x *= (2.F / (float) xPixels);
      UnitVectorProjectionX.y *= (2.F / (float) xPixels);
      UnitVectorProjectionX.z *= (2.F / (float) xPixels);

      UnitVectorProjectionY.x *= (2.F / (float) yPixels);
      UnitVectorProjectionY.y *= (2.F / (float) yPixels);
      UnitVectorProjectionY.z *= (2.F / (float) yPixels);
    }

    void Screen::Resize(unsigned int newPixelsX, unsigned int newPixelsY)
    {
      if (newPixelsX * newPixelsY <= COLOURED_PIXELS_MAX)
      {
        xPixels = newPixelsX;
        yPixels = newPixelsY;
      }
    }

    const Vector3D<float>& Screen::GetVtx() const
    {
      return mVtx;
    }
    const Vector3D<float>& Screen::GetUnitVectorProjectionX() const
    {
      return UnitVectorProjectionX;
    }
    const Vector3D<float>& Screen::GetUnitVectorProjectionY() const
    {
      return UnitVectorProjectionY;
    }
    int Screen::GetPixelsX() const
    {
      return xPixels;
    }
    int Screen::GetPixelsY() const
    {
      return yPixels;
    }
  }
}
