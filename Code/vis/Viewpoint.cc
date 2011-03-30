#include "vis/Viewpoint.h"

namespace hemelb
{
  namespace vis
  {
    void Viewpoint::RotateToViewpoint(float iXIn,
                                      float iYIn,
                                      float iZIn,
                                      float* oXOut,
                                      float* oYOut,
                                      float* oZOut)
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

      Rotate(SinXRotation, CosXRotation, SinYRotation, CosYRotation, iXIn, iYIn, iZIn, oXOut,
             oYOut, oZOut);
    }

    void Viewpoint::Rotate(float sinX,
                           float cosX,
                           float sinY,
                           float cosY,
                           float xIn,
                           float yIn,
                           float zIn,
                           float* xOut,
                           float* yOut,
                           float* zOut)
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

      const float lTemp = zIn * cosX - yIn * sinX;

      *xOut = lTemp * sinY + xIn * cosY;
      *yOut = zIn * sinX + yIn * cosX;
      *zOut = lTemp * cosY - xIn * sinY;
    }

    void Viewpoint::Project(const float p1[], float p2[])
    {
      float x1[3], x2[3];

      for (int l = 0; l < 3; l++)
      {
        x1[l] = p1[l] - x[l];
      }

      Rotate(-SinYRotation, CosYRotation, -SinXRotation, CosXRotation, x1[1], x1[0], x1[2], &x2[1],
             &x2[0], &x2[2]);

      float temp = dist / (p2[2] = -x2[2]);

      p2[0] = temp * x2[0];
      p2[1] = temp * x2[1];
    }

  }
}
