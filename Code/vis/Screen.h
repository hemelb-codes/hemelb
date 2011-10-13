#ifndef HEMELB_VIS_SCREEN_H
#define HEMELB_VIS_SCREEN_H

#include "io/Writer.h"
#include "vis/Vector3D.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    class Screen
    {
      public:
        static const unsigned int COLOURED_PIXELS_MAX = 2048 * 2048;

        Screen();
        ~Screen();

        void Set(float maxX,
                 float maxY,
                 int pixelsX,
                 int pixelsY,
                 float rad,
                 const Viewpoint* viewpoint);

        void Resize(unsigned int pixelsX, unsigned int pixelsY);

        /**
         * Does a transform from input array into output array. This function
         * will still work if the two arrays point to the same location in
         * memory. It only operates on the first two elements of input.
         *
         * @param input
         * @param output
         */
        template<typename T>
        void Transform(const float& inputX, const float& inputY, T& outputX, T& outputY) const
        {
          outputX = (T) (ScaleX * (inputX + MaxXValue));
          outputY = (T) (ScaleY * (inputY + MaxYValue));
        }

        const Vector3D<float>& GetVtx() const;
        const Vector3D<float>& GetUnitVectorProjectionX() const;
        const Vector3D<float>& GetUnitVectorProjectionY() const;
        int GetPixelsX() const;
        int GetPixelsY() const;

        bool MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress);

      private:
        float ScaleX, ScaleY;
        int xPixels, yPixels;
        float MaxXValue, MaxYValue;
        Vector3D<float> mVtx;

        // Projection of unit vectors along screen axes into normal space.
        Vector3D<float> UnitVectorProjectionX;
        Vector3D<float> UnitVectorProjectionY;
    };
  }
}

#endif /* HEMELB_VIS_SCREEN_H */
