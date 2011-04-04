#ifndef HEMELB_VIS_SCREEN_H
#define HEMELB_VIS_SCREEN_H

#include "vis/ColPixel.h"

namespace hemelb
{
  namespace vis
  {
    class Screen
    {
      public:
        void AddPixel(const ColPixel* newPixel, lb::StressTypes iStressType, int mode);
        void RenderLine(const float endPoint1[3],
                        const float endPoint2[3],
                        lb::StressTypes iStressType,
                        int mode);

        /**
         * Does a transform from input array into output array. This function
         * will still work if the two arrays point to the same location in
         * memory. It only operates on the first two elements of input.
         *
         * @param input
         * @param output
         */
        template<typename T>
        void Transform(float* input, T output[2]) const
        {
          output[0] = (T) (ScaleX * (input[0] + MaxXValue));
          output[1] = (T) (ScaleY * (input[1] + MaxYValue));
        }

        static const unsigned int COLOURED_PIXELS_MAX = 2048 * 2048;

        float vtx[3];

        // Projection of unit vectors along screen axes into normal space.
        float UnitVectorProjectionX[3];
        float UnitVectorProjectionY[3];

        float MaxXValue, MaxYValue;
        float ScaleX, ScaleY;

        int PixelsX, PixelsY;

        // number of ColPixels.
        unsigned int col_pixels;
        // Array of pixel ids.
        int *col_pixel_id;
        ColPixel localPixels[COLOURED_PIXELS_MAX];

      private:
        template<bool xLimited> void RenderLineHelper(int x,
                                                      int y,
                                                      int incE,
                                                      int incNE,
                                                      int limit,
                                                      lb::StressTypes stressType,
                                                      int mode);
    };
  }
}

#endif /* HEMELB_VIS_SCREEN_H */
