#ifndef HEMELB_VIS_SCREEN_H
#define HEMELB_VIS_SCREEN_H

#include "io/Writer.h"
#include "vis/ColPixel.h"
#include "vis/Location.h"
#include "vis/ScreenPixels.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    class Screen
    {
        friend class Control;

      public:
        Screen();
        ~Screen();

        void AddPixel(const ColPixel* newPixel, const VisSettings* visSettings);
        void RenderLine(const Location<float>& endPoint1,
			const Location<float>& endPoint2,
                        const VisSettings* visSettings);

        void Set(float maxX,
                 float maxY,
                 int pixelsX,
                 int pixelsY,
                 float rad,
                 const Viewpoint* viewpoint);

        void Resize(unsigned int pixelsX, unsigned int pixelsY);

        void Reset();

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

        const Location<float>& GetVtx() const;
        const Location<float>& GetUnitVectorProjectionX() const;
        const Location<float>& GetUnitVectorProjectionY() const;
        int GetPixelsX() const;
        int GetPixelsY() const;

        ScreenPixels* SwapBuffers(ScreenPixels*);
        const ScreenPixels* GetPixels() const;

        bool MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress);

        unsigned int GetPixelCount() const;

      private:
        float ScaleX, ScaleY;
        float MaxXValue, MaxYValue;
        Location<float> mVtx;

        // Projection of unit vectors along screen axes into normal space.
	Location<float>  UnitVectorProjectionX;
        Location<float> UnitVectorProjectionY;

        ScreenPixels* pixels;
    };
  }
}

#endif /* HEMELB_VIS_SCREEN_H */
