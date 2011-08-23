#ifndef HEMELB_VIS_SCREEN_H
#define HEMELB_VIS_SCREEN_H

#include "io/Writer.h"
#include "vis/ColPixel.h"
#include "vis/Vector3D.h"
#include "vis/ScreenPixels.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"
#include "vis/XYCoordinates.h"

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
        void RenderLine(const XYCoordinates<float>& endPoint1,
                            const XYCoordinates<float>& endPoint2,
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
	  XYCoordinates<T> TransformScreenToPixelCoordinates
	  (const XYCoordinates<float>& iXYIn) const
        {
          return XYCoordinates<T>(
	    static_cast<T>(ScaleX * (iXYIn.x + MaxXValue)),
	    static_cast<T>(ScaleY * (iXYIn.y + MaxYValue))
	    );
        }

        const Vector3D<float>& GetVtx() const;
        const Vector3D<float>& GetUnitVectorProjectionX() const;
        const Vector3D<float>& GetUnitVectorProjectionY() const;
        int GetPixelsX() const;
        int GetPixelsY() const;

        ScreenPixels* SwapBuffers(ScreenPixels*);
        const ScreenPixels* GetPixels() const;

        bool MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress);

        unsigned int GetPixelCount() const;

      private:
        float ScaleX, ScaleY;

	// 
        float MaxXValue, MaxYValue;
        Vector3D<float> mVtx;

        // Projection of unit vectors along screen axes into normal space.
	Vector3D<float>  UnitVectorProjectionX;
        Vector3D<float> UnitVectorProjectionY;

        ScreenPixels* mPixels;
    };
  }
}

#endif /* HEMELB_VIS_SCREEN_H */
