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
	    static_cast<T>(mPixelsPerUnitX * (iXYIn.x + mMaxXValue)),
	    static_cast<T>(mPixelsPerUnitY * (iXYIn.y + mMaxYValue))
	    );
        }

        const Vector3D<float>& GetCameraToBottomLeftOfScreenVector() const;
        const Vector3D<float>& GetPixelUnitVectorProjectionX() const;
        const Vector3D<float>& GetPixelUnitVectorProjectionY() const;
        int GetPixelsX() const;
        int GetPixelsY() const;

        ScreenPixels* SwapBuffers(ScreenPixels*);
        const ScreenPixels* GetPixels() const;

        bool MouseIsOverPixel(int mouseX, int mouseY, float* density, float* stress);

        unsigned int GetPixelCount() const;

      private:
        //The number of pixels per unit of X or Y in screen coordinates
	float mPixelsPerUnitX;
	float mPixelsPerUnitY;

	//The extent of the screen in screen coordintes 
	//(from -mmMaxXValue to mmMaxXValue)
        float mMaxXValue;
	float mMaxYValue;
	
        Vector3D<float> mCameraToBottomLeftOfScreen;

        // Projection of unit vectors along screen axes into normal space.
	Vector3D<float>  mPixelUnitVectorProjectionX;
        Vector3D<float> mPixelUnitVectorProjectionY;

        ScreenPixels* mPixels;
    };
  }
}

#endif /* HEMELB_VIS_SCREEN_H */
