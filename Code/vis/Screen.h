#ifndef HEMELB_VIS_SCREEN_H
#define HEMELB_VIS_SCREEN_H

#include "io/Writer.h"
#include "util/Vector3D.h"
#include "vis/ColPixel.h"
#include "vis/rayTracer/RayDataEnhanced.h"
#include "vis/rayTracer/RayDataNormal.h"
#include "vis/rayTracer/ClusterNormal.h"
#include "vis/rayTracer/ClusterWithWallNormals.h"
#include "vis/ScreenPixels.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"
#include "vis/XYCoordinates.h"

namespace hemelb
{
  namespace vis
  {
    //The ray data type and cluster type used for the ray tracer
    //typedef raytracer::RayDataNormal RayDataType_t;
    typedef raytracer::RayDataEnhanced<raytracer::DepthCuing::MIST> RayDataType_t;
    typedef raytracer::ClusterWithWallNormals ClusterType_t;

    class Screen
    {
        friend class Control;

      public:
        Screen();
        ~Screen();

        void AddPixel(const ColPixel<RayDataType_t>& newPixel,
		      const VisSettings& iVisSettings);

	void AddRayData( const XYCoordinates<int>& iPixelCoordinates, 
			 const RayDataType_t& iRayData,
			 const VisSettings& iVisSettings);
 
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

        const util::Vector3D<float>& GetCameraToBottomLeftOfScreenVector() const;
        const util::Vector3D<float>& GetPixelUnitVectorProjectionX() const;
        const util::Vector3D<float>& GetPixelUnitVectorProjectionY() const;
        int GetPixelsX() const;
        int GetPixelsY() const;

        ScreenPixels<RayDataType_t>* SwapBuffers(ScreenPixels<RayDataType_t>*);
        const ScreenPixels<RayDataType_t>* GetPixels() const;

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
	
        util::Vector3D<float> mCameraToBottomLeftOfScreen;

        // Projection of unit vectors along screen axes into normal space.
	util::Vector3D<float>  mPixelUnitVectorProjectionX;
        util::Vector3D<float> mPixelUnitVectorProjectionY;

        ScreenPixels<RayDataType_t>* mPixels;
    };
  }
}

#endif /* HEMELB_VIS_SCREEN_H */
