// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_SCREEN_H
#define HEMELB_VIS_SCREEN_H

#include "io/writers/Writer.h"
#include "util/Vector3D.h"
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

        void Set(float maxX, float maxY, int pixelsX, int pixelsY, float rad,
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
        XYCoordinates<T> TransformScreenToPixelCoordinates(const XYCoordinates<float>& iXYIn) const
        {
          return XYCoordinates<T>(static_cast<T>(mPixelsPerUnitX * (iXYIn.x + MaxXValue)),
                                  static_cast<T>(mPixelsPerUnitY * (iXYIn.y + MaxYValue)));
        }

        const util::Vector3D<float>& GetCameraToBottomLeftOfScreenVector() const;
        const util::Vector3D<float>& GetPixelUnitVectorProjectionX() const;
        const util::Vector3D<float>& GetPixelUnitVectorProjectionY() const;
        int GetPixelsX() const;
        int GetPixelsY() const;

      private:
        int xPixels, yPixels;

        //The number of pixels per unit of X or Y in screen coordinates
        float mPixelsPerUnitX;
        float mPixelsPerUnitY;

        //The extent of the screen in screen coordintes
        //(from -MaxXValue to MaxXValue)
        float MaxXValue;
        float MaxYValue;

        util::Vector3D<float> mCameraToBottomLeftOfScreen;

        // Projection of unit vectors along screen axes into normal space.
        util::Vector3D<float> mPixelUnitVectorProjectionX;
        util::Vector3D<float> mPixelUnitVectorProjectionY;
    };
  }
}

#endif /* HEMELB_VIS_SCREEN_H */
