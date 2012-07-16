// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_RESULTPIXEL_H
#define HEMELB_VIS_RESULTPIXEL_H

#include <cmath>

#include "util/utilityFunctions.h"
#include "vis/BasicPixel.h"
#include "vis/rayTracer/RayDataNormal.h"
#include "vis/streaklineDrawer/StreakPixel.h"
#include "vis/VisSettings.h"
#include "vis/DomainStats.h"

namespace hemelb
{
  namespace vis
  {
    class ResultPixel : public BasicPixel
    {
      public:
        ResultPixel(const BasicPixel* glyph);

        ResultPixel(const raytracer::RayDataNormal* ray);

        ResultPixel(const streaklinedrawer::StreakPixel* streak);

        const raytracer::RayDataNormal* GetRayPixel() const;

        void Combine(const ResultPixel& other);

        void WritePixel(unsigned* pixel_index,
                        unsigned char rgb_data[12],
                        const DomainStats& iDomainStats,
                        const VisSettings& visSettings) const;

        /*
         * Debugging function to output details about the pixel to the logger.
         */
        void LogDebuggingInformation() const;

      private:

        static void PickColour(float value, float colour[3]);

        static void MakePixelColour(int rawRed, int rawGreen, int rawBlue, unsigned char* dest);

        bool hasGlyph;
        const raytracer::RayDataNormal* normalRayPixel;
        const streaklinedrawer::StreakPixel* streakPixel;
    };
  }
}

#endif /* HEMELB_VIS_RESULTPIXEL_H */
