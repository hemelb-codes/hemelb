// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_RAYTRACER_HSLTORGBCONVERTER_H
#define HEMELB_VIS_RAYTRACER_HSLTORGBCONVERTER_H

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      class HSLToRGBConverter
      {
        public:
          /**
           * Converts a colour in HSL coordinates to
           * RGB coordinates (between 0 and 255)
           * iHue must be between 0.0F and 360.0F in degrees
           * and the other two between 0.0F and 1.0F
           */
          static void Convert(float iHue, float iSaturation, float iLightness,
                              unsigned char oRGBColour[3]);
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_HSLTORGBCONVERTER_H
