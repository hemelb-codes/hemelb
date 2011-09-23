#include <assert.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits>

//#include "vis/raytracer/HSLToRGBConverter.h"
#include "HSLToRGBConverter.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      //See http://en.wikipedia.org/wiki/HSL_and_HSV for the
      //formulae used
      //TODO: Possibly replace truncation with rounding?
      void HSLToRGBConverter::Convert(float iHue,
                                      float iSaturation,
                                      float iLightness,
                                      unsigned char oRGBColour[3])
      {
        assert(iHue >=0.0F && iHue < 360.0F);
        assert(iSaturation >= 0.0F && iSaturation <=1.0F);
        assert(iLightness >= 0.0F && iLightness <=1.0F);

        //All values are stored as 32-bit unsigned integers
        //stored within 16-bits, to allow multiplication between
        //two values 

        //All values are scaled up to be between 0 and the
        //maximum value of a 16-bit interger
        static const uint32_t DegreesScalar = std::numeric_limits<uint16_t>::max() / 360;

        static const uint32_t OtherScalar = std::numeric_limits<uint16_t>::max();

        uint32_t lHue = static_cast<uint32_t>(iHue * static_cast<float>(DegreesScalar));

uint32_t        lSaturation = static_cast<uint32_t>(iSaturation * static_cast<float>(OtherScalar));

uint32_t        lLightness = static_cast<uint32_t>(iLightness * static_cast<float>(OtherScalar));

int32_t        lTemp = 2 * static_cast<int32_t>(lLightness) - static_cast<int32_t>(OtherScalar);

        uint32_t lChroma = ( (OtherScalar - abs(lTemp)) * lSaturation) / OtherScalar;

        uint32_t lHuePrime = lHue / 60;

        lTemp = static_cast<int>(lHuePrime % (2 * DegreesScalar)) - static_cast<int>(DegreesScalar);

        uint32_t lIntermediate = (lChroma * (DegreesScalar - abs(lTemp))) / DegreesScalar;

        uint32_t lHueInt = lHuePrime / DegreesScalar;
        uint32_t lRed;
        uint32_t lGreen;
        uint32_t lBlue;

        switch (lHueInt)
        {
          case 0:
            lRed = lChroma;
            lGreen = lIntermediate;
            lBlue = 0.0F;
            break;

          case 1:
            lRed = lIntermediate;
            lGreen = lChroma;
            lBlue = 0.0F;
            break;

          case 2:
            lRed = 0.0F;
            lGreen = lChroma;
            lBlue = lIntermediate;
            break;

          case 3:
            lRed = 0.0F;
            lGreen = lIntermediate;
            lBlue = lChroma;
            break;

          case 4:
            lRed = lIntermediate;
            lGreen = 0.0F;
            lBlue = lChroma;
            break;

          case 5:
            lRed = lChroma;
            lGreen = 0.0F;
            lBlue = lIntermediate;
            break;

          default:
            assert(false);
        }

        static const uint32_t CharScalar = OtherScalar / 255;

        int32_t lMatcher = static_cast<int>(lLightness) - static_cast<int>(lChroma) / 2;

        oRGBColour[0] = static_cast<unsigned char>( (lRed + lMatcher) / CharScalar);oRGBColour
        [1] = static_cast<unsigned char>( (lGreen + lMatcher) / CharScalar);oRGBColour
        [2] = static_cast<unsigned char>( (lBlue + lMatcher) / CharScalar);}}
    }
  }
