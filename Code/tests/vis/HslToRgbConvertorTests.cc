// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "vis/rayTracer/HSLToRGBConverter.h"

namespace hemelb
{
  namespace tests
  {
    namespace {
      void DoColourTest(float iHue,
			float iSaturation,
			float iLightness,
			unsigned char iExpectedRed,
			unsigned char iExpectedGreen,
			unsigned char iExpectedBlue)
      {
	unsigned char lRGB[3];

	vis::raytracer::HSLToRGBConverter::Convert(iHue, iSaturation, iLightness, lRGB);

	UNSCOPED_INFO("HslRgbConvertor test differed on red channel");
	REQUIRE(iExpectedRed == lRGB[0]);
	UNSCOPED_INFO("HslRgbConvertor test differed on green channel");
	REQUIRE(iExpectedGreen == lRGB[1]);
	UNSCOPED_INFO("HslRgbConvertor test differed on blue channel");
	REQUIRE(iExpectedBlue == lRGB[2]);
      }
    }

    TEST_CASE("HslToRgbConvertorTests") {
      DoColourTest(0, 0, 0, 0, 0, 0);
      DoColourTest(1, 1, 1, 255, 255, 255);
      DoColourTest(250, 0.5, 0.5, 84, 63, 191);
    }

  }
}

