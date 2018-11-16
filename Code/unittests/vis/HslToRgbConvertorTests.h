
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_VIS_HSLTORGBCONVERTORTESTS_H
#define HEMELB_UNITTESTS_VIS_HSLTORGBCONVERTORTESTS_H

#include <cppunit/TestFixture.h>

namespace hemelb
{
  namespace unittests
  {
    namespace vis
    {
      class HslToRgbConvertorTests : public CppUnit::TestFixture
      {
        CPPUNIT_TEST_SUITE( HslToRgbConvertorTests );
        CPPUNIT_TEST( TestColours );
        CPPUNIT_TEST_SUITE_END();
        public:
          void TestColours()
          {
            DoColourTest(0, 0, 0, 0, 0, 0);
            DoColourTest(1, 1, 1, 255, 255, 255);
            DoColourTest(250, 0.5, 0.5, 84, 63, 191);
          }

        private:
          void DoColourTest(float iHue,
                            float iSaturation,
                            float iLightness,
                            unsigned char iExpectedRed,
                            unsigned char iExpectedGreen,
                            unsigned char iExpectedBlue)
          {
            unsigned char lRGB[3];

            hemelb::vis::raytracer::HSLToRGBConverter::Convert(iHue, iSaturation, iLightness, lRGB);

            CPPUNIT_ASSERT_EQUAL_MESSAGE("HslRgbConvertor test differed on red channel",
                                         iExpectedRed,
                                         lRGB[0]);
            CPPUNIT_ASSERT_EQUAL_MESSAGE("HslRgbConvertor test differed on green channel",
                                         iExpectedGreen,
                                         lRGB[1]);
            CPPUNIT_ASSERT_EQUAL_MESSAGE("HslRgbConvertor test differed on blue channel",
                                         iExpectedBlue,
                                         lRGB[2]);
          }

      };
      CPPUNIT_TEST_SUITE_REGISTRATION( HslToRgbConvertorTests );
    }
  }
}

#endif
