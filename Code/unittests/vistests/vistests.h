#ifndef HEMELB_UNITTESTS_VISTESTS_VISTESTS_H
#define HEMELB_UNITTESTS_VISTESTS_VISTESTS_H

#include <cppunit/TestCaller.h>
#include <cppunit/TestFixture.h>

#include "unittests/vistests/HslToRgbConvertorTests.h"

namespace hemelb
{
  namespace unittests
  {
    namespace vistests
    {

      /**
       * Class containing tests for the functionality of the vis system.
       */
      class VisTestSuite : public CppUnit::TestSuite
      {
        public:
          VisTestSuite()
          {
            addTest(new CppUnit::TestCaller<HslToRgbConvertorTests>("TestColours",
                                                                    &HslToRgbConvertorTests::TestColours));
          }
      };

    }
  }
}

#endif /* HEMELB_UNITTESTS_VISTESTS_VISTESTS_H */
