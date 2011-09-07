#ifndef HEMELB_UNITTESTS_LBTESTS_KERNELTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_KERNELTESTS_H

#include <cppunit/TestFixture.h>

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {

      /**
       * Class containing tests for the functionality of the lattice-Boltzmann kernels.
       */
      class KernelTests : public CppUnit::TestFixture
      {
        public:
          void setUp()
          {

          }
          void tearDown()
          {

          }

          void runTest()
          {
            CPPUNIT_ASSERT(1 == 0);
          }
      };

    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_KERNELTESTS_H */
