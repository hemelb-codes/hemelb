#ifndef HEMELB_UNITTESTS_LBTESTS_LBTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_LBTESTS_H

#include <cppunit/TestCaller.h>
#include <cppunit/TestFixture.h>

#include "unittests/lbtests/KernelTests.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {

      /**
       * Class containing tests for the functionality of the lattice-Boltzmann kernels.
       */
      class LbTestSuite : public CppUnit::TestSuite
      {
        public:
          LbTestSuite()
          {
            addTest(new CppUnit::TestCaller<KernelTests>("runTest", &KernelTests::runTest));
          }

        private:
      };

    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_LBTESTS_H */
