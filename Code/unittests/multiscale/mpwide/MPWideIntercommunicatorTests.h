#ifndef HEMELB_UNITTESTS_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATORTESTS_H
#define HEMELB_UNITTESTS_MULTISCALE_MPWIDE_MPWIDEINTERCOMMUNICATORTESTS_H
#include "MPWide.h"
namespace hemelb
{
  namespace unittests
  {
    namespace multiscale
    {
      /***
       * Mock intercommunicating entity which looks a bit like a HemeLB conceptually
       * It has an input, and an output, and the flow rate depends on the difference in pressures.
       */
      namespace mpwide
    {
    

      class MPWideIntercommunicatorTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE(MPWideIntercommunicatorTests);
          CPPUNIT_TEST(TestMPWidePresent);
    CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {

          }
          void tearDown()
          {
          }
        private:
          void TestMPWidePresent()
          {
             char *hostname;
             hostname = MPW_DNSResolve("localhost");
          }
      };
      //class
      CPPUNIT_TEST_SUITE_REGISTRATION(MPWideIntercommunicatorTests);
        }
     } //multiscale
  } //unittests
} //hemelb
#endif  //HEMELB_UNITTEST_MULTISCALE_MOCKINTERCOMMUNICATORTESTS_H
