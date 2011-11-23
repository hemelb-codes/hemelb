#ifndef HEMELB_UNITTESTS_REPORTING_FILEMANAGER_H
#define HEMELB_UNITTESTS_REPORTING_FILEMANAGER_H
#include "reporting/FileManager.h"
#include "unittests/helpers/FolderTestFixture.h"
namespace hemelb
{
  namespace unittests
  {
    namespace reporting
    {
      using namespace hemelb::reporting;

      class FileManagerTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE(FileManagerTests);
          CPPUNIT_TEST(TestCreate);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            FolderTestFixture::setUp();
            int argc=9;
            char* argv[9];
            int processorCount=5;
            argv[0]="hemelb";
            argv[1]="-in";
            argv[2]="config.xml";
            argv[3]="-i";
            argv[4]="1";
            argv[5]="-s";
            argv[6]="1";
            argv[7]="-ss";
            argv[8]="1111";
            configuration::CommandLine cl=configuration::CommandLine(argc,argv);
            fileManager=new FileManager(cl,true,processorCount);
          }

          void tearDown()
          {
            FolderTestFixture::tearDown();
            delete fileManager;
          }

          void TestCreate(){
            AssertPresent("results");
            AssertPresent("results/Images");
            AssertPresent("results/Snapshots");
          }

        private:

          FileManager *fileManager;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(FileManagerTests);
    }
  }
}
#endif // ONCE
