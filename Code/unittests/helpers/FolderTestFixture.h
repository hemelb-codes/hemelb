#ifndef HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#include <cppunit/TestFixture.h>
#include <ctime>
namespace hemelb{
  namespace unittests{
    class FolderTestFixture: public CppUnit::TestFixture {

      public:
      FolderTestFixture(){
        std::stringstream temp_path_stream;
        temp_path_stream<<"HemeLBTest"<<time(NULL) << std::flush;
        temp_path=temp_path_stream.str();
      }
      void setUp(){
        // store current location
        origin=util::GetCurrentDir();
        // move to a temporary folder
        util::ChangeDirectory(util::GetTemporaryDir());
        // create a folder to work in
        util::MakeDirAllRXW(temp_path);
        util::ChangeDirectory(temp_path);

      }
      void tearDown(){
        // return to origin
        // doesn't matter not to clean up in tempdir.
        util::ChangeDirectory(origin);

      }
      protected:
        void AssertPresent(const std::string &fname){
          CPPUNIT_ASSERT(util::DoesDirectoryExist(fname.c_str()));
        }
      private:
      std::string origin;
      std::string temp_path;

    };
  }
}
#endif // ONCE
