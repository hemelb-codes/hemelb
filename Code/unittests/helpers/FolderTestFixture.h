#ifndef HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#include <cppunit/TestFixture.h>
#include <cmath>
#include <iomanip>
namespace hemelb{
  namespace unittests{
    class FolderTestFixture: public CppUnit::TestFixture {

      public:
      void setUp(){
        std::stringstream temp_path_stream;
        // next line is a hack to get the build working again
        // I will try to find a portable uuid solution
        temp_path_stream<< util::GetTemporaryDir()<<"/"<<"HemeLBTest"<< std::fixed << floor(util::myClock()*100000) << std::flush;
        temp_path=temp_path_stream.str();
        // store current location
        origin=util::GetCurrentDir();

        // create a folder to work in
        util::MakeDirAllRXW(temp_path);
        MoveToTempdir();
      }

      void tearDown(){
        ReturnToOrigin();
        // doesn't matter not to clean up in tempdir.
      }

      protected:
      void ReturnToOrigin(){
        // return to origin
        util::ChangeDirectory(origin);
      }

      void MoveToTempdir(){
        util::ChangeDirectory(GetTempdir());
      }
      void AssertPresent(const std::string &fname){
        CPPUNIT_ASSERT(util::DoesDirectoryExist(fname.c_str()));
      }
      const std::string & GetTempdir(){
        return temp_path;
      }
      private:
      std::string origin;
      std::string temp_path;
    };
  }
}
#endif // ONCE
