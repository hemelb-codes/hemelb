#ifndef HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#include <cppunit/TestFixture.h>
#include <ctime>
namespace hemelb{
  namespace unittests{
    class FolderTestFixture: public CppUnit::TestFixture {

      public:
      void setUp(){
        std::stringstream temp_path_stream;
        temp_path_stream<< util::GetTemporaryDir()<<"/"<<"HemeLBTest"<< util::GetUUID() << std::flush;
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
