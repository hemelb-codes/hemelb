#ifndef HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#include <cppunit/TestFixture.h>
#include <cmath>
#include <iomanip>
#include "unittests/resources/Resource.h"
namespace hemelb{
  namespace unittests{
    namespace helpers{
      class FolderTestFixture: public CppUnit::TestFixture {

        public:
        void setUp(){
          std::stringstream temp_path_stream;
          // next line is a hack to get the build working again
          // TODO: find a portable uuid solution. BOOST?
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
        void CopyResourceToTempdir(const std::string & resource){
          bool ok = util::FileCopy(resources::Resource(resource).Path().c_str(),(temp_path+"/"+resource).c_str());
          CPPUNIT_ASSERT(ok);
        }
        void MoveToTempdir(){
          util::ChangeDirectory(GetTempdir());
        }
        void AssertPresent(const std::string &fname){
          // "does directory exist" actually works for files too.
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
}
#endif // ONCE
