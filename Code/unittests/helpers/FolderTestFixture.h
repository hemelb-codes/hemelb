
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H

#include <cmath>
#include <iomanip>
#include "resources/Resource.h"
#include "util/utilityFunctions.h"
#include "unittests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      class FolderTestFixture : public HasCommsTestFixture
      {

        public:
          void setUp()
          {
            HasCommsTestFixture::setUp();

            std::stringstream tempPathStream;
            // next line is a hack to get the build working again
            // TODO: find a portable uuid solution. BOOST?
            tempPathStream << util::GetTemporaryDir() << "/" << "HemeLBTest" << std::fixed
                << floor(util::myClock() * 100000) << std::flush;
            tempPath = tempPathStream.str();
            // store current location
            origin = util::GetCurrentDir();

            // create a folder to work in
            util::MakeDirAllRXW(tempPath);
            MoveToTempdir();
          }

          void tearDown()
          {
            ReturnToOrigin();
            // doesn't matter not to clean up in tempdir.
            HasCommsTestFixture::tearDown();
          }

        protected:
          void ReturnToOrigin()
          {
            // return to origin
            util::ChangeDirectory(origin);
          }
          void CopyResourceToTempdir(const std::string & resource)
          {
            // TODO this should use a filesystem-independent path join (BOOST)
            bool ok = util::FileCopy(resources::Resource(resource).Path().c_str(),
                                     (tempPath + "/" + resource).c_str());
            CPPUNIT_ASSERT(ok);
          }
          void MoveToTempdir()
          {
            util::ChangeDirectory(GetTempdir());
          }
          void AssertPresent(const std::string &fname)
          {
            // "does directory exist" actually works for files too.
            CPPUNIT_ASSERT(util::DoesDirectoryExist(fname.c_str()));
          }
          const std::string & GetTempdir()
          {
            return tempPath;
          }
        private:
          std::string origin;
          std::string tempPath;
      };
    }
  }
}
#endif // ONCE
