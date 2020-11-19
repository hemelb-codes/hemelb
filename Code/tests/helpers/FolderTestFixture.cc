// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/FolderTestFixture.h"

#include <sstream>
#include <cmath>
#include <iomanip>

#include <catch2/catch.hpp>

#include "resources/Resource.h"
#include "util/utilityFunctions.h"


namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {

      FolderTestFixture::FolderTestFixture()
      {
	{
	  std::stringstream tempPathStream;
	  // next line is a hack to get the build working again
	  // TODO: find a portable uuid solution. BOOST?
	  tempPathStream << util::GetTemporaryDir() << "/" << "HemeLBTest" << std::fixed
			 << std::floor(util::myClock() * 100000) << std::flush;
	  tempPath = tempPathStream.str();
	}
	// store current location
	origin = util::GetCurrentDir();

	// create a folder to work in
	util::MakeDirAllRXW(tempPath);
	MoveToTempdir();
      }

      FolderTestFixture::~FolderTestFixture() {
	ReturnToOrigin();
	util::DeleteDirTree(tempPath);
      }

      void FolderTestFixture::ReturnToOrigin()
      {
	// return to origin
	util::ChangeDirectory(origin);
      }

      void FolderTestFixture::CopyResourceToTempdir(const std::string & resource)
      {
	// TODO this should use a filesystem-independent path join (BOOST)
	bool ok = util::FileCopy(resources::Resource(resource).Path().c_str(),
				 (tempPath + "/" + resource).c_str());
	
	REQUIRE(ok);
      }

      void FolderTestFixture::MoveToTempdir()
      {
	util::ChangeDirectory(GetTempdir());
      }

      void FolderTestFixture::AssertPresent(const std::string &fname)
      {
	// "does directory exist" actually works for files too.
	REQUIRE(util::DoesDirectoryExist(fname.c_str()));
      }
      const std::string & FolderTestFixture::GetTempdir()
      {
	return tempPath;
      }
    }
  }
}

