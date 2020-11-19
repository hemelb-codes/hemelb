// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_TESTS_HELPERS_FOLDERTESTFIXTURE_H

#include <string>

#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      class FolderTestFixture : public HasCommsTestFixture
      {
      private:
	std::string origin;
	std::string tempPath;

      public:
	FolderTestFixture();
	~FolderTestFixture();

      protected:
	void ReturnToOrigin();
	void CopyResourceToTempdir(const std::string & resource);
	void MoveToTempdir();
	void AssertPresent(const std::string &fname);
	const std::string & GetTempdir();
      };
    }
  }
}
#endif // ONCE
