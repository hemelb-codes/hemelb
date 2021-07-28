// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/FolderTestFixture.h"

#include <sstream>
#include <fstream>
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
      //! \brief Modify XML document by deleting an element if it exists
      //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
      //! these hoops to test it.
      //! \param[in] document: Document to modify
      //! \param[in] elements: hierarchy of elements, last item will be removed.
      //!   Should not include "hemelbsettings"
      //! \param[in] value: Value to set the attribute to
      void DeleteXMLInput(TiXmlDocument &document, std::vector<std::string> const& elements)
      {
        std::string const attribute = elements.back();
        auto element = document.FirstChildElement("hemelbsettings");
        for (std::string const &name : elements)
        {
          auto next_child = element->FirstChildElement(name);
          if(next_child  == nullptr)
          {
            return;
          }
          element = next_child;
        }
        element->Parent()->RemoveChild(element);
      }

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

      void FolderTestFixture::CopyResourceToTempdir(const std::string & resource) const
      {
	// TODO this should use a filesystem-independent path join (BOOST)
	bool ok = util::FileCopy(resources::Resource(resource).Path().c_str(),
				 (tempPath + "/" + resource).c_str());
	
	REQUIRE(ok);
      }

      void FolderTestFixture::DeleteXMLInput(std::string const &resource, std::vector<std::string> const& elements) const
      {
	std::string const filename = tempPath + "/" + resource;
	TiXmlDocument document(filename.c_str());
	document.LoadFile();
	helpers::DeleteXMLInput(document, elements);
	std::ofstream output(filename);
	output << document;
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

