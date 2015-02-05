// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_UNITTESTS_HELPERS_FOLDERTESTFIXTURE_H

#include <cppunit/TestAssert.h>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "resources/Resource.h"
#include "util/utilityFunctions.h"
#include "unittests/helpers/HasCommsTestFixture.h"
#include <tinyxml.h>

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

          //! \brief Modify XML resource
          //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
          //! these hoops to test it.
          //! \param[in] resource: Name of the resource.
          //!   Should exist in the temp dir prior to call.
          //! \param[in] elements: hierarchy of elements + attribute (last item)
          //!   Should not include "hemelbsettings"
          //! \param[in] value: Value to set the attribute to
          template<class T>
            void ModifyXMLInput(
                std::string const &resource , std::vector<std::string>&& elements, T const &_value)
            {
              std::string const filename = tempPath + "/" + resource;
              std::string const attribute = elements.back();
              elements.pop_back();
              TiXmlDocument document(filename.c_str());
              document.LoadFile();
              auto child = document.FirstChildElement("hemelbsettings");
              for(std::string const &name: elements)
              {
                child = child->FirstChildElement(name);
              }
              std::ostringstream attr_value;
              attr_value << _value;
              child->SetAttribute(attribute, attr_value.str().c_str());

              std::ofstream output(filename);
              output << document;
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
