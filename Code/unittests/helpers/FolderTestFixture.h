// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
#include "net/MpiCommunicator.h"
#include <tinyxml.h>

namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      //! \brief Modify XML document
      //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
      //! these hoops to test it.
      //! \param[in] document: Document to modify
      //! \param[in] elements: hierarchy of elements + attribute (last item)
      //!   Should not include "hemelbsettings"
      //! \param[in] value: Value to set the attribute to
      template<class T>
      void ModifyXMLInput(TiXmlDocument &document, std::vector<std::string>&& elements,
                          T const &_value)
      {
        std::string const attribute = elements.back();
        elements.pop_back();
        auto child = document.FirstChildElement("hemelbsettings");
        for (std::string const &name : elements)
        {
          auto next_child = child->FirstChildElement(name);
          if(next_child  == nullptr)
          {
            next_child = new TiXmlElement(name);
            child->LinkEndChild(next_child);
          }
          child = next_child;
        }
        std::ostringstream attr_value;
        attr_value << _value;
        child->SetAttribute(attribute, attr_value.str().c_str());
      }

      //! \brief Modify XML document by deleting an element if it exists
      //! \details HemeLB parameters cannot be modified programmatically, so we have to jump
      //! these hoops to test it.
      //! \param[in] document: Document to modify
      //! \param[in] elements: hierarchy of elements, last item will be removed.
      //!   Should not include "hemelbsettings"
      //! \param[in] value: Value to set the attribute to
      inline void DeleteXMLInput(TiXmlDocument &document, std::vector<std::string>&& elements)
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
            SyncTempPath(Comms());

            // create a folder to work in
            if(Comms().Rank() == 0)
            {
              util::MakeDirAllRXW(tempPath);
            }
            HEMELB_MPI_CALL(MPI_Barrier, (Comms()));
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
          void CopyResourceToTempdir(const std::string & resource) const
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
          void ModifyXMLInput(std::string const &resource, std::vector<std::string>&& elements,
                              T const &_value) const
          {
            std::string const filename = tempPath + "/" + resource;
            TiXmlDocument document(filename.c_str());
            document.LoadFile();
            helpers::ModifyXMLInput(document, std::move(elements), _value);
            std::ofstream output(filename);
            output << document;
          }

          void DeleteXMLInput(
              std::string const &resource, std::vector<std::string>&& elements) const
          {
            std::string const filename = tempPath + "/" + resource;
            TiXmlDocument document(filename.c_str());
            document.LoadFile();
            helpers::DeleteXMLInput(document, std::move(elements));
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

          void SyncTempPath(net::MpiCommunicator const &comm)
          {
            int n = tempPath.size();
            comm.Broadcast(n, 0);
            // no broadcast for strings defined...
            std::vector<int> value(n);
            if(comm.Rank() == 0)
            {
              for(int i(0); i < n; ++i)
              {
                value[i] = static_cast<int>(tempPath[i]);
              }
            }
            tempPath.resize(n);
            comm.Broadcast(value, 0);
            for(int i(0); i < n; ++i)
            {
              tempPath[i] = static_cast<char>(value[i]);
            }
          }
        private:
          std::string origin;
          std::string tempPath;
      };
    }
  }
}
#endif // ONCE
