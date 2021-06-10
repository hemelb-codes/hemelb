// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_FOLDERTESTFIXTURE_H
#define HEMELB_TESTS_HELPERS_FOLDERTESTFIXTURE_H

#include <string>
#include <tinyxml.h>

#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace tests
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
      void ModifyXMLInput(TiXmlDocument &document, std::vector<std::string> const& elements,
                          T const &_value)
      {
        std::string const& attribute = elements.back();
	// Point to the *actual last elem*
	auto end = --elements.end();
	auto child = document.FirstChildElement("hemelbsettings");
	for (auto iter = elements.begin(); iter != end; ++iter) {
	  auto& name = *iter;
          auto next_child = child->FirstChildElement(name);
          if(next_child  == nullptr) {
            next_child = new TiXmlElement(name);
            child->LinkEndChild(next_child);
          }
          child = next_child;
        }
        std::ostringstream attr_value;
        attr_value << _value;
        child->SetAttribute(attribute, attr_value.str().c_str());
      }

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
