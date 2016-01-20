
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_IO_XML_H
#define HEMELB_UNITTESTS_IO_XML_H

#include "io/xml/XmlAbstractionLayer.h"
#include "unittests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace io
    {
      namespace xml = hemelb::io::xml;

      class XmlTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE( XmlTests);
          CPPUNIT_TEST( TestRead);
          CPPUNIT_TEST( TestSiblings);
          CPPUNIT_TEST( TestGetParent);
          CPPUNIT_TEST( TestGetChildNull);
          CPPUNIT_TEST_EXCEPTION( TestGetChildThrows, xml::ChildError);
          CPPUNIT_TEST( TestGetParentNull);
          CPPUNIT_TEST_EXCEPTION(TestGetParentThrows, xml::ParentError);
          CPPUNIT_TEST( TestAttributeConversion0);
          CPPUNIT_TEST( TestAttributeConversion1);
          CPPUNIT_TEST( TestAttributeConversion2);
          CPPUNIT_TEST( TestAttributeConversion3);
          CPPUNIT_TEST( TestAttributeConversion4);
          CPPUNIT_TEST( TestAttributeConversion5);
          CPPUNIT_TEST( TestAttributeConversion6);
          CPPUNIT_TEST( TestAttributeConversionFails0);
          CPPUNIT_TEST( TestAttributeConversionFails1);
          CPPUNIT_TEST( TestAttributeConversionFails2);
          CPPUNIT_TEST( TestAttributeConversionFails3);
          CPPUNIT_TEST( TestAttributeConversionFails4);
          CPPUNIT_TEST( TestAttributeConversionFails5);
          CPPUNIT_TEST( TestAttributeConversionFails6);
          CPPUNIT_TEST( TestAttributeConversionFails7);
          CPPUNIT_TEST( TestAttributeConversionFails8);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            helpers::FolderTestFixture::setUp();
            const std::string testFile("xmltest.xml");
            CopyResourceToTempdir(testFile);
            xmlDoc = new xml::Document(testFile);
          }

          void tearDown()
          {
            delete xmlDoc;
            helpers::FolderTestFixture::tearDown();
          }

          void TestRead()
          {
            xml::Element html = xmlDoc->GetRoot();

            // Check that the root elem is html
            CPPUNIT_ASSERT(html != xml::Element::Missing());
            CPPUNIT_ASSERT_EQUAL(std::string("html"), html.GetName());

            // Check a bunch of elements are present.
            {
              xml::Element head = html.GetChildOrThrow("head");
              {
                xml::Element title = head.GetChildOrThrow("title");
              }
              xml::Element body = html.GetChildOrThrow("body");
              {
                xml::Element div_banner = body.GetChildOrThrow("div");
                CPPUNIT_ASSERT_EQUAL(std::string("banner"), div_banner.GetAttributeOrThrow("id"));
                {
                  xml::Element div_header = div_banner.GetChildOrThrow("div");
                  CPPUNIT_ASSERT_EQUAL(std::string("header"), div_header.GetAttributeOrThrow("id"));

                }
              }
            }
          }
          void TestSiblings()
          {
            xml::Element banner = xmlDoc->GetRoot().GetChildOrThrow("body").GetChildOrThrow("div");
            int n = 0;
            std::string expectedIds[] = { "header", "metanav" };

            for (xml::Element div = banner.GetChildOrThrow("div"); div != xml::Element::Missing(); div
                = div.NextSiblingOrNull("div"))
            {
              CPPUNIT_ASSERT_EQUAL(expectedIds[n], div.GetAttributeOrThrow("id"));
              ++n;
            }
            CPPUNIT_ASSERT_EQUAL(2, n);
          }

          void TestGetChildNull()
          {
            xml::Element html = xmlDoc->GetRoot();
            CPPUNIT_ASSERT(html.GetChildOrNull("NoSuchElement") == xml::Element::Missing());
          }

          void TestGetChildThrows()
          {
            xml::Element html = xmlDoc->GetRoot();
            html.GetChildOrThrow("NoSuchElement");
          }

          void TestGetParent()
          {
            xml::Element html = xmlDoc->GetRoot();
            xml::Element body = html.GetChildOrThrow("body");

            xml::Element shouldBeHtml = body.GetParentOrThrow();
            CPPUNIT_ASSERT(html == shouldBeHtml);
          }

          void TestGetParentNull()
          {
            xml::Element html = xmlDoc->GetRoot();
            CPPUNIT_ASSERT(xml::Element::Missing() == html.GetParentOrNull());
          }

          void TestGetParentThrows()
          {
            xml::Element html = xmlDoc->GetRoot();
            html.GetParentOrThrow();
          }
	  
          // Save a load of typing with this macro
#define MAKE_TESTATTR(i, TYPE, EXPECTED)		     \
          void TestAttributeConversion##i () \
          { \
            xml::Element shouldWork = xmlDoc->GetRoot().GetChildOrThrow("conversiontests").GetChildOrThrow("shouldwork"); \
            xml::Element datum = GetDatum(shouldWork, i); \
            std::string type = datum.GetAttributeOrThrow("type"); \
            CPPUNIT_ASSERT_EQUAL(std::string(#TYPE), type); \
            TYPE value; \
            datum.GetAttributeOrThrow("value", value); \
            CPPUNIT_ASSERT_EQUAL(EXPECTED, value); \
          }
          // <datum type="int" value="120" />
          MAKE_TESTATTR(0, int, 120);
          // <datum type="int" value="-24324" />
          MAKE_TESTATTR(1, int, -24324);
          // <datum type="double" value="1.0" />
          MAKE_TESTATTR(2, double, 1.0);
          // <datum type="double" value="1.6e-3" />
          MAKE_TESTATTR(3, double, 1.6e-3);
          // <datum type="hemelb::util::Vector3D<double>" value="(-1.4,11e7,42)" />
          MAKE_TESTATTR(4, hemelb::util::Vector3D<double>, hemelb::util::Vector3D<double>(-1.4, 11e7, 42));
          // <datum type="hemelb::util::Vector3D<int>" value="(-1,11,42)" />
          MAKE_TESTATTR(5, hemelb::util::Vector3D<int>, hemelb::util::Vector3D<int>(-1, 11, 42));
          // <datum type="unsigned" value="42" />
          MAKE_TESTATTR(6, unsigned, 42U);
	  	  
#define MAKE_TESTATTRFAIL(i, TYPE)		  \
          void TestAttributeConversionFails##i () \
          {								\
	          xml::Element shouldFail = xmlDoc->GetRoot().GetChildOrThrow("conversiontests").GetChildOrThrow("shouldfail"); \
	          xml::Element datum = GetDatum(shouldFail, i);		\
	          std::string type = datum.GetAttributeOrThrow("type");	\
	          CPPUNIT_ASSERT_EQUAL(std::string(#TYPE), type);		\
	          TYPE value;							\
	          bool didThrow = false;					\
	          try								\
	          {								\
	        	  datum.GetAttributeOrThrow("value", value);		\
	          }								\
	          catch (hemelb::io::xml::ParseError& e)			\
	          {								\
	        	  std::cout << e.what() << std::endl;			\
	        	  didThrow = true;						\
	          }								\
	          CPPUNIT_ASSERT(didThrow);					\
          }

          MAKE_TESTATTRFAIL(0, double)
          MAKE_TESTATTRFAIL(1, double)
          MAKE_TESTATTRFAIL(2, int)
          MAKE_TESTATTRFAIL(3, hemelb::util::Vector3D<double>)
          MAKE_TESTATTRFAIL(4, hemelb::util::Vector3D<double>)
          MAKE_TESTATTRFAIL(5, hemelb::util::Vector3D<int>)
          MAKE_TESTATTRFAIL(6, hemelb::util::Vector3D<double>)
          MAKE_TESTATTRFAIL(7, hemelb::util::Vector3D<int>)
          MAKE_TESTATTRFAIL(8, unsigned);

        private:
          xml::Document* xmlDoc;
          xml::Element GetDatum(xml::Element& parent, unsigned iRequired)
          {
            unsigned i = 0;
            for (xml::Element datum = parent.GetChildOrThrow("datum");
                datum != xml::Element::Missing();
                datum = datum.NextSiblingOrNull("datum"))
              {
              if (i == iRequired)
                return datum;
              ++i;
              }
            throw Exception() << "Cannot find element 'datum' with required index = " << iRequired;
          }

      };

      CPPUNIT_TEST_SUITE_REGISTRATION(XmlTests);
    }
  }
}
#endif // HEMELB_UNITTESTS_IO_XML_H
