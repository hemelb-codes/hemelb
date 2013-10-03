//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

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
          CPPUNIT_TEST( TestAttributeConversion);
          CPPUNIT_TEST( TestAttributeConversionFails);
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

#define IF_TYPE_EQ_THEN_CHECK(DTYPE) \
    if (type == #DTYPE) \
    { \
      DTYPE value; \
      datum.GetAttributeOrThrow("value", value); \
    }
          void TestAttributeConversion()
          {
            xml::Element shouldWork = xmlDoc->GetRoot().GetChildOrThrow("conversiontests").GetChildOrThrow("shouldwork");

            for (xml::Element datum = shouldWork.GetChildOrThrow("datum");
                datum != xml::Element::Missing();
                datum = datum.NextSiblingOrNull("datum"))
            {
              std::string type = datum.GetAttributeOrThrow("type");
              IF_TYPE_EQ_THEN_CHECK(int);
              IF_TYPE_EQ_THEN_CHECK(double);
              IF_TYPE_EQ_THEN_CHECK(hemelb::util::Vector3D<double>);
              //IF_TYPE_EQ_THEN_CHECK(hemelb::util::Vector3D<int>);
              IF_TYPE_EQ_THEN_CHECK(unsigned);
            }
          }

          void TestAttributeConversionFails()
          {
            xml::Element shouldFail = xmlDoc->GetRoot().GetChildOrThrow("conversiontests").GetChildOrThrow("shouldfail");
            for (xml::Element datum = shouldFail.GetChildOrThrow("datum");
                datum != xml::Element::Missing();
                datum = datum.NextSiblingOrNull("datum"))
            {
              std::string type = datum.GetAttributeOrThrow("type");
              bool didThrow = false;
              try
              {
                IF_TYPE_EQ_THEN_CHECK(int);
                IF_TYPE_EQ_THEN_CHECK(double);
                IF_TYPE_EQ_THEN_CHECK(hemelb::util::Vector3D<double>);
                IF_TYPE_EQ_THEN_CHECK(hemelb::util::Vector3D<int>);
                IF_TYPE_EQ_THEN_CHECK(unsigned);
              } catch (std::exception& e) {
                std::cout << e.what() << std::endl;
                didThrow = true;
              }
              CPPUNIT_ASSERT(didThrow);
            }
          }
        private:
          xml::Document* xmlDoc;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(XmlTests);
    }
  }
}
#endif // HEMELB_UNITTESTS_IO_XML_H
