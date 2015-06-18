//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLINSERTER_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLINSERTER_TESTS_H

#include <cppunit/TestFixture.h>
#include "io/xml/XmlAbstractionLayer.h"
#include "redblood/FlowExtension.h"
#include "redblood/io.h"
#include "redblood/RBCInserter.h"
#include "unittests/redblood/Fixtures.h"
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellInserterTests : public FlowExtensionFixture {
        CPPUNIT_TEST_SUITE (CellInserterTests);
        CPPUNIT_TEST (testNoPeriodicInsertion);
        CPPUNIT_TEST (testCellOutsideFlowExtension);
        CPPUNIT_TEST (testPeriodicInsertion);
        CPPUNIT_TEST_SUITE_END();

        public:

          void setUp()
          {
            converter.reset(new util::UnitConverter(0.5, 0.6, 0e0));
            every = 10;
            offset = 5;
          }


          TiXmlDocument getDocument(
              PhysicalDistance radius = 1e0, PhysicalDistance scale = 1e0, bool noInserter=false)
          {
            std::ostringstream sstr;
            sstr <<
               "<parent>"
               "<inlets><inlet>"
               "  <normal units=\"dimensionless\" value=\"(0.0,1.0,1.0)\" />"
               "  <position units=\"m\" value=\"(0.1,0.2,0.3)\" />"
               "  <flowextension>"
               "    <length units=\"m\" value=\"0.5\" />"
               "    <radius units=\"m\" value=\"" << radius << "\" />"
               "    <fadelength units=\"m\" value=\"0.4\" />"
               "  </flowextension>";
            if(not noInserter)
            {
              sstr << "  <insertcell>"
                "    <every units=\"s\" value=\"" << every << "\"/>"
                "    <offset units=\"s\" value=\"" << offset << "\"/>"
                "  </insertcell>";
            }
            sstr << "</inlet></inlets>"
               "<redbloodcells>"
               "  <cell>"
               "    <shape mesh_path=\"" << resources::Resource("red_blood_cell.txt").Path()
                  << "\" />"
               "    <scale units=\"m\" value=\"" << scale << "\" />"
               "  </cell>"
               "</redbloodcells>"
               "</parent>";
            TiXmlDocument doc;
            doc.Parse(sstr.str().c_str());
            return doc;
          }
          void testNoPeriodicInsertion()
          {
            auto doc = getDocument(1e0, 0.1e0, true);
            auto const inserter = read_rbcinserter(doc.FirstChildElement("parent"), *converter);
            CPPUNIT_ASSERT(not inserter);
          }

          void testCellOutsideFlowExtension()
          {
            auto doc = getDocument(1e0, 1e0, false);
            CPPUNIT_ASSERT_THROW(
                read_rbcinserter(doc.FirstChildElement("parent"), *converter), hemelb::Exception);
          }
          void testPeriodicInsertion()
          {
            // Creates an inserter and checks it exists
            auto doc = getDocument(1, 0.1, false);
            auto const inserter = read_rbcinserter(doc.FirstChildElement("parent"), *converter);
            CPPUNIT_ASSERT(inserter);

            // all calls up to offset result in node added cell
            bool was_called = false;
            auto addCell = [&was_called](CellContainer::value_type)
            {
              // This function is called by inserter with a new cell
              // It is meant to actually do the job of adding cells to the simulation
              was_called = true;
            };
            auto const dt = converter->ConvertTimeToPhysicalUnits(1e0);
            for(int i(0); i < int(std::floor(offset / dt)); ++i)
            {
              (*inserter)(addCell);
              CPPUNIT_ASSERT(not was_called);
            }

            // Now first call should result in adding a cell
            (*inserter)(addCell);
            CPPUNIT_ASSERT(was_called);
            was_called = false;

            // Now should not output a cell until "every" time has gone by
            // Do it thrice for good measure
            for(int k = 0; k < 3; ++k)
            {
              for(int i(0); i < int(std::floor(every / dt)) - 1; ++i)
              {
                (*inserter)(addCell);
                CPPUNIT_ASSERT(not was_called);
              }
              // And should call it!
              (*inserter)(addCell);
              CPPUNIT_ASSERT(was_called);
              was_called = false;
            }
          }
        private:
          std::unique_ptr<util::UnitConverter> converter;
          PhysicalTime every, offset;
      };


      CPPUNIT_TEST_SUITE_REGISTRATION (CellInserterTests);
    } // namespace: redblood
  } // namespace: unittests
} // namespace: hemelb

#endif // HEMELB_UNITTESTS_REDBLOOD_FLOWEXTENSION_TESTS_H
