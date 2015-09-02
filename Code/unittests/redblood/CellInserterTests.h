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
      class CellInserterTests : public FlowExtensionFixture
      {
          CPPUNIT_TEST_SUITE (CellInserterTests);
          CPPUNIT_TEST (testNoPeriodicInsertion);
          CPPUNIT_TEST (testCellOutsideFlowExtension);
          CPPUNIT_TEST (testPeriodicInsertion);
          CPPUNIT_TEST (testTranslation);
          CPPUNIT_TEST_SUITE_END();

        public:

          void setUp()
          {
            converter.reset(new util::UnitConverter(0.5, 0.6, 0e0));
            every = 10;
            offset = 5;
            cells.emplace("joe", std::make_shared<Cell>(tetrahedron()));
          }

          TiXmlDocument getDocument(LatticeDistance radius = 1e0, int numInserters = 0)
          {
            std::ostringstream sstr;
            sstr << "<parent>"
                "<inlets><inlet>"
                "  <normal units=\"dimensionless\" value=\"(0.0,1.0,1.0)\" />"
                "  <position units=\"m\" value=\"(0.1,0.2,0.3)\" />"
                "  <flowextension>"
                "    <length units=\"m\" value=\"0.5\" />"
                "    <radius units=\"m\" value=\"" << radius << "\" />"
                "    <fadelength units=\"m\" value=\"0.4\" />"
                "  </flowextension>";
            for (int i = 0; i < numInserters; ++i)
            {
              sstr << "  <insertcell template=\"joe\">"
                  "    <every units=\"s\" value=\"" << every << "\"/>"
                  "    <offset units=\"s\" value=\"" << offset << "\"/>"
                  "  </insertcell>";
            }
            sstr << "</inlet></inlets>"
                "</parent>";
            TiXmlDocument doc;
            doc.Parse(sstr.str().c_str());
            return doc;
          }
          void testNoPeriodicInsertion()
          {
            auto doc = getDocument(1e0, 0);
            *cells["joe"] *= 0.1e0;
            auto const inserter = readRBCInserters(doc.FirstChildElement("parent"),
                                                   *converter,
                                                   cells);
            CPPUNIT_ASSERT(not inserter);
          }

          void testCellOutsideFlowExtension()
          {
            auto doc = getDocument(1e0, 1);
            CPPUNIT_ASSERT_THROW(readRBCInserters(doc.FirstChildElement("parent"),
                                                  *converter,
                                                  cells),
                                 hemelb::Exception);
          }
          void testPeriodicInsertion()
          {
            // Creates an inserter and checks it exists
            auto doc = getDocument(1, 2);
            *cells["joe"] *= 0.1e0;
            auto const inserter = readRBCInserters(doc.FirstChildElement("parent"),
                                                   *converter,
                                                   cells);
            CPPUNIT_ASSERT(inserter);

            // all calls up to offset result in node added cell
            int num_calls = 0;
            auto addCell = [&num_calls](CellContainer::value_type)
            {
              // This function is called by inserter with a new cell
              // It is meant to actually do the job of adding cells to the simulation
              ++num_calls;
            };
            auto const dt = converter->ConvertTimeToPhysicalUnits(1e0);
            auto const N = static_cast<int>(std::floor(offset / dt));
            for (int i(0); i < N; ++i)
            {
              inserter(addCell);
              CPPUNIT_ASSERT_EQUAL(0, num_calls);
            }

            // Now first call should result in adding a cell
            inserter(addCell);
            CPPUNIT_ASSERT_EQUAL(2, num_calls);
            num_calls = 0;

            // Now should not output a cell until "every" time has gone by
            // Do it thrice for good measure
            for (int k = 0; k < 3; ++k)
            {
              for (int i(0); i < int(std::floor(every / dt)) - 1; ++i)
              {
                inserter(addCell);
                CPPUNIT_ASSERT_EQUAL(0, num_calls);
              }
              // And should call it!
              inserter(addCell);
              CPPUNIT_ASSERT_EQUAL(2, num_calls);
              num_calls = 0;
            }
          }

          void testTranslation()
          {
            auto const identity = rotationMatrix(LatticePosition(0, 0, 1), LatticePosition(0, 0, 1));
            RBCInserterWithPerturbation inserter
            (
              []() { return true; }, cells["joe"]->clone(),
              identity, 0e0, 0e0,
              LatticePosition(2, 0, 0), LatticePosition(0, 4, 0)
            );

            auto const barycenter = cells["joe"]->GetBarycenter();
            for(size_t i(0); i < 500; ++i)
            {
              auto const cell = inserter.drop();
              auto const n = cell->GetBarycenter();
              CPPUNIT_ASSERT(std::abs(n.x - barycenter.x) <= 2e0);
              CPPUNIT_ASSERT(std::abs(n.y - barycenter.y) <= 4e0);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenter.z, n.z, 1e-8);
            }
          }

        private:
          std::unique_ptr<util::UnitConverter> converter;
          LatticeTime every, offset;
          TemplateCellContainer cells;
      };


      CPPUNIT_TEST_SUITE_REGISTRATION (CellInserterTests);
    } // namespace: redblood
  } // namespace: unittests
} // namespace: hemelb

#endif // HEMELB_UNITTESTS_REDBLOOD_FLOWEXTENSION_TESTS_H