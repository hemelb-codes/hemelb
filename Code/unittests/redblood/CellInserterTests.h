// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLINSERTERTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLINSERTERTESTS_H

#include <cppunit/TestFixture.h>
#include "io/xml/XmlAbstractionLayer.h"
#include "redblood/FlowExtension.h"
#include "redblood/xmlIO.h"
#include "redblood/RBCInserter.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/FolderTestFixture.h"
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
          CPPUNIT_TEST (testSpatialOffset);CPPUNIT_TEST_SUITE_END();

        public:

          void setUp()
          {
            converter.reset(new util::UnitConverter(0.5, 0.6, PhysicalPosition::Zero(), 1000.0, 0.0));
            every = 10;
            offset = 5;
            cells.emplace("joe", std::make_shared<Cell>(tetrahedron()));
          }

          TiXmlDocument getDocument(LatticeDistance radius = 1e0, int numInserters = 0)
          {
            std::ostringstream sstr;
            sstr << "<hemelbsettings>"
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
                "</hemelbsettings>";
            TiXmlDocument doc;
            doc.Parse(sstr.str().c_str());
            return doc;
          }
          void testNoPeriodicInsertion()
          {
            auto doc = getDocument(1e0, 0);
            *cells["joe"] *= 0.1e0;
            auto const inserter = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
                                                   *converter,
                                                   cells);
            CPPUNIT_ASSERT(not inserter);
          }

          void testCellOutsideFlowExtension()
          {
            auto doc = getDocument(1e0, 1);
            CPPUNIT_ASSERT_THROW(readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
                                                  *converter,
                                                  cells),
                                 hemelb::Exception);
          }
          void testPeriodicInsertion()
          {
            // Creates an inserter and checks it exists
            auto doc = getDocument(1, 2);
            *cells["joe"] *= 0.1e0;
            auto const inserter = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
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
            auto const identity = rotationMatrix(LatticePosition(0, 0, 1),
                                                 LatticePosition(0, 0, 1));
            RBCInserterWithPerturbation inserter([]()
            { return true;},
                                                 cells["joe"]->clone(), identity, 0e0, 0e0, LatticePosition(2,
                                                                                                            0,
                                                                                                            0),
                                                 LatticePosition(0, 4, 0));

            auto const barycenter = cells["joe"]->GetBarycenter();
            for (size_t i(0); i < 500; ++i)
            {
              auto const cell = inserter.drop();
              auto const n = cell->GetBarycenter();
              CPPUNIT_ASSERT(std::abs(n.x - barycenter.x) <= 2e0);
              CPPUNIT_ASSERT(std::abs(n.y - barycenter.y) <= 4e0);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenter.z, n.z, 1e-8);
            }
          }

          void testSpatialOffset()
          {
            auto doc = getDocument(1e0, 1);
            helpers::ModifyXMLInput(doc,
                                    { "inlets", "inlet", "insertcell", "offset", "value" },
                                    0.0);
            helpers::ModifyXMLInput(doc,
                                    { "inlets", "inlet", "insertcell", "every", "value" },
                                    0.4);
            *cells["joe"] *= 0.1e0;
            CellContainer::value_type current_cell;
            auto addCell = [&current_cell](CellContainer::value_type cell)
            {
              current_cell = cell;
            };
            auto const inserter = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
                                                   *converter,
                                                   cells);
            CPPUNIT_ASSERT(inserter);
            inserter(addCell);
            CPPUNIT_ASSERT(current_cell);
            auto const cell0 = current_cell;

            helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "x", "units" }, "m");
            helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "x", "value" }, 0.1);
            helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "y", "units" }, "m");
            helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "y", "value" }, 0.1);
            auto const insertTranslated = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
                                                           *converter,
                                                           cells);
            CPPUNIT_ASSERT(insertTranslated);
            insertTranslated(addCell);
            CPPUNIT_ASSERT(current_cell);
            auto const trans = cell0->GetBarycenter() - current_cell->GetBarycenter();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, trans.Dot(LatticePosition(0, 1, 1)), 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1 / 0.6 * 0.1 / 0.6 * 2e0,
                                         trans.GetMagnitudeSquared(),
                                         1e-8);

            helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "z", "units" }, "m");
            helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "z", "value" }, 0.1);
            auto const insertWithZ = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
                                                      *converter,
                                                      cells);
            CPPUNIT_ASSERT(insertWithZ);
            insertWithZ(addCell);
            CPPUNIT_ASSERT(current_cell);
            auto const transZ = cell0->GetBarycenter() - current_cell->GetBarycenter();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1 / 0.6 * std::sqrt(2),
                                         transZ.Dot(LatticePosition(0, 1, 1)),
                                         1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1 / 0.6 * 0.1 / 0.6 * 3e0,
                                         transZ.GetMagnitudeSquared(),
                                         1e-8);
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

#endif // HEMELB_UNITTESTS_REDBLOOD_CELLINSERTERTESTS_H
