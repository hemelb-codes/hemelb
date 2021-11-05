// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include <tinyxml.h>

#include "io/xml/XmlAbstractionLayer.h"
#include "redblood/FlowExtension.h"
#include "redblood/xmlIO.h"
#include "redblood/RBCInserter.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/FolderTestFixture.h"
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    TEST_CASE_METHOD(FlowExtensionFixture, "CellInserterTests", "[redblood]") {
      auto converter = std::make_unique<util::UnitConverter>(0.5, 0.6, PhysicalPosition::Zero(), 1000.0, 0.0);
      LatticeTime every = 10, offset = 5;
      TemplateCellContainer cells;
      cells.emplace("joe", std::make_shared<Cell>(tetrahedron()));

      auto getDocument = [&](LatticeDistance radius = 1e0, int numInserters = 0) -> TiXmlDocument {
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
      };

      SECTION("testNoPeriodicInsertion") {
	auto doc = getDocument(1e0, 0);
	*cells["joe"] *= 0.1e0;
	auto const inserter = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
					       *converter,
					       cells);
	REQUIRE(not inserter);
      }

      SECTION("testCellOutsideFlowExtension") {
	auto doc = getDocument(1e0, 1);
	REQUIRE_THROWS_AS(readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
					      *converter,
					      cells),
			  hemelb::Exception);
      }

      SECTION("testPeriodicInsertion") {
	// Creates an inserter and checks it exists
	auto doc = getDocument(1, 2);
	*cells["joe"] *= 0.1e0;
	auto const inserter = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
					       *converter,
					       cells);
	REQUIRE(inserter);

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
	for (int i(0); i < N; ++i) {
	  inserter(addCell);
	  REQUIRE(0 == num_calls);
	}

	// Now first call should result in adding a cell
	inserter(addCell);
	REQUIRE(2 == num_calls);
	num_calls = 0;

	// Now should not output a cell until "every" time has gone by
	// Do it thrice for good measure
	for (int k = 0; k < 3; ++k) {
	  for (int i(0); i < int(std::floor(every / dt)) - 1; ++i) {
	    inserter(addCell);
	    REQUIRE(0 == num_calls);
	  }
	  // And should call it!
	  inserter(addCell);
	  REQUIRE(2 == num_calls);
	  num_calls = 0;
	}
      }

      SECTION("testTranslation") {
	auto const identity = rotationMatrix(LatticePosition(0, 0, 1),
					     LatticePosition(0, 0, 1));
	RBCInserterWithPerturbation inserter([]()
					     { return true;},
					     cells["joe"]->clone(), identity, 0e0, 0e0, LatticePosition(2,
													0,
													0),
					     LatticePosition(0, 4, 0));

	auto const barycenter = cells["joe"]->GetBarycenter();
	for (size_t i(0); i < 500; ++i) {
	  auto const cell = inserter.drop();
	  auto const n = cell->GetBarycenter();
	  REQUIRE(std::abs(n.x - barycenter.x) <= 2e0);
	  REQUIRE(std::abs(n.y - barycenter.y) <= 4e0);
	  REQUIRE(Approx(barycenter.z).margin(1e-8) == n.z);
	}
      }

      SECTION("testSpatialOffset") {
	auto doc = getDocument(1e0, 1);
	helpers::ModifyXMLInput(doc,
				{ "inlets", "inlet", "insertcell", "offset", "value" },
				0.0);
	helpers::ModifyXMLInput(doc,
				{ "inlets", "inlet", "insertcell", "every", "value" },
				0.4);
	*cells["joe"] *= 0.1e0;
	CellContainer::value_type current_cell;
	auto addCell = [&current_cell](CellContainer::value_type cell) {
	  current_cell = cell;
	};
	auto const inserter = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
					       *converter,
					       cells);
	REQUIRE(inserter);
	inserter(addCell);
	REQUIRE(current_cell);
	auto const cell0 = current_cell;

	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "x", "units" }, "m");
	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "x", "value" }, 0.1);
	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "y", "units" }, "m");
	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "y", "value" }, 0.1);
	auto const insertTranslated = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
						       *converter,
						       cells);
	REQUIRE(insertTranslated);
	insertTranslated(addCell);
	REQUIRE(current_cell);
	auto const trans = cell0->GetBarycenter() - current_cell->GetBarycenter();
	auto approx = Approx(0.0).margin(1e-8);
	REQUIRE(approx(0) == trans.Dot(LatticePosition(0, 1, 1)));
	REQUIRE(approx(0.1 / 0.6 * 0.1 / 0.6 * 2e0) == trans.GetMagnitudeSquared());

	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "z", "units" }, "m");
	helpers::ModifyXMLInput(doc, { "inlets", "inlet", "insertcell", "z", "value" }, 0.1);
	auto const insertWithZ = readRBCInserters(doc.FirstChildElement("hemelbsettings")->FirstChildElement("inlets"),
						  *converter,
						  cells);
	REQUIRE(insertWithZ);
	insertWithZ(addCell);
	REQUIRE(current_cell);
	auto const transZ = cell0->GetBarycenter() - current_cell->GetBarycenter();
	REQUIRE(approx(0.1 / 0.6 * std::sqrt(2)) == transZ.Dot(LatticePosition(0, 1, 1)));
	REQUIRE(approx(0.1 / 0.6 * 0.1 / 0.6 * 3e0) == transZ.GetMagnitudeSquared());
      }

    }

  }
}
