// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "configuration/SimConfig.h"
#include "redblood/MeshIO.h"
#include "redblood/RBCInserter.h"
#include "redblood/xmlIO.h"
#include "redblood/FaderCell.h"

#include "tests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace redblood
  {
    // Currently, this is private... So declare here.
    Cell::Moduli readModuli(io::xml::Element const& node, util::UnitConverter const &converter);
  }
  namespace tests
  {
    using namespace redblood;

    class UninitialisedSimConfig : public configuration::SimConfig {
    public:
      UninitialisedSimConfig(const std::string& path) : configuration::SimConfig(path) {
      }
    };

    TEST_CASE_METHOD(helpers::FolderTestFixture, "CellIOTests", "[redblood]") {
      TiXmlDocument doc;
      LatticeDistance scale;

      auto converter = std::make_unique<util::UnitConverter>(0.5, 0.6, LatticePosition(1, 2, 3), 1000.0, 0.0);

      CopyResourceToTempdir("red_blood_cell.txt");
      CopyResourceToTempdir("empty_for_relative_paths.xml");
      auto config = std::make_unique<UninitialisedSimConfig>("empty_for_relative_paths.xml");

      {
	// It seems TiXML might take care of deallocation
	auto const parent = new TiXmlElement("parent");
	doc.LinkEndChild(parent);
	auto const cell = new TiXmlElement("cell");
	auto const shape = new TiXmlElement("shape");
	shape->SetAttribute("mesh_path", "red_blood_cell.txt");
	shape->SetAttribute("mesh_format", "Krueger");
	cell->LinkEndChild(shape);
	auto const scaleXML = new TiXmlElement("scale");
	scale = 1.5;
	scaleXML->SetDoubleAttribute("value", scale);
	scaleXML->SetAttribute("units", "m");
	cell->LinkEndChild(scaleXML);

	auto const moduli = new TiXmlElement("moduli");
	auto add_stuff = [moduli](std::string name, std::string units, Dimensionless value) {
	  auto elem = new TiXmlElement(name);
	  elem->SetDoubleAttribute("value", value);
	  elem->SetAttribute("units", units);
	  moduli->LinkEndChild(elem);
	};
	add_stuff("surface", "lattice", 2e0);
	add_stuff("dilation", "lattice", 0.58);
	add_stuff("bending", "Nm", 2e-18);
	cell->LinkEndChild(moduli);

	parent->LinkEndChild(cell);
      }
    
      auto approx = Approx(0.0).margin(1e-12);

      // Reads cell with minimum stuff
      SECTION("testReadCellWithDefaults") {
	// Remove moduli, so we get default behavior
	doc.FirstChildElement("parent")->FirstChildElement("cell")->RemoveChild(doc.FirstChildElement("parent")->FirstChildElement("cell")->FirstChildElement("moduli"));

	auto cellbase = readCell(doc.FirstChildElement("parent"), *config, *converter);
	std::unique_ptr<Cell const> const cell(static_cast<Cell const*>(cellbase.release()));
	auto const kruegerIO = redblood::KruegerMeshIO{};
	auto const data = kruegerIO.readFile("red_blood_cell.txt", true);
	REQUIRE(static_cast<site_t>(data->vertices.size()) == cell->GetNumberOfNodes());
	REQUIRE(approx(converter->ConvertToLatticeUnits("m", scale)) == cell->GetScale());
	REQUIRE(approx(1e0) == cell->moduli.volume);
	REQUIRE(approx(1e0) == cell->moduli.surface);
	REQUIRE(approx(converter->ConvertToLatticeUnits("Nm", 2e-19)) == cell->moduli.bending);
	REQUIRE(approx(0.75) == cell->moduli.dilation);
	REQUIRE(approx(converter->ConvertToLatticeUnits("N/m", 5e-6)) == cell->moduli.strain);
      }

      SECTION("testReadCellModuli") {
	auto const moduli =
	  readModuli(doc.FirstChildElement("parent")->FirstChildElement("cell"), *converter);
	REQUIRE(approx(1e0) == moduli.volume);
	REQUIRE(approx(2e0) == moduli.surface);
	REQUIRE(approx(converter->ConvertToLatticeUnits("Nm", 2e-18)) == moduli.bending);
	REQUIRE(approx(0.58) == moduli.dilation);
	REQUIRE(approx(converter->ConvertToLatticeUnits("N/m", 5e-6)) == moduli.strain);
      }

      SECTION("testReadMeshTemplates") {
	const char* xml_text = "<parent>"
	  "  <inlets>"
	  "   <inlet>"
	  "     <normal units=\"dimensionless\" value=\"(0.0,0.0,1.0)\" />"
	  "     <position units=\"m\" value=\"(0.0,0.0,-0.024)\" />"
	  "     <flowextension>"
	  "       <length units=\"m\" value=\"0.1\" />"
	  "       <radius units=\"m\" value=\"0.01\" />"
	  "       <fadelength units=\"m\" value=\"0.05\" />"
	  "     </flowextension>"
	  "   </inlet>"
	  "  </inlets>"
	  "  <outlets>"
	  "    <outlet>"
	  "      <normal units=\"dimensionless\" value=\"(0.0,0.0,-1.0)\" />"
	  "      <position units=\"m\" value=\"(0.0,0.0,0.024)\" />"
	  "      <flowextension>"
	  "        <length units=\"m\" value=\"0.1\" />"
	  "        <radius units=\"m\" value=\"0.01\" />"
	  "        <fadelength units=\"m\" value=\"0.05\" />"
	  "      </flowextension>"
	  "    </outlet>"
	  "  </outlets>"
	  "  <redbloodcells>"
	  "    <cells>"
	  "      <cell>"
	  "        <shape mesh_path=\"red_blood_cell.txt\" mesh_format=\"Krueger\" />"
	  "        <scale units=\"m\" value=\"0.6\"/>"
	  "      </cell>"
	  "     <cell name=\"joe\">"
	  "       <shape mesh_path=\"red_blood_cell.txt\" mesh_format=\"Krueger\" />"
	  "       <scale units=\"m\" value=\"0.5\"/>"
	  "     </cell>"
	  "   </cells>"
	  "  </redbloodcells>"
	  "</parent>";
	TiXmlDocument document;
	document.Parse(xml_text);
	auto const cells = readTemplateCells(document.FirstChildElement("parent"), *config, *converter);
	REQUIRE(size_t(2) == cells->size());
	REQUIRE(size_t(1) == cells->count("default"));
	REQUIRE(size_t(1) == cells->count("joe"));
	auto const default_ = std::static_pointer_cast<FaderCell>( (*cells)["default"]);
	auto const joe = std::static_pointer_cast<FaderCell>( (*cells)["joe"]);
	REQUIRE(default_->GetTemplateName() == "default");
	REQUIRE(joe->GetTemplateName() == "joe");
	REQUIRE(Approx(converter->ConvertToLatticeUnits("m", 0.6)).margin(1e-8)
		== default_->GetScale());
	REQUIRE(Approx(converter->ConvertToLatticeUnits("m", 0.5)).margin(1e-8)
		== joe->GetScale());
	REQUIRE(static_cast<bool>(default_->GetIOlets()));
	REQUIRE(static_cast<bool>(joe->GetIOlets()));
	REQUIRE(default_->GetIOlets() == joe->GetIOlets());
	REQUIRE(size_t(2) == joe->GetIOlets()->size());
      }
    }
  }
}
