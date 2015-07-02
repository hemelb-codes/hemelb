//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_IOTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_IOTESTS_H

#include <cppunit/TestFixture.h>
#include "unittests/helpers/FolderTestFixture.h"
#include "redblood/RBCInserter.h"
#include "redblood/io.h"
#include "redblood/FaderCell.h"

namespace hemelb
{
  namespace redblood
  {
    // Currently, this is private... So declare here.
    Cell::Moduli readModuli(io::xml::Element const& node, util::UnitConverter const &converter);
  }
  namespace unittests
  {
    namespace redblood
    {
      class CellIOTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (CellIOTests);
          CPPUNIT_TEST (testReadCellWithDefaults);
          CPPUNIT_TEST (testReadCellModuli);
          CPPUNIT_TEST (testReadMeshTemplates);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            converter.reset(new util::UnitConverter(0.5, 0.6, PhysicalPosition(1, 2, 3)));
            // It seems TiXML might take care of deallocation
            auto const parent = new TiXmlElement("parent");
            doc.LinkEndChild(parent);
            auto const cell = new TiXmlElement("cell");
            auto const shape = new TiXmlElement("shape");
            shape->SetAttribute("mesh_path", resources::Resource("red_blood_cell.txt").Path());
            cell->LinkEndChild(shape);
            auto const scaleXML = new TiXmlElement("scale");
            scale = 1.5;
            scaleXML->SetDoubleAttribute("value", scale);
            scaleXML->SetAttribute("units", "m");
            cell->LinkEndChild(scaleXML);

            auto const moduli = new TiXmlElement("moduli");
            auto add_stuff = [moduli](std::string name, std::string units, Dimensionless value)
            {
              auto elem = new TiXmlElement(name);
              elem->SetDoubleAttribute("value", value);
              elem->SetAttribute("units", units);
              moduli->LinkEndChild(elem);
            };
            add_stuff("surface", "LB", 2e0);
            add_stuff("dilation", "N/m", 1e0);
            add_stuff("bending", "Nm", 2e-18);
            cell->LinkEndChild(moduli);

            parent->LinkEndChild(cell);
          }

          // Reads cell with minimum stuff
          void testReadCellWithDefaults()
          {
            // Remove moduli, so we get default behavior
            doc.FirstChildElement("parent")->FirstChildElement("cell")->RemoveChild(doc.FirstChildElement("parent")->FirstChildElement("cell")->FirstChildElement("moduli"));

            auto cellbase = readCell(doc.FirstChildElement("parent"), *converter);
            std::unique_ptr<Cell const> const cell(static_cast<Cell const*>(cellbase.release()));
            auto const data = readMesh(resources::Resource("red_blood_cell.txt").Path());
            CPPUNIT_ASSERT_EQUAL(static_cast<site_t>(data->vertices.size()),
                                 cell->GetNumberOfNodes());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("m", scale),
                                         cell->GetScale(),
                                         1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, cell->moduli.volume, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, cell->moduli.surface, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("Nm", 2e-19),
                                         cell->moduli.bending,
                                         1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("N/m", 5e-1),
                                         cell->moduli.dilation,
                                         1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("N/m", 5e-6),
                                         cell->moduli.strain,
                                         1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, cell->nodeWall.intensity, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, cell->nodeWall.cutoff, 1e-12);
            CPPUNIT_ASSERT_EQUAL(site_t(2), site_t(cell->nodeWall.exponent));
          }

          void testReadCellModuli()
          {
            auto const moduli =
                readModuli(doc.FirstChildElement("parent")->FirstChildElement("cell"), *converter);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, moduli.volume, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(2e0, moduli.surface, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("Nm", 2e-18),
                                         moduli.bending,
                                         1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("N/m", 1e0),
                                         moduli.dilation,
                                         1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("N/m", 5e-6),
                                         moduli.strain,
                                         1e-12);
          }

          void testReadMeshTemplates()
          {
            auto const path = resources::Resource("red_blood_cell.txt").Path();
            std::ostringstream sstr;
            sstr << "<parent>"
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
                "        <shape mesh_path=\"" << path << "\"/>"
                "        <scale units=\"m\" value=\"0.6\"/>"
                "      </cell>"
                "     <cell name=\"joe\">"
                "       <shape mesh_path=\"" << path << "\"/>"
                "       <scale units=\"m\" value=\"0.5\"/>"
                "     </cell>"
                "   </cells>"
                "  </redbloodcells>"
                "</parent>";
            TiXmlDocument document;
            document.Parse(sstr.str().c_str());
            auto const cells = readTemplateCells(document.FirstChildElement("parent"), *converter);
            CPPUNIT_ASSERT_EQUAL(size_t(2), cells->size());
            CPPUNIT_ASSERT_EQUAL(size_t(1), cells->count("default"));
            CPPUNIT_ASSERT_EQUAL(size_t(1), cells->count("joe"));
            auto const default_ = std::static_pointer_cast<FaderCell>( (*cells)["default"]);
            auto const joe = std::static_pointer_cast<FaderCell>( (*cells)["joe"]);
            CPPUNIT_ASSERT(default_->GetTemplateName() == "default");
            CPPUNIT_ASSERT(joe->GetTemplateName() == "joe");
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("m", 0.6),
                                         default_->GetScale(),
                                         1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(converter->ConvertToLatticeUnits("m", 0.5),
                                         joe->GetScale(),
                                         1e-8);
            CPPUNIT_ASSERT(static_cast<bool>(default_->GetIOlets()));
            CPPUNIT_ASSERT(static_cast<bool>(joe->GetIOlets()));
            CPPUNIT_ASSERT(default_->GetIOlets() == joe->GetIOlets());
            CPPUNIT_ASSERT_EQUAL(size_t(2), joe->GetIOlets()->size());
          }
        private:
          TiXmlDocument doc;
          std::unique_ptr<util::UnitConverter> converter;
          PhysicalDistance scale;
      };


      CPPUNIT_TEST_SUITE_REGISTRATION (CellIOTests);
    }
  }
}

#endif  // ONCE
