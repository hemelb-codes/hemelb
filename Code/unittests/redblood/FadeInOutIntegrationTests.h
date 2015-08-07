#ifndef HEMELB_UNITTESTS_REDBLOOD_FADEINOUTINTEGRATION_H
#define HEMELB_UNITTESTS_REDBLOOD_FADEINOUTINTEGRATION_H

#include <cppunit/extensions/HelperMacros.h>
#include <boost/uuid/uuid_io.hpp>
#include <memory>

#include "SimulationMaster.h"
#include "lb/BuildSystemInterface.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "unittests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class FadeInOutIntegrationTests : public hemelb::unittests::helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (FadeInOutIntegrationTests);
          CPPUNIT_TEST (testIntegration);CPPUNIT_TEST_SUITE_END()
          ;

          typedef Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type Traits;
          typedef hemelb::redblood::CellController<Traits::Kernel> CellControl;
          typedef SimulationMaster<Traits> MasterSim;

        public:
          void setUp()
          {
            FolderTestFixture::setUp();
            CopyResourceToTempdir("large_cylinder_rbc.xml");
            CopyResourceToTempdir("large_cylinder.gmy");
            CopyResourceToTempdir("red_blood_cell.txt");

            std::vector<std::string> intel;
            intel.push_back("simulation");
            intel.push_back("steps");
            intel.push_back("value");
            ModifyXMLInput("large_cylinder_rbc.xml", std::move(intel), 1000);

            argv[0] = "hemelb";
            argv[1] = "-in";
            argv[2] = "large_cylinder_rbc.xml";
            argv[3] = "-i";
            argv[4] = "1";
            argv[5] = "-ss";
            argv[6] = "1111";
            options = std::make_shared<hemelb::configuration::CommandLine>(argc, argv);

            master = std::make_shared<MasterSim>(*options, Comms());
          }

          void testIntegration()
          {
            auto const & converter = master->GetUnitConverter();
            auto const volumeFactor = std::pow(converter.ConvertToLatticeUnits("m", 1e0), -3)*1e12;
            auto checkVolume = [volumeFactor]( const hemelb::redblood::CellContainer & cells)
            {
              static PhysicalVolume expected = -1e0;
              if(cells.empty())
              {
                return;
              }
              CPPUNIT_ASSERT_EQUAL(int(1), int(cells.size()));
              auto cell = *cells.begin();
              auto const volume = cell->GetVolume() * volumeFactor;
              if(expected < 0e0)
              {
                expected = volume;
              }
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, volume, 1e-8);
            };
            auto checkPosition = []( const hemelb::redblood::CellContainer & cells)
            {
              if(cells.empty())
              {
                return;
              }
              CPPUNIT_ASSERT_EQUAL(int(1), int(cells.size()));
              auto cell = *cells.begin();
              static LatticePosition first, current, tenth;
              static int iter = 0;
              if(iter == 0)
              {
                first = cell->GetBarycenter();
                current = first;
                tenth = first;
              }
              auto const position = cell->GetBarycenter();
              CPPUNIT_ASSERT_DOUBLES_EQUAL(first.x, position.x, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(first.y, position.y, 1e-8);
              CPPUNIT_ASSERT(current.z <= position.z);
              current = position;
              ++iter;
              if(iter % 10 == 0)
              {
                CPPUNIT_ASSERT(tenth.z < position.z);
                tenth = position;
              }
            };
            auto printData = []( const hemelb::redblood::CellContainer & cells)
            {
              static int iter = -1;
              ++iter;
              if(cells.empty())
              {
                return;
              }
              auto cell = *cells.begin();
              auto const tag = cell->GetTag();
              auto const b = cell->GetBarycenter();
              auto const v = cell->GetVolume();
              auto const e = (*cell)();
              HEMELB_CAPTURE5(iter, tag, b, v, e);
            };
            auto printCells = [&converter](const hemelb::redblood::CellContainer &cells)
            {
              static int iter = 0;
              if(cells.empty())
              {
                return;
              }
              auto cell = *cells.begin();
              if(iter% 1000 == 0)
              {
                std::ostringstream sstr;
                sstr << "/tmp/cell-" << cell->GetTag() << "." << iter<< ".vtk";
                writeVTKMesh(sstr.str(), cell, converter);
              }
              ++iter;
            };

            int iter = 0;
            auto iterate = [&iter](const hemelb::redblood::CellContainer&)
            {
              ++iter;
            };


            CPPUNIT_ASSERT(master);
            auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
            CPPUNIT_ASSERT(controller);
            controller->AddCellChangeListener(checkVolume);
            controller->AddCellChangeListener(checkPosition);
            controller->AddCellChangeListener(iterate);
            // controller->AddCellChangeListener(printData);
            // controller->AddCellChangeListener(printCells);

            // run the simulation
            master->RunSimulation();
            master->Finalise();

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
            CPPUNIT_ASSERT(iter > 0);
          }

        private:
          std::shared_ptr<MasterSim> master;
          std::shared_ptr<hemelb::configuration::CommandLine> options;
          int const argc = 7;
          char const * argv[7];

      };



      CPPUNIT_TEST_SUITE_REGISTRATION (FadeInOutIntegrationTests);
    } // namespace redblood
  } // namespace unittests
} // namespace hemelb

#endif
