#ifndef HEMELB_UNITTESTS_REDBLOOD_FADEINOUTINTEGRATION_H
#define HEMELB_UNITTESTS_REDBLOOD_FADEINOUTINTEGRATION_H

#include <cppunit/extensions/HelperMacros.h>
#include <memory>

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
        CPPUNIT_TEST_SUITE(FadeInOutIntegrationTests);
          CPPUNIT_TEST(testIntegration);CPPUNIT_TEST_SUITE_END()
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
            ModifyXMLInput("large_cylinder_rbc.xml", std::move(intel), 100);

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

          void tearDown()
          {
            master->Finalise();
            master.reset();
          }

          void testIntegration()
          {
            // add callback to put cell positions in a vector
            std::function<void(const hemelb::redblood::CellContainer &)> output_callback =
                [](const hemelb::redblood::CellContainer & cells)
                {
                  for (auto cell: cells)
                    std::cout << "hemelb::redblood::Cell@" << std::addressof(cell) << ": " << cell->GetBarycenter() << std::endl;
                };
            std::shared_ptr<CellControl> controller =
                std::static_pointer_cast<CellControl>(master->GetCellController());
            controller->AddCellChangeListener(output_callback);

            // run the simulation
            master->RunSimulation();
          }

        private:
          std::shared_ptr<MasterSim> master;
          std::shared_ptr<hemelb::configuration::CommandLine> options;
          int const argc = 7;
          char const * argv[7];

      };
      // class FadeInOutIntegrationTests

      CPPUNIT_TEST_SUITE_REGISTRATION(FadeInOutIntegrationTests);

    } // namespace redblood
  } // namespace unittests
} // namespace hemelb

#endif
