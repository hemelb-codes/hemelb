//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_FIXTURE_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_FIXTURE_TESTS_H

#include <cppunit/TestFixture.h>

#include <algorithm>
#include <random>
#include <memory>
#include <iterator>

#include "redblood/parallel/IntegrateVelocities.h"
#include "redblood/parallel/CellParallelization.h"
#include "redblood/parallel/NodeCharacterizer.h"
#include "configuration/CommandLine.h"
#include "SimulationMaster.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/LatticeDataAccess.h"
#include "unittests/helpers/FolderTestFixture.h"
#include "unittests/redblood/parallel/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      //! Parallel and Sequential move in lock step
      class ParallelFixtureTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (ParallelFixtureTests);
          CPPUNIT_TEST (testMidDomainNodes);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();

          void testMidDomainNodes();
        protected:
          std::shared_ptr<hemelb::configuration::CommandLine> options;

          //! Meta-function to create simulation type
          template<class STENCIL>
          struct MasterSim
          {
              typedef ::hemelb::Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type LBTraits;
              typedef typename LBTraits::ChangeStencil<STENCIL>::Type Traits;
              typedef OpenedSimulationMaster<Traits> Type;
          };

          //! Creates a master simulation
          template<class STENCIL>
          std::shared_ptr<typename MasterSim<STENCIL>::Type> CreateMasterSim(
              net::MpiCommunicator const &comm) const
          {
            typedef typename MasterSim<STENCIL>::Type MasterSim;
            return std::make_shared<MasterSim>(*options, comm);
          }
      };


      void ParallelFixtureTests::setUp()
      {
         FolderTestFixture::setUp();

         // Have everything ready to creates simulations
         if (Comms().Rank() == 0)
         {
           CopyResourceToTempdir("large_cylinder_rbc.xml");
           CopyResourceToTempdir("large_cylinder.gmy");
           CopyResourceToTempdir("red_blood_cell.txt");

           ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 20000);
         }
         Comms().Barrier();

         options = std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
                                                               "-in",
                                                               "large_cylinder_rbc.xml",
                                                               "-i",
                                                               "1",
                                                               "-ss",
                                                               "1111",
                                                               "-out",
                                                               "" });
      }

      void ParallelFixtureTests::testMidDomainNodes()
      {
        typedef hemelb::redblood::stencil::TwoPoint Stencil;

        auto const world = Comms();
        if(world.Size() < 2)
        {
          return;
        }

        auto master = CreateMasterSim<Stencil>(world);
        CPPUNIT_ASSERT(master);

        auto &latDat = master->GetLatticeData();
        helpers::ZeroOutForces(latDat);
        auto const positions = GatherSpecialPositions(latDat, 3, 0, world);

        CPPUNIT_ASSERT_EQUAL(std::size_t(3 * (world.Size() - 1)), positions.size());
        for(std::size_t i(0); i < positions.size(); ++i)
        {
          auto const procs
            = hemelb::redblood::parallel::details::positionAffectsProcs<Stencil>(latDat, positions[i]);
          std::ostringstream sstr;
          sstr << world.Rank() << " - " << i << " - " << positions[i] << ": " ;
          for(auto proc: procs)
            sstr << proc << " ";
          std::cout << sstr.str() << std::endl;
          // CPPUNIT_ASSERT(procs.size() >= 1);
          // CPPUNIT_ASSERT_EQUAL(std::size_t(1), procs.count(proc_t(i / 3 + 1)));
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (ParallelFixtureTests);
    }
  }
}

#endif  // ONCE
