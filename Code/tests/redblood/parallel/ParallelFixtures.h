// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_REDBLOOD_PARALLEL_PARALLELFIXTURES_H
#define HEMELB_TESTS_REDBLOOD_PARALLEL_PARALLELFIXTURES_H

#include "redblood/Cell.h"
#include "tests/redblood/Fixtures.h"
#include "SimulationMaster.h"

namespace hemelb
{
  namespace tests
  {
    //! \brief gathers mid-domain and egde positions from all procs
    //! \details If there are insufficient number of edges, mid-domains are used instead.
    //! erase removes the components from the first process.
    std::vector<LatticeVector> GatherSpecialPositions(geometry::Domain const & domain,
						      size_t mid, size_t edges,
						      net::MpiCommunicator const &c);

    //! Make some functionality available
    template<class TRAITS>
    class OpenedSimulationMaster : public SimulationMaster<TRAITS>
    {
    public:
      using SimulationMaster<TRAITS>::SimulationMaster;
      using SimulationMaster<TRAITS>::Finalise;

      void DoTimeStep()
      {
	SimulationMaster<TRAITS>::DoTimeStep();
      }

      std::shared_ptr<hemelb::configuration::SimConfig> GetSimConfig()
      {
	return this->simConfig;
      }
    };

    class DummyCell : public NodeCell {
    public:
      LatticeForce force;

      DummyCell(LatticePosition const&position, LatticeForce f = 0e0,
		std::string const &templateName = "nope");
      DummyCell(std::vector<LatticePosition> const &positions, LatticeForce f = 0e0,
		std::string const &templateName = "nope");
      template<class ITER>
      DummyCell(ITER first, ITER last, LatticeForce f = 0e0, std::string const &templateName =
		"nope") :
	NodeCell(first, last, templateName), force(f)
      {
      }

      LatticeEnergy operator()(std::vector<LatticeForceVector> &f) const override;
      std::unique_ptr<CellBase> cloneImpl() const override;
    };

    //! Creates list of cells for each set of positions from each process
    std::vector<hemelb::redblood::CellContainer::value_type> CreateCellsFromSpecialPositions(
											     geometry::Domain const & domain, size_t mid, size_t edges,
											     net::MpiCommunicator const &c, size_t nCells);

    net::MpiCommunicator CreateDumbGraphComm(net::MpiCommunicator const &comm);
  }
}

#endif  // ONCE
