//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_FIXTURES_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_FIXTURES_H

#include "redblood/Cell.h"
#include "redblood/parallel/CellParallelization.h"
#include "SimulationMaster.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
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
      };

      class DummyCell : public NodeCell
      {
        public:
          LatticeForce force;

          DummyCell(
              LatticePosition const&position, LatticeForce f = 0e0,
              std::string const &templateName = "nope")
            : DummyCell(std::vector<LatticePosition>{position}, f, templateName)
          {
          }
          DummyCell(
              std::vector<LatticePosition> const &positions, LatticeForce f = 0e0,
              std::string const &templateName = "nope")
            : NodeCell(positions, templateName), force(f)
          {
          }
          template<class ITER>
          DummyCell(ITER first, ITER last,
              LatticeForce f = 0e0, std::string const &templateName = "nope")
            : NodeCell(first, last, templateName), force(f)
          {
          }

          LatticeEnergy operator()(std::vector<LatticeForceVector> &f) const override
          {
            f.resize(GetNumberOfNodes());
            std::fill(f.begin(), f.end(), force);
            return 0e0;
          }
          std::unique_ptr<CellBase> cloneImpl() const override
          {
            return std::unique_ptr<DummyCell>{new DummyCell(*this)};
          }
      };

    }
  }
}

#endif  // ONCE
