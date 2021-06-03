// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_REDBLOOD_FIXTURES_H
#define HEMELB_TESTS_REDBLOOD_FIXTURES_H

#include <memory>

#include "redblood/types_fwd.h"
#include "redblood/Mesh.h"

namespace hemelb {
  namespace redblood {
    class Cell;
  }

  namespace tests {
    template<class CELLTYPE = redblood::Cell>
    redblood::CellContainer TwoPancakeSamosas(LatticeDistance cutoff)
    {
      redblood::CellContainer cells;
      redblood::Mesh pancake = redblood::pancakeSamosa();
      pancake += LatticePosition(1, 1, 1) * cutoff * 0.5;
      // safer to clone so cells has its own copy
      cells.insert(std::make_shared<CELLTYPE>(pancake.clone()));
      pancake += LatticePosition(3, 0, 1) * cutoff;
      cells.insert(std::make_shared<CELLTYPE>(pancake.clone()));

      return cells;
    }

    class BasisFixture {
    public:
      BasisFixture();
    protected:
      hemelb::redblood::MeshData mesh;
    };

  }
}

#endif
