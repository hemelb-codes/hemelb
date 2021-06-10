// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_FIXTURES_H
#define HEMELB_UNITTESTS_REDBLOOD_FIXTURES_H

#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/Comparisons.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/FlowExtension.h"
#include "redblood/types.h"

#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace unittests
  {
    class BasisFixture : public CppUnit::TestFixture
    {
      public:
        void setUp()
        {
          // facets at something degrees from one another
          mesh.vertices.push_back(LatticePosition(0, 0, 0));
          mesh.vertices.push_back(LatticePosition(1, 0, 0));
          mesh.vertices.push_back(LatticePosition(0, 1, 0));
          mesh.vertices.push_back(LatticePosition(0, 0, 1));

          redblood::MeshData::Facet indices;
          indices[0] = 0;
          indices[1] = 1;
          indices[2] = 2;
          mesh.facets.push_back(indices);
          indices[0] = 0;
          indices[1] = 2;
          indices[2] = 3;
          mesh.facets.push_back(indices);
          indices[0] = 0;
          indices[1] = 3;
          indices[2] = 1;
          mesh.facets.push_back(indices);
          indices[0] = 1;
          indices[1] = 3;
          indices[2] = 2;
          mesh.facets.push_back(indices);
          redblood::orientFacets(mesh);
        }

      protected:
        hemelb::redblood::MeshData mesh;
    };

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

    //! Fake cell that contains a single node
    class NodeCell : public hemelb::redblood::Cell
    {
      public:
        NodeCell(LatticePosition const&position, std::string const &templateName = "nope") :
            NodeCell(std::vector<LatticePosition> { position }, templateName)
        {
        }
        template<class ITER>
        NodeCell(ITER first, ITER last, std::string const &templateName = "nope") :
            NodeCell(std::vector<LatticePosition> { first, last }, templateName)
        {
        }
        NodeCell(std::vector<LatticePosition> const &positions, std::string const &templateName =
                     "nope") :
                hemelb::redblood::Cell(positions,
                                           hemelb::redblood::Mesh(std::make_shared<
                                                                      hemelb::redblood::MeshData>(hemelb::redblood::MeshData { positions,
                                                                                                                               { } }),
                                                                  std::make_shared<
                                                                      hemelb::redblood::MeshTopology>()),
                                           1e0,
                                           templateName)
        {
        }

        LatticeEnergy operator()() const override
        {
          return 0e0;
        }
        LatticeEnergy operator()(std::vector<LatticeForceVector> &) const override
        {
          return 0e0;
        }
        std::unique_ptr<CellBase> cloneImpl() const override
        {
          return std::unique_ptr<NodeCell> { new NodeCell(*this) };
        }
    };
  }
}
#endif
