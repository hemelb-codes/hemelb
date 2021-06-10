// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/redblood/Fixtures.h"
namespace hemelb {
  namespace tests {
    
    BasisFixture::BasisFixture() {
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

    SquareDuctTetrahedronFixture::SquareDuctTetrahedronFixture(redblood::Mesh const & initial_mesh, size_t refinement) :
      FourCubeBasedTestFixture{}, mesh{redblood::refine(initial_mesh, refinement)} {
      mesh *= Dimensionless(cubeSizeWithHalo - 3) * 0.5;
      mesh += LatticePosition{Dimensionless(cubeSizeWithHalo) * 0.5};
    }

    FlowExtensionFixture::FlowExtensionFixture() {
      flowExt.normal = util::Vector3D<LatticeDistance>(1.0, 0.0, 0.0);
      flowExt.origin = util::Vector3D<LatticeDistance>(0.0, 0.0, 0.0);
      flowExt.length = 10.0;
      flowExt.radius = 1.0;
    }

  }
}
