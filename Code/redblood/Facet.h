// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_REDBLOOD_FACET_H
#define HEMELB_REDBLOOD_FACET_H

#include <cmath>
#include <vector>

#include "constants.h"
#include "redblood/Mesh.h"

namespace hemelb
{
  namespace redblood
  {
    // Helper class to avoid explicit indexing over vertices
    class Facet
    {
    public:
      // References nodes of a facet
      LatticePosition const *nodes[3];
      // Indices of nodes in original array
      MeshData::Facet const &indices;
      Facet(MeshData const &mesh, size_t index);
      Facet(MeshData::Vertices const &vertices, MeshData::Facets const &facets, size_t index);
      Facet(MeshData::Vertices const &vertices, MeshData::Facet const &indices);

      // returns an edge nodes[i] - nodes[j]
      LatticePosition operator()(size_t i, size_t j) const;
      // returns node i
      LatticePosition const &operator()(size_t i) const;
      // Returns edges
      LatticePosition edge(size_t i) const;
      // Returns edge length
      LatticeDistance length(size_t i) const;
      // Returns angle cosine
      Dimensionless cosine() const;
      // Returns angle sine
      Dimensionless sine() const;
      // Computes vector normal to facet
      // Order of edges determines direction of normal
      // This is taken straight from T. Krueger's code
      LatticePosition normal() const;
      // Unit vector normal to facet
      LatticePosition unitNormal() const;
      // Area of the facet
      LatticeArea area() const;
    };

    // Facet that also includes forces
    class ForceFacet : public Facet
    {
    public:
      // References forces on a node
      LatticeForceVector *forces[3];
      ForceFacet(MeshData::Vertices const &vertices, MeshData::Facet const &indices,
		 std::vector<LatticeForceVector> &forcesIn);
      ForceFacet(MeshData const &mesh, size_t index, std::vector<LatticeForceVector> &forcesIn);
      ForceFacet(MeshData::Vertices const &vertices, MeshData::Facets const &facets,
		 size_t index, std::vector<LatticeForceVector> &forcesIn);
      LatticeForceVector &GetForce(size_t i) const;
    };

    // Computes common nodes for neighboring facets
    // This routine will report nonsense if facets are not neighbors
    typedef std::pair<size_t, size_t> IndexPair;
    IndexPair commonNodes(Facet const &a, Facet const &b);

    // Figures out nodes that are not in common
    // Returns non-sense if the nodes are not neighbors.
    IndexPair singleNodes(Facet const &a, Facet const &b);

    // Computes angle between two facets
    Angle angle(LatticePosition const &a, LatticePosition const &b);
    Angle angle(Facet const &a, Facet const &b);

    Angle orientedAngle(Facet const &a, Facet const &b);

    // Returns Dxx, Dyy, Dxy packed in vector
    LatticePosition displacements(Facet const &deformed, Facet const &ref,
				  Dimensionless origMesh_scale = 1e0);

    // Returns Gxx, Gyy, Gxy packed in vector
    LatticePosition squaredDisplacements(LatticePosition const &disp);
    LatticePosition squaredDisplacements(Facet const &deformed, Facet const &ref,
					 Dimensionless origMesh_scale = 1e0);

    // Strain invariants I1 and I2
    std::pair<Dimensionless, Dimensionless> strainInvariants(LatticePosition const &squaredDisp);

    std::pair<Dimensionless, Dimensionless> strainInvariants(Facet const &deformed,
							     Facet const &ref,
							     Dimensionless origMesh_scale = 1e0);
  }
}

#endif
