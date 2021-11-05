// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_CELLENERGY_H
#define HEMELB_REDBLOOD_CELLENERGY_H

#include <vector>

#include "units.h"
#include "redblood/Mesh.h"

namespace hemelb
{
  namespace redblood
  {
    class Facet;
    class ForceFacet;

    // Facet bending energy between two facets
    LatticeEnergy facetBending(Facet const &facetA, Facet const &facetB, Facet const &facetA_eq,
			       Facet const &facetB_eq, LatticeModulus intensity);
    // Facet bending energy and force between neighboring facets
    LatticeEnergy facetBending(ForceFacet const &facetA, ForceFacet const &facetB,
			       Facet const &facetA_eq, Facet const &facetB_eq,
			       LatticeModulus intensity);

    LatticeEnergy facetBending(MeshData::Vertices const &vertices, MeshData const &orig,
			       size_t facetIndex, size_t neighborIndex, LatticeModulus intensity,
			       std::vector<LatticeForceVector> &forces);
    LatticeEnergy facetBending(MeshData::Vertices const &vertices, MeshData const &orig,
			       size_t facetIndex, size_t neighborIndex, LatticeModulus intensity);

    LatticeEnergy volumeEnergy(MeshData::Vertices const &vertices, MeshData const &orig,
			       LatticeModulus intensity, Dimensionless origMesh_scale = 1e0);

    LatticeEnergy volumeEnergy(MeshData::Vertices const &vertices, MeshData const &orig,
			       LatticeModulus intensity, std::vector<LatticeForceVector> &forces,
			       Dimensionless origMesh_scale = 1e0);

    LatticeEnergy surfaceEnergy(MeshData::Vertices const &vertices, MeshData const &orig,
				LatticeModulus intensity, Dimensionless origMesh_scale = 1e0);
    LatticeEnergy surfaceEnergy(MeshData::Vertices const &vertices, MeshData const &orig,
				LatticeModulus intensity, std::vector<LatticeForceVector> &forces,
				Dimensionless origMesh_scale = 1e0);

    LatticeEnergy strainEnergyDensity(std::pair<Dimensionless, Dimensionless> const &strainParams,
				      LatticeModulus shearModulus, LatticeModulus dilationModulus);
    LatticeEnergy strainEnergy(Facet const &deformed, Facet const &undeformed,
			       LatticeModulus shearModulus, LatticeModulus dilationModulus,
			       Dimensionless origMesh_scale = 1e0);

    LatticeEnergy strainEnergy(ForceFacet const &deformed, Facet const &undeformed,
			       LatticeModulus shearModulus, LatticeModulus dilationModulus,
			       Dimensionless origMesh_scale = 1e0);

    LatticeEnergy strainEnergy(MeshData::Vertices const &vertices, MeshData const &origin,
			       LatticeModulus shearModulus, LatticeModulus dilationModulus,
			       Dimensionless origMesh_scale = 1e0);
    LatticeEnergy strainEnergy(MeshData::Vertices const &vertices, MeshData const &origin,
			       LatticeModulus shearModulus, LatticeModulus dilationModulus,
			       std::vector<LatticeForceVector> &forces,
			       Dimensionless origMesh_scale = 1e0);
  }
}

#endif
