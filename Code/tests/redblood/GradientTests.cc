// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "resources/Resource.h"
#include "redblood/CellEnergy.h"
#include "tests/redblood/Fixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    struct EnergyFunctional {
      EnergyFunctional(MeshData const &original) :
	original(original)
      {
      }
      MeshData const &original;
    };

#define HEMELB_OP_MACRO(CLASS, FUNCTION)				\
  struct CLASS : public EnergyFunctional                                                           \
  {                                                                                                \
    CLASS(MeshData const &original) : EnergyFunctional(original)                                   \
    {                                                                                              \
    }                                                                                              \
    LatticeEnergy operator()(MeshData const &mesh) const                                           \
    {                                                                                              \
      return FUNCTION(mesh.vertices, original, 1.0);                                               \
    }                                                                                              \
    LatticeEnergy operator()(MeshData const &mesh, std::vector<LatticeForceVector> &forces) const  \
    {                                                                                              \
      return FUNCTION(mesh.vertices, original, 1.0, forces);                                       \
    }                                                                                              \
  }
    HEMELB_OP_MACRO(VolumeEnergy, volumeEnergy);
    HEMELB_OP_MACRO(SurfaceEnergy, surfaceEnergy);
#undef HEMELB_OP_MACRO

    struct BendingEnergy : public EnergyFunctional
    {
      BendingEnergy(MeshData const &original) :
	EnergyFunctional(original)
      {
      }
      LatticeEnergy operator()(MeshData const &mesh) const
      {
	return facetBending(mesh.vertices, original, 0, 1, 1.0)
	  + facetBending(mesh.vertices, original, 0, 2, 1.0)
	  + facetBending(mesh.vertices, original, 0, 3, 1.0)
	  + facetBending(mesh.vertices, original, 1, 2, 1.0)
	  + facetBending(mesh.vertices, original, 1, 3, 1.0)
	  + facetBending(mesh.vertices, original, 2, 3, 1.0);
      }
      LatticeEnergy operator()(MeshData const &mesh,
			       std::vector<LatticeForceVector> &forces) const
      {
	return facetBending(mesh.vertices, original, 0, 1, 1.0, forces)
	  + facetBending(mesh.vertices, original, 0, 2, 1.0, forces)
	  + facetBending(mesh.vertices, original, 0, 3, 1.0, forces)
	  + facetBending(mesh.vertices, original, 1, 2, 1.0, forces)
	  + facetBending(mesh.vertices, original, 1, 3, 1.0, forces)
	  + facetBending(mesh.vertices, original, 2, 3, 1.0, forces);
      }
    };

    struct StrainEnergy : public EnergyFunctional
    {
      StrainEnergy(MeshData const &original) :
	EnergyFunctional(original)
      {
      }
      LatticeEnergy operator()(MeshData const &mesh) const
      {
	return strainEnergy(mesh.vertices, original, 1.0, 2.0);
      }
      LatticeEnergy operator()(MeshData const &mesh,
			       std::vector<LatticeForceVector> &forces) const
      {
	return strainEnergy(mesh.vertices, original, 1.0, 2.0, forces);
      }
    };

    TEST_CASE_METHOD(EnergyVsGradientFixture, "GradientTests", "[redblood]") {
      MeshData original = mesh;
      mesh.vertices[0] += LatticePosition(-0.01, 0.02342, 0.03564);
      mesh.vertices[1] += LatticePosition(0.0837, -0.012632, 0.0872935);
      mesh.vertices[2] += LatticePosition(0.02631, -0.00824223, -0.098362);
    
      SECTION("testBending") {
	energyVsForces(BendingEnergy(original));
      }
      SECTION("testSurface") {
	energyVsForces(SurfaceEnergy(original));
      }
      SECTION("testVolume") {
	energyVsForces(VolumeEnergy(original));
      }
      SECTION("testStrain") {
	energyVsForces(StrainEnergy(original));
      }
    }
  }
}
