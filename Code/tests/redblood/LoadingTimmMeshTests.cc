// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <sstream>

#include <catch2/catch.hpp>

#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"
#include "tests/redblood/Fixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    // Return number of flips
    unsigned checkMeshOrientation(MeshData const &mesh)
    {
      MeshData copy(mesh);
      return orientFacets(copy);
    }

    std::shared_ptr<MeshData> readMeshFromFileNoFix(const std::string filename) {
      auto fullPath = resources::Resource(filename).Path();
      REQUIRE(std::filesystem::exists(fullPath));
      return redblood::KruegerMeshIO{}.readFile(fullPath, false);
    }
    
    TEST_CASE_METHOD(BasisFixture,
                     "Test we can load .msh files output by Timm's code",
                     "[redblood]")
    {

      SECTION("A MSH file known to have inconsistent face orientations needs some faces flipped") {
        // This is one of Timm's original ico files and has a face
        // with the wrong vertex orientation leading to an inward
        // pointing normal. Loading it without fixing face orientation
        // should lead to inconsistent orientation
        auto mesh = readMeshFromFileNoFix("rbc_ico_720.msh");
        REQUIRE(checkMeshOrientation(*mesh) > 0);
      }

      SECTION("The mesh file with faces corrected has no flips") {
        // This is the same ico file with the face fixed
        // Loading it without fixing face orientation should be fine
        auto mesh = readMeshFromFileNoFix("rbc_ico_720_correct.msh");
        REQUIRE(checkMeshOrientation(*mesh) == 0);
      }

      SECTION("An output mesh from Timm's code has faces oriented consistently (all inward for this case)") {
        // This is a mesh taken from a simulation in Timm's code,
	// which should have correct face orientation by definition
	auto mesh = readMeshFromFileNoFix("992Particles_rank3_26_t992.msh");
	REQUIRE(checkMeshOrientation(*mesh) );
      }
    }    
  }
}

