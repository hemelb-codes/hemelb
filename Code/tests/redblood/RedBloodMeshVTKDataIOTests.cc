// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <sstream>

#include <catch2/catch.hpp>
#include <vtkPolyData.h>

#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"
#include "resources/Resource.h"
#include "util/UnitConverter.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    bool any(std::array<redblood::IdType, 3> const& v, redblood::IdType const value) {
      return std::find(v.begin(), v.end(), value) != v.end();
    }

    TEST_CASE("Test reading and writing RBC meshes in VTK format", "[redblood]") {

      redblood::VTKMeshIO io = {};

      SECTION("testVTPReadMesh") {
	std::string filename = resources::Resource("rbc_ico_720.vtp").Path();
	std::shared_ptr<MeshData> mesh = io.readFile(filename, true);

	REQUIRE(mesh);
	REQUIRE(mesh->vertices.size() == 362ul);
	REQUIRE(mesh->facets.size() == 720ul);

	using Vector3d = util::Vector3D<double>;
	Vector3d vfirst(0.850650808352040, 0.525731112119134, 0.000000000000000);
	Vector3d vlast(0.289241015911102, -0.180331975221633, -0.199399739503860);
	REQUIRE(mesh->vertices.front() == ApproxV(vfirst));
	REQUIRE(mesh->vertices.back() == ApproxV(vlast));

	REQUIRE(any(mesh->facets.front(), 0));
	REQUIRE(any(mesh->facets.front(), 12));
	REQUIRE(any(mesh->facets.front(), 17));
	REQUIRE(any(mesh->facets.back(), 67));
	REQUIRE(any(mesh->facets.back(), 361));
	REQUIRE(any(mesh->facets.back(), 157));
      }

      SECTION("testOrientOriginalRBCMesh") {
	// This file has 684 out of 720 faces oriented inwards
	std::string filename = resources::Resource("rbc_ico_720.vtp").Path();

	std::shared_ptr<MeshData> meshData;
	vtkSmartPointer<vtkPolyData> polyData;
	std::tie(meshData, polyData) = io.readUnoriented(redblood::VTKMeshIO::Mode::file, filename);

	auto numSwaps = orientFacets(*meshData, *polyData);

	REQUIRE(numSwaps == 684u);
	REQUIRE(volume(*meshData) > 0.0);
      }

      SECTION("testOrientTimmSimRBCMesh") {
	// This file has all 720 faces oriented inwards
	std::string filename = resources::Resource("992Particles_rank3_26_t992.vtp").Path();

	std::shared_ptr<MeshData> meshData;
	vtkSmartPointer<vtkPolyData> polyData;
	std::tie(meshData, polyData) = io.readUnoriented(redblood::VTKMeshIO::Mode::file, filename);

	auto numSwaps = orientFacets(*meshData, *polyData);

	REQUIRE(numSwaps == 0u);
	REQUIRE(volume(*meshData) > 0.0);
      }
    }
  }
}
