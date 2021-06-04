// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <sstream>
#include <catch2/catch.hpp>

#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"
#include "resources/Resource.h"
//#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    TEST_CASE("TopologyTests") {
      std::string filename = resources::Resource("red_blood_cube.txt").Path();
      std::shared_ptr<MeshData> mesh = redblood::KruegerMeshIO{}.readFile(filename, true);
      // Checks the mesh input makes sense
      REQUIRE(mesh->vertices.size() == 8);
      REQUIRE(mesh->facets.size() == 12);

      // Creates topology
      std::shared_ptr<MeshTopology> topo = std::make_shared<MeshTopology>(*mesh);
      
      SECTION("testNodeToVertex") {
	REQUIRE(topo->vertexToFacets.size() == 8);

	// expected[vertex] = {nfacets, facet indices}
	unsigned int expected[8][6] = { { 4, 0, 1, 6, 9 }, { 5, 0, 2, 3, 8, 9 }, { 4,
										   0,
										   1,
										   3,
										   10 },
					{ 5, 1, 6, 7, 10, 11 }, { 5, 4, 6, 7, 8, 9 }, { 4,
											2,
											4,
											5,
											8 },
					{ 5, 2, 3, 5, 10, 11 }, { 4, 4, 5, 7, 11 }, };

	for (unsigned vertex(0); vertex < 8; ++vertex)
	  {
	    std::set<size_t> facets = topo->vertexToFacets[vertex];
	    REQUIRE(facets.size() == expected[vertex][0]);
	    
	    for (size_t facet(1); facet <= expected[vertex][0]; ++facet)
              {
                REQUIRE(facets.count(expected[vertex][facet]) == 1);
              }
	  }
      }
      
      SECTION("testFacetNeighbors") {
	REQUIRE(topo->facetNeighbors.size() == 12);

	// expected[facet] = {neighbor indices}
	redblood::IdType expected[12][3] = { { 1, 3, 9 }, { 0, 6, 10 }, { 3, 5, 8 }, { 0, 2, 10 }, { 5,
												     7,
												     8 },
					     { 4, 2, 11 }, { 1, 7, 9 }, { 6, 4, 11 }, { 9, 4, 2 }, { 0,
												     6,
												     8 },
					     { 1, 3, 11 }, { 10, 5, 7 }, };

	for (unsigned facet(0); facet < 12; ++facet)
	  {
	    auto const& neighs = topo->facetNeighbors[facet];
	    REQUIRE(neighs.size() == 3);

	    REQUIRE(std::is_permutation(neighs.begin(), neighs.end(), expected[facet]));
	  }
      }
    }
  }
}
