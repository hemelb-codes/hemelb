// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "redblood/Borders.h"

namespace hemelb
{
  namespace ttests
  {
    using namespace redblood;

    size_t index(size_t x, size_t y, size_t z)
    {
      return (x + 1) + (y + 1) * 3 + (z + 1) * 9;
    }

    std::vector<size_t> visitme(BorderBoxIterator iterator)
    {
      std::vector<size_t> result(3 * 3 * 3, 0);
      for (; iterator; ++iterator) {
	++result[index(iterator->x, iterator->y, iterator->z)];
      }
      return result;
    }

    // Tests functionality that is *not* part of the HemeLB API
    // Checks that we know how to compute geometric properties between facets
    // However, HemeLB only cares about energy and forces
    TEST_CASE("BordersTests", "[redblood]") {
      
      SECTION("testNothingToIterate") {
	auto const visited = visitme(BorderBoxIterator(0));
	for (auto const v : visited) {
	  REQUIRE(not v);
	}
      }

      SECTION("testIterate") {
	auto testIterate = [](size_t border, std::vector<LatticeVector> const &indices) {
	  auto visited = visitme(BorderBoxIterator(border));
	  for (auto const &id : indices) {
	    REQUIRE(size_t(1) == visited[index(id.x, id.y, id.z)]);
	    visited[index(id.x, id.y, id.z)] = 0;
	  }

	  for (auto const v : visited) {
	    REQUIRE(size_t(0) == v);
	  }
	};
      
	testIterate((size_t) Borders::CENTER, { { 0, 0, 0 } });
	testIterate((size_t) Borders::BOTTOM, { { -1, 0, 0 } });
	testIterate((size_t) Borders::TOP, { { 1, 0, 0 } });
	testIterate((size_t) Borders::SOUTH, { { 0, -1, 0 } });
	testIterate((size_t) Borders::NORTH, { { 0, 1, 0 } });
	testIterate((size_t) Borders::WEST, { { 0, 0, -1 } });
	testIterate((size_t) Borders::EAST, { { 0, 0, 1 } });
	testIterate((size_t) Borders::BOTTOM + (size_t) Borders::TOP, { { -1, 0, 0 },
	      { 1, 0, 0 } });
	testIterate((size_t) Borders::BOTTOM + (size_t) Borders::CENTER,
		    { { -1, 0, 0 }, { 0, 0, 0 } });
	testIterate((size_t) Borders::BOTTOM + (size_t) Borders::NORTH, { { -1, 0, 0 }, { 0,
		1,
		0 },
											  { -1, 1, 0 } });
	testIterate((size_t) Borders::BOTTOM + (size_t) Borders::NORTH + (size_t) Borders::EAST,
		    { { -1, 0, 0 }, { 0, 1, 0 }, { -1, 1, 0 }, { 0, 0, 1 }, { -1, 0, 1 }, { 0,
										1,
										1 },
											    { -1, 1, 1 } });
      }
    }

  }
}
