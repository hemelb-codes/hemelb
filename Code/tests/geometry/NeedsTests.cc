// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <vector>

#include <catch2/catch.hpp>

#include "geometry/needs/Needs.h"

#include "tests/helpers/MockNetHelper.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::geometry;
    using namespace hemelb::tests::helpers;
    TEST_CASE_METHOD(MockNetHelper, "NeedsTests") {
      proc_t readingCores;
      site_t blockCount;
      proc_t size;
      proc_t rank;
      std::vector<bool> inputNeededBlocks;
      Needs *mockedNeeds = nullptr;

      auto SetupMocks = [&](const site_t block_count,
			    const proc_t reading_cores,
			    const proc_t core_count,
			    const proc_t current_core) {
	blockCount = block_count;
	readingCores = reading_cores;
	rank = current_core;
	size = core_count;
	setUp(core_count, current_core);

	inputNeededBlocks = std::vector<bool>(block_count);

	for (site_t i = 0; i < block_count; i++) {
	  // Mock with a tridiagonal needs example
	  inputNeededBlocks[i] = (current_core == i || current_core + 1 == i || current_core - 1 == i);
	}
      };

      auto ShareMockNeeds = [&]() {
	mockedNeeds = new Needs(blockCount,
				inputNeededBlocks,
				readingCores,
				*netMock,
				false);
      };

      SECTION("TestReadingOne") {
	SetupMocks(6, 2, 5, 0);
	REQUIRE(communicatorMock->Size() == 5);
	// Start to record the expected communications calls.
	// First will come, sending to the reading cores, each of the lengths.
	// I would expect to send to the other reading core, my count of needed cores
	int core_0_requires_from_0_count = 1;
	int core_0_requires_from_1_count = 1;

	netMock->RequireSend(&core_0_requires_from_0_count, 1, 0, "Count");
	netMock->RequireSend(&core_0_requires_from_1_count, 1, 1, "Count");

	// And I would expect the reading core to post a receive from each of the other cores,
	// asking for its count of needed blocks from this reading core
	int core_0_requires_count = 1;
	int core_1_requires_count = 2;
	int core_2_requires_count = 1;
	int core_3_requires_count = 2;
	int core_4_requires_count = 1;

	netMock->RequireReceive(&core_0_requires_count, 1, 0, "Count");
	netMock->RequireReceive(&core_1_requires_count, 1, 1, "Count");
	netMock->RequireReceive(&core_2_requires_count, 1, 2, "Count");
	netMock->RequireReceive(&core_3_requires_count, 1, 3, "Count");
	netMock->RequireReceive(&core_4_requires_count, 1, 4, "Count");

	// Then, I would expect to send to myself, my needs
	std::vector<site_t> core_0_requires_from_0;
	core_0_requires_from_0.push_back(0);
	netMock->RequireSend(&core_0_requires_from_0[0], 1, 0, "Needs");

	// Then, I would expect to send to the other reading core, my needs
	std::vector<site_t> core_0_requires_from_1;
	core_0_requires_from_1.push_back(1);
	netMock->RequireSend(&core_0_requires_from_1[0], 1, 1, "Needs");

	// Then, I should receieve from myself, my own requirements
	// Then, I would expect to receive the lists of needed blocks themselves.
	// Core 0, I expect it to need blocks 0,1
	std::vector<site_t> core_0_requires;
	core_0_requires.push_back(0);
	netMock->RequireReceive(&core_0_requires[0], 1, 0, "Needs");

	// Then, I would expect to receive the lists of needed blocks themselves.
	// From core 1, the other reading core, I expect it to need blocks 0,1,2
	std::vector<site_t> core_1_requires;
	core_1_requires.push_back(0);

	core_1_requires.push_back(2);
	netMock->RequireReceive(&core_1_requires[0], 2, 1, "Needs");
	// From core 2, I expect it to need blocks 1,2,3
	std::vector<site_t> core_2_requires;

	core_2_requires.push_back(2);

	netMock->RequireReceive(&core_2_requires[0], 1, 2, "Needs");
	// Then, I would expect to receive the lists of needed blocks themselves.
	// From core 3,  I expect it to need blocks 2,3,4
	std::vector<site_t> core_3_requires;
	core_3_requires.push_back(2);

	core_3_requires.push_back(4);
	netMock->RequireReceive(&core_3_requires[0], 2, 3, "Needs");
	// Then, I would expect to receive the lists of needed blocks themselves.
	// From core 4, I expect it to need blocks 3,4,5
	std::vector<site_t> core_4_requires;

	core_4_requires.push_back(4);

	netMock->RequireReceive(&core_4_requires[0], 1, 4, "Needs");
	ShareMockNeeds();
	// Finally, I would expect the resulting array of needs on core one to be as planned:
	std::vector<proc_t> needing_block_0;
	needing_block_0.push_back(0);
	needing_block_0.push_back(1);
	std::vector<proc_t> needing_block_1;

	std::vector<proc_t> needing_block_2;
	needing_block_2.push_back(1);
	needing_block_2.push_back(2);
	needing_block_2.push_back(3);
	std::vector<proc_t> needing_block_3;

	std::vector<proc_t> needing_block_4;
	needing_block_4.push_back(3);
	needing_block_4.push_back(4);
	std::vector<proc_t> needing_block_5;

	REQUIRE(needing_block_0 == mockedNeeds->ProcessorsNeedingBlock(0));
	REQUIRE(needing_block_1 == mockedNeeds->ProcessorsNeedingBlock(1));
	REQUIRE(needing_block_2 == mockedNeeds->ProcessorsNeedingBlock(2));
	REQUIRE(needing_block_3 == mockedNeeds->ProcessorsNeedingBlock(3));
	REQUIRE(needing_block_4 == mockedNeeds->ProcessorsNeedingBlock(4));
	REQUIRE(needing_block_5 == mockedNeeds->ProcessorsNeedingBlock(5));

	netMock->ExpectationsAllCompleted();
      }

      SECTION("TestNonReading") {
	SetupMocks(6, 2, 5, 2);

	// Start to record the expected communications calls.
	// First will come, sending to the reading cores, each of the lengths.
	// So I would expect the non-reading core to post a send to each of the reading cores, its count of needed blocks
	int core_2_requires_from_0_count = 1;
	int core_2_requires_from_1_count = 2;

	netMock->RequireSend(&core_2_requires_from_0_count, 1, 0);
	netMock->RequireSend(&core_2_requires_from_1_count, 1, 1);

	// Then, I would expect to send my list of needed blocks
	// From core 2, the other reading core, I expect it to need blocks 1,2,3
	std::vector<site_t> core_2_requires_from_0;
	std::vector<site_t> core_2_requires_from_1;

	core_2_requires_from_1.push_back(1);
	core_2_requires_from_0.push_back(2);
	core_2_requires_from_1.push_back(3);

	netMock->RequireSend(&core_2_requires_from_0[0], 1, 0);
	netMock->RequireSend(&core_2_requires_from_1[0], 2, 1);

	ShareMockNeeds();
	// Finally, I would expect the resulting array of needs to be empty
	std::vector<proc_t> empty_needs_array;

	REQUIRE(empty_needs_array == mockedNeeds->ProcessorsNeedingBlock(0));
	REQUIRE(empty_needs_array == mockedNeeds->ProcessorsNeedingBlock(1));
	REQUIRE(empty_needs_array == mockedNeeds->ProcessorsNeedingBlock(2));
	REQUIRE(empty_needs_array == mockedNeeds->ProcessorsNeedingBlock(3));
	REQUIRE(empty_needs_array == mockedNeeds->ProcessorsNeedingBlock(4));
	REQUIRE(empty_needs_array == mockedNeeds->ProcessorsNeedingBlock(5));

	netMock->ExpectationsAllCompleted();
      }

      delete mockedNeeds;
    }

  }
}
