
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_GEOMETRY_NEEDSTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEEDSTESTS_H
#include "geometry/needs/Needs.h"
#include "unittests/helpers/MockNetHelper.h"
#include "unittests/helpers/CppUnitCompareVectors.h"
#include <cppunit/TestFixture.h>

namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {
      using namespace hemelb::geometry;
      using namespace hemelb::unittests::helpers;
      class NeedsTests : public CppUnit::TestFixture, MockNetHelper
      {
          CPPUNIT_TEST_SUITE (NeedsTests);
          CPPUNIT_TEST (TestReadingOne);
          CPPUNIT_TEST (TestNonReading);CPPUNIT_TEST_SUITE_END();

        public:
          NeedsTests() :
              mockedNeeds(NULL)
          {
          }

          void setUp()
          {
          }

          void tearDown()
          {
            delete mockedNeeds;
            MockNetHelper::tearDown();
          }

          void TestReadingOne()
          {
            SetupMocks(6, 2, 5, 0);
            CPPUNIT_ASSERT_EQUAL(communicatorMock->Size(),5);
            // Start to record the expected communications calls.
	    // First we gather the counts of blocks needed on each reading core
	    auto mock = std::dynamic_pointer_cast<MockMpiCommunicator>(communicatorMock);
	    const std::vector<int> block_size_reading_core_0 = {1, 2, 1, 2, 1};
	    // We don't really care about what the other ranks need from RC1
	    const std::vector<int> block_size_reading_core_1 = {1, 0, 2, 0, 0};
	    
	    mock->AddGatherResult(block_size_reading_core_0);
	    mock->AddGatherResult(block_size_reading_core_1);
	    
	    // Now we list the blocks needed from reading core 0 (us)
	    // Note counts match above
	    const std::vector<site_t> blocks_needed_from_reading_core_0 = {
	      0,
	      0, 2,
	      2,
	      2, 4,
	      4
	    };
	    mock->AddGatherVResult(blocks_needed_from_reading_core_0);
	    // Blocks needed from RC1
	    const std::vector<site_t> blocks_needed_from_reading_core_1 = {1};
	    mock->AddGatherVResult(blocks_needed_from_reading_core_1);
	    
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

            CPPUNIT_ASSERT_EQUAL(needing_block_0, mockedNeeds->ProcessorsNeedingBlock(0));
            CPPUNIT_ASSERT_EQUAL(needing_block_1, mockedNeeds->ProcessorsNeedingBlock(1));
            CPPUNIT_ASSERT_EQUAL(needing_block_2, mockedNeeds->ProcessorsNeedingBlock(2));
            CPPUNIT_ASSERT_EQUAL(needing_block_3, mockedNeeds->ProcessorsNeedingBlock(3));
            CPPUNIT_ASSERT_EQUAL(needing_block_4, mockedNeeds->ProcessorsNeedingBlock(4));
            CPPUNIT_ASSERT_EQUAL(needing_block_5, mockedNeeds->ProcessorsNeedingBlock(5));
	    // I guess this means that no comms went through the net
            netMock->ExpectationsAllCompleted();
          }

          void TestNonReading()
          {
	    // We are core 2 - a no reading one
            SetupMocks(6, 2, 5, 2);
	    
	    CPPUNIT_ASSERT_EQUAL(communicatorMock->Size(),5);

	    // Record the expected communications calls.
	    auto mock = std::dynamic_pointer_cast<MockMpiCommunicator>(communicatorMock);
	    
	    // First we gather the counts of blocks needed on each reading core
	    // For RC0 we know this from above
	    const std::vector<int> block_size_reading_core_0 = {1, 2, 1, 2, 1};
	    // We don't really care about what the other ranks need from RC1
	    const std::vector<int> block_size_reading_core_1 = {0, 0, 2, 0, 0};
	    
	    mock->AddGatherResult(block_size_reading_core_0);
	    mock->AddGatherResult(block_size_reading_core_1);
	    
	    // Now we list the blocks needed from RC0
	    // Note counts match above
	    const std::vector<site_t> blocks_needed_from_reading_core_0 = {
	      0,
	      0, 2,
	      2,
	      2, 4,
	      4
	    };
	    mock->AddGatherVResult(blocks_needed_from_reading_core_0);
	    // Blocks needed from RC1
	    const std::vector<site_t> blocks_needed_from_reading_core_1 = {1 ,3};
	    mock->AddGatherVResult(blocks_needed_from_reading_core_1);

            ShareMockNeeds();
            // Finally, I would expect the resulting array of needs to be empty
            std::vector<proc_t> empty_needs_array;

            CPPUNIT_ASSERT_EQUAL(empty_needs_array, mockedNeeds->ProcessorsNeedingBlock(0));
            CPPUNIT_ASSERT_EQUAL(empty_needs_array, mockedNeeds->ProcessorsNeedingBlock(1));
            CPPUNIT_ASSERT_EQUAL(empty_needs_array, mockedNeeds->ProcessorsNeedingBlock(2));
            CPPUNIT_ASSERT_EQUAL(empty_needs_array, mockedNeeds->ProcessorsNeedingBlock(3));
            CPPUNIT_ASSERT_EQUAL(empty_needs_array, mockedNeeds->ProcessorsNeedingBlock(4));
            CPPUNIT_ASSERT_EQUAL(empty_needs_array, mockedNeeds->ProcessorsNeedingBlock(5));

            netMock->ExpectationsAllCompleted();
          }

          void SetupMocks(const site_t block_count,
                          const proc_t reading_cores,
                          const proc_t core_count,
                          const proc_t current_core)
          {
            blockCount = block_count;
            readingCores = reading_cores;
            rank = current_core;
            size = core_count;
            MockNetHelper::setUp(core_count,current_core);

            inputNeededBlocks = std::vector<bool>(block_count);

            for (site_t i = 0; i < block_count; i++)
            {
              // Mock with a tridiagonal needs example
              inputNeededBlocks[i] = (current_core == i || current_core + 1 == i || current_core - 1 == i);
            }
          }

          void ShareMockNeeds()
          {
            mockedNeeds = new Needs(blockCount,
                                          inputNeededBlocks,
                                          readingCores,
				          communicatorMock,
                                          false);
          }

        private:
          proc_t readingCores;
          site_t blockCount;
          proc_t size;
          proc_t rank;
          std::vector<bool> inputNeededBlocks;
          Needs *mockedNeeds;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (NeedsTests);
    }
  }
}
#endif // ONCE
