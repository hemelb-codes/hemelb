#ifndef HEMELBSETUPTOOL_TEST_SEGMENTFACTORY_HPP
#define HEMELBSETUPTOOL_TEST_SEGMENTFACTORY_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "SegmentFactory.h"
#include "Neighbours.h"

class SegmentFactoryTests : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(SegmentFactoryTests);
  CPPUNIT_TEST(SinglePoints);
  CPPUNIT_TEST(SingleSegs);
  CPPUNIT_TEST(Four);
  CPPUNIT_TEST_SUITE_END();

  typedef SegmentFactory::Seg Seg;
public:

  void SinglePoints() {
	  // Check that for a single point we create everything correctly
	  SegmentFactory sf(1);
	  auto dirs = Neighbours::GetDisplacements();

	  for (auto dir: dirs) {
		  auto points = sf.MakeStartPoints(dir);
		  // Only one
		  CPPUNIT_ASSERT(points.size() == 1);
		  // Should start in the opposite direction
		  CPPUNIT_ASSERT(points[0] == -dir);
	  }
  }
  void SingleSegs() {
	  // Check that for a single point we create everything correctly
	  SegmentFactory sf(1);
	  auto dirs = Neighbours::GetDisplacements();

	  for (auto dir: dirs) {
		  auto segs = sf.MakeSegments(dir);
		  // Only one
		  CPPUNIT_ASSERT(segs.size() == 1);
		  // Should start in the opposite direction
		  CPPUNIT_ASSERT(segs[0].first == -dir);
		  CPPUNIT_ASSERT(segs[0].second == dir);
	  }
  }

  void Four() {
	  SegmentFactory sf(4);
	  auto dirs = Neighbours::GetDisplacements();
	  Index box_min(0);
	  Index box_max(3);

	  for (auto dir: dirs) {

		  auto segs = sf.MakeSegments(dir);

		  // can be 1,2,3 based on direction type
		  int mag = dir.GetMagnitudeSquared();
		  // Use the type to figure out how many there should be
		  switch (mag) {
		  case 1:
			  //<100> type
			  CPPUNIT_ASSERT(segs.size() == 16);
			  break;
		  case 2:
			  // <110> type
			  CPPUNIT_ASSERT(segs.size() == (16+12));
			  break;
		  case 3:
			  // <111> type
			  CPPUNIT_ASSERT(segs.size() == (16+12+9));
			  break;
		  default:
			  CPPUNIT_FAIL("Direction of unknown type!");
			  break;
		  }

		  // Check that all points between first and second are in the box
		  for (auto seg: segs) {
			  auto start = seg.first;
			  auto end = seg.second;
			  for(auto pt = start+dir; !(pt == end); pt += dir) {
				  CPPUNIT_ASSERT(pt.IsInRange(box_min, box_max));
			  }
		  }
	  }

  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(SegmentFactoryTests);

#endif
