// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <catch2/catch.hpp>

#include "extraction/LbDataSourceIterator.h"
#include "extraction/StraightLineGeometrySelector.h"
#include "extraction/PlaneGeometrySelector.h"
#include "extraction/WholeGeometrySelector.h"
#include "extraction/GeometrySurfaceSelector.h"
#include "extraction/SurfacePointSelector.h"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    TEST_CASE_METHOD(helpers::HasCommsTestFixture, "GeometrySelector")
    {
      
      const float VoxelSize{0.01f};
      const int CubeSize{10};
      const float CentreCoordinate{ (CubeSize - 1.0f) / 2.0f};


      auto latticeData = std::unique_ptr<FourCubeLatticeData>{FourCubeLatticeData::Create(Comms(), CubeSize + 2, 1)};

      auto simState = lb::SimulationState{60.0 / (70.0 * 5000.0), 1000};
      auto propertyCache = lb::MacroscopicPropertyCache(simState, *latticeData);
      auto unitConverter = util::UnitConverter(simState.GetTimeStepLength(),
					       VoxelSize,
					       PhysicalPosition::Zero(),
					       DEFAULT_FLUID_DENSITY_Kg_per_m3,
					       0.0);
      auto dataSourceIterator = extraction::LbDataSourceIterator(propertyCache,
								 *latticeData,
								 0,
								 unitConverter);

      // Helper to ensure all required sites are marked for extraction
      auto TestExpectedIncludedSites = [&] (const extraction::GeometrySelector& geometrySelector,
					    const std::vector<util::Vector3D<site_t>>& includedSites) {
	// For every point in the LatticeData
	dataSourceIterator.Reset();
	while (dataSourceIterator.ReadNext()) {
	  // Is the current site in includedSites?
	  bool expectedIncluded =
	    std::find(includedSites.begin(), includedSites.end(), dataSourceIterator.GetPosition()) != includedSites.end();
	  bool actualIncluded =
	    geometrySelector.Include(dataSourceIterator, dataSourceIterator.GetPosition());

	  INFO("Site at ("
	       << dataSourceIterator.GetPosition().x << ", "
	       << dataSourceIterator.GetPosition().y << ", "
	       << dataSourceIterator.GetPosition().z << ")"
	       << (expectedIncluded ? " was" : " was not")
	       << " expected to be included but actually"
	       << (actualIncluded ? " was." : " was not."))
	  REQUIRE(expectedIncluded == actualIncluded);
	}
      };
      
      auto TestOutOfGeometrySites = [&](const hemelb::extraction::GeometrySelector& geometrySelector) {
	util::Vector3D<site_t> invalidLocation{0};
	REQUIRE(!geometrySelector.Include(dataSourceIterator, invalidLocation));

	invalidLocation = util::Vector3D<site_t>{CubeSize + 1};
	REQUIRE(!geometrySelector.Include(dataSourceIterator, invalidLocation));
      };


      SECTION("StraightLineGeometrySelector") {
	const util::Vector3D<float> lineEndPoint1{CentreCoordinate * VoxelSize};
	const util::Vector3D<float> lineEndPoint2{ (CubeSize + 1) * VoxelSize};

	auto straightLineGeometrySelector = extraction::StraightLineGeometrySelector{lineEndPoint1, lineEndPoint2};

	TestOutOfGeometrySites(straightLineGeometrySelector);

	// The line runs from the centre to the (size,size,size) point.
	// This includes anything with all three coordinates the same, above
	// a minimum.
	const site_t coordMin = (CubeSize) / 2;
	
	// Gather the list of coordinates we expect to be on the line.
	std::vector<util::Vector3D<site_t> > includedCoords;

	for (site_t xCoord = coordMin; xCoord <= CubeSize; ++xCoord) {
	  includedCoords.push_back(util::Vector3D<site_t>(xCoord));
	}

	TestExpectedIncludedSites(straightLineGeometrySelector, includedCoords);
      }

      SECTION("PlaneGeometrySelector") {
	const util::Vector3D<float> planeNormal{1.0};
	const util::Vector3D<float> planePosition{CentreCoordinate * VoxelSize};
	const util::Vector3D<float> centrePoint{CentreCoordinate};
	const float planeRadius = {CubeSize * VoxelSize / 3.0f};

	auto planeGeometrySelector = extraction::PlaneGeometrySelector{planePosition, planeNormal};
	auto planeGeometrySelectorWithRadius = extraction::PlaneGeometrySelector{planePosition,
										 planeNormal,
										 planeRadius};
	TestOutOfGeometrySites(planeGeometrySelector);
	TestOutOfGeometrySites(planeGeometrySelectorWithRadius);

	// Gather the list of coordinates we expect to be on the plane with and
	// without a radius used.
	std::vector<util::Vector3D<site_t> > includedCoordsWithRadius;
	std::vector<util::Vector3D<site_t> > includedCoordsWithoutRadius;

	// The plane has normal (1, 1, 1), is centred in the cube and has a radius of
	// CubeSize / 3 lattice units (when the radius is in use).
	for (site_t xCoord = 1; xCoord <= CubeSize; ++xCoord) {
	  for (site_t yCoord = 1; yCoord <= CubeSize; ++yCoord) {
	    for (site_t zCoord = 1; zCoord <= CubeSize; ++zCoord) {
	      const auto x = util::Vector3D<site_t>(xCoord, yCoord, zCoord);
	      // Use that p.n = x.n for x on the same plane.
	      // I.e. the current point's coordinate dotted with the
	      // normal must be roughly equal to the centre point's
	      // coordinate dotted with the normal for the current
	      // point to be on the plane.
	      if (std::abs(planeNormal.Dot(x) - 3 * CentreCoordinate) > 0.5f) {
		continue;
	      }

	      includedCoordsWithoutRadius.push_back(x);

	      // Compute the distance from the centre point, include the site if it is within the radius.
	      if ( (x - centrePoint).GetMagnitudeSquared() < CubeSize*CubeSize / 9.0f) {
		includedCoordsWithRadius.push_back(x);
	      }
	    }
	  }
	}
	
	TestExpectedIncludedSites(planeGeometrySelector, includedCoordsWithoutRadius);
	TestExpectedIncludedSites(planeGeometrySelectorWithRadius, includedCoordsWithRadius);
      }

      SECTION("WholeGeometrySelector") {
	auto wholeGeometrySelector = extraction::WholeGeometrySelector{};

	TestOutOfGeometrySites(wholeGeometrySelector);

	// Variables to count the number of sites encountered, and store the iterated-over
	// positions, and properties.
	int count = 0;

	dataSourceIterator.Reset();
	while (dataSourceIterator.ReadNext()) {
	  ++count;
	  // Every site returned should be included.
	  REQUIRE(wholeGeometrySelector.Include(dataSourceIterator, dataSourceIterator.GetPosition()));
	}

	// The number of sites passed by the iterator should be the whole cube.
	REQUIRE(count == CubeSize * CubeSize * CubeSize);
      }

      SECTION("GeometrySurfaceSelector") {
	auto geometrySurfaceSelector = extraction::GeometrySurfaceSelector{};

	TestOutOfGeometrySites(geometrySurfaceSelector);

	// Gather the list of coordinates we expect to be on the geometry surface, we do not
	// consider inlets or outlets as walls.
	std::vector<util::Vector3D<site_t> > includedCoordsOnTheSurface;

	for (site_t xCoord = 1; xCoord <= CubeSize; ++xCoord) {
	  for (site_t yCoord = 1; yCoord <= CubeSize; ++yCoord) {
	    for (site_t zCoord = 1; zCoord <= CubeSize; ++zCoord) {
	      if (xCoord == 1 || xCoord == CubeSize || yCoord == 1 || yCoord == CubeSize) {
		includedCoordsOnTheSurface.emplace_back(xCoord,
							yCoord,
							zCoord);
	      }
	    }
	  }
	}
	
	TestExpectedIncludedSites(geometrySurfaceSelector, includedCoordsOnTheSurface);
      }

      SECTION("SurfacePointSelector") {
	const util::Vector3D<float> surfacePoint{(CubeSize + 0.999f) * VoxelSize};
	auto surfacePointSelector = extraction::SurfacePointSelector(surfacePoint);
	TestOutOfGeometrySites(surfacePointSelector);

	// The only site we expect to be included is the cube corner
	// (CubeSize-1, CubeSize-1, CubeSize-1) which is exactly
	// sqrt(3) times voxel size away from surfacePoint
	
	std::vector<util::Vector3D<site_t> > includedCoordsOnTheSurface;
	includedCoordsOnTheSurface.push_back(util::Vector3D<site_t>{CubeSize});

	TestExpectedIncludedSites(surfacePointSelector, includedCoordsOnTheSurface);
      }

      SECTION("SurfacePointSelectorMultipleHits") {
	// surfacePointMultipleHits is one voxel away from the cube
	// corner (CubeSize-1, CubeSize-1, CubeSize-1) in the x
	// direction, the geometry selector will pick three sites in a
	// sqrt(3) times voxel size radius.
	const util::Vector3D<float> surfacePointMultipleHits{
	  (CubeSize + 0.999f) * VoxelSize,
	  CubeSize * VoxelSize,
	  CubeSize * VoxelSize
	};

	auto surfacePointSelectorMultipleHits = extraction::SurfacePointSelector(surfacePointMultipleHits);
	TestOutOfGeometrySites(surfacePointSelectorMultipleHits);

	std::vector<util::Vector3D<site_t> > includedCoordsOnTheSurface;
	includedCoordsOnTheSurface.push_back(util::Vector3D<site_t>(CubeSize));
	includedCoordsOnTheSurface.push_back(util::Vector3D<site_t>(CubeSize,
								    CubeSize,
								    CubeSize - 1));
	includedCoordsOnTheSurface.push_back(util::Vector3D<site_t>(CubeSize,
								    CubeSize - 1,
								    CubeSize));
	includedCoordsOnTheSurface.push_back(util::Vector3D<site_t>(CubeSize,
								    CubeSize - 1,
								    CubeSize - 1));

	TestExpectedIncludedSites(surfacePointSelectorMultipleHits, includedCoordsOnTheSurface);
      }

    }
  }
}
