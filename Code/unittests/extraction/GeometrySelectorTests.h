
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_EXTRACTION_GEOMETRYSELECTORTESTS_H
#define HEMELB_UNITTESTS_EXTRACTION_GEOMETRYSELECTORTESTS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <cppunit/TestFixture.h>
#include "extraction/LbDataSourceIterator.h"
#include "extraction/StraightLineGeometrySelector.h"
#include "extraction/PlaneGeometrySelector.h"
#include "extraction/WholeGeometrySelector.h"
#include "extraction/GeometrySurfaceSelector.h"
#include "unittests/FourCubeLatticeData.h"
#include "unittests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace extraction
    {
      class GeometrySelectorTests : public helpers::HasCommsTestFixture
      {
          CPPUNIT_TEST_SUITE ( GeometrySelectorTests);
          CPPUNIT_TEST ( TestStraightLineGeometrySelector);
          CPPUNIT_TEST ( TestPlaneGeometrySelector);
          CPPUNIT_TEST ( TestWholeGeometrySelector);
          CPPUNIT_TEST ( TestGeometrySurfaceSelector);
          CPPUNIT_TEST ( TestSurfacePointSelector);
          CPPUNIT_TEST ( TestSurfacePointSelectorMultipleHits);CPPUNIT_TEST_SUITE_END();

        public:
          GeometrySelectorTests() :
                VoxelSize(0.01),
                CubeSize(10),
                CentreCoordinate( ((distribn_t) CubeSize - 1.0) / 2.0),
                planeNormal(1.0),

                planePosition(CentreCoordinate * VoxelSize),
                planeRadius(distribn_t(CubeSize) * VoxelSize / 3.0),
                lineEndPoint1(CentreCoordinate * VoxelSize),
                lineEndPoint2( (CubeSize + 1) * VoxelSize),
                surfacePoint( (CubeSize + 1) * VoxelSize),

                surfacePointMultipleHits( (CubeSize + 1) * VoxelSize,
                                         CubeSize * VoxelSize,
                                         CubeSize * VoxelSize)
          {

          }

          void setUp()
          {
            helpers::HasCommsTestFixture::setUp();
            latticeData = FourCubeLatticeData::Create(Comms(), CubeSize + 2, 1);
            simState = new lb::SimulationState(60.0 / (70.0 * 5000.0), 1000);
            propertyCache = new hemelb::lb::MacroscopicPropertyCache(*simState, *latticeData);
            unitConverter = new hemelb::util::UnitConverter(simState->GetTimeStepLength(),
                                                            VoxelSize,
                                                            PhysicalPosition::Zero());
            dataSourceIterator = new hemelb::extraction::LbDataSourceIterator(*propertyCache,
                                                                              *latticeData,
                                                                              0,
                                                                              *unitConverter);

            planeGeometrySelector = new hemelb::extraction::PlaneGeometrySelector(planePosition,
                                                                                  planeNormal);
            planeGeometrySelectorWithRadius
                = new hemelb::extraction::PlaneGeometrySelector(planePosition,
                                                                planeNormal,
                                                                planeRadius);
            straightLineGeometrySelector
                = new hemelb::extraction::StraightLineGeometrySelector(lineEndPoint1, lineEndPoint2);
            wholeGeometrySelector = new hemelb::extraction::WholeGeometrySelector();

            geometrySurfaceSelector = new hemelb::extraction::GeometrySurfaceSelector();

            surfacePointSelector = new hemelb::extraction::SurfacePointSelector(surfacePoint);
            surfacePointSelectorMultipleHits
                = new hemelb::extraction::SurfacePointSelector(surfacePointMultipleHits);
          }

          void tearDown()
          {
            delete planeGeometrySelector;
            delete planeGeometrySelectorWithRadius;
            delete straightLineGeometrySelector;
            delete wholeGeometrySelector;
            delete geometrySurfaceSelector;

            delete dataSourceIterator;
            delete unitConverter;
            delete propertyCache;
            delete simState;
            delete latticeData;
            helpers::HasCommsTestFixture::tearDown();
          }

          void TestStraightLineGeometrySelector()
          {
            TestOutOfGeometrySites(straightLineGeometrySelector);

            // The line runs from the centre to the (size,size,size) point.
            // This includes anything with all three coordinates the same, above
            // a minimum.
            const site_t coordMin = (CubeSize) / 2;

            // Gather the list of coordinates we expect to be on the line.
            std::vector<util::Vector3D<site_t> > includedCoords;

            for (site_t xCoord = coordMin; xCoord <= CubeSize; ++xCoord)
            {
              includedCoords.push_back(util::Vector3D<site_t>(xCoord));
            }

            TestExpectedIncludedSites(straightLineGeometrySelector, includedCoords);
          }

          void TestPlaneGeometrySelector()
          {
            TestOutOfGeometrySites(planeGeometrySelector);
            TestOutOfGeometrySites(planeGeometrySelectorWithRadius);

            // Gather the list of coordinates we expect to be on the plane with and
            // without a radius used.
            std::vector<util::Vector3D<site_t> > includedCoordsWithRadius;
            std::vector<util::Vector3D<site_t> > includedCoordsWithoutRadius;

            // The plane has normal (1, 1, 1), is centred in the cube and has a radius of
            // CubeSize / 3 lattice units (when the radius is in use).
            for (site_t xCoord = 1; xCoord <= CubeSize; ++xCoord)
            {
              for (site_t yCoord = 1; yCoord <= CubeSize; ++yCoord)
              {
                for (site_t zCoord = 1; zCoord <= CubeSize; ++zCoord)
                {
                  // Use that p.n = x.n for x on the same plane.
                  // I.e. the current point's coordinate dotted with the normal must be roughly equal
                  // to the centre point's coordinate dotted with the normal for the current point
                  // to be on the plane.
                  if (std::abs(distribn_t(xCoord + yCoord + zCoord) - 3 * CentreCoordinate) > 0.5)
                  {
                    continue;
                  }

                  includedCoordsWithoutRadius.push_back(util::Vector3D<site_t>(xCoord,
                                                                               yCoord,
                                                                               zCoord));

                  // Compute the distance from the centre point, include the site if it is within the radius.
                  if ( (util::Vector3D<distribn_t>(xCoord, yCoord, zCoord) - util::Vector3D<
                      distribn_t>(CentreCoordinate)).GetMagnitude() < distribn_t(CubeSize) / 3.0)
                  {
                    includedCoordsWithRadius.push_back(util::Vector3D<site_t>(xCoord,
                                                                              yCoord,
                                                                              zCoord));
                  }
                }
              }
            }

            TestExpectedIncludedSites(planeGeometrySelector, includedCoordsWithoutRadius);
            TestExpectedIncludedSites(planeGeometrySelectorWithRadius, includedCoordsWithRadius);
          }

          void TestWholeGeometrySelector()
          {
            TestOutOfGeometrySites(wholeGeometrySelector);

            // Variables to count the number of sites encountered, and store the iterated-over
            // positions, and properties.
            int count = 0;

            dataSourceIterator->Reset();
            while (dataSourceIterator->ReadNext())
            {
              ++count;
              // Every site returned should be included.
              CPPUNIT_ASSERT(wholeGeometrySelector->Include(*dataSourceIterator,
                                                            dataSourceIterator->GetPosition()));
            }

            // The number of sites passed by the iterator should be the whole cube.
            CPPUNIT_ASSERT_EQUAL(CubeSize * CubeSize * CubeSize, count);
          }

          void TestGeometrySurfaceSelector()
          {
            TestOutOfGeometrySites(geometrySurfaceSelector);

            // Gather the list of coordinates we expect to be on the geometry surface, we do not
            // consider inlets or outlets as walls.
            std::vector<util::Vector3D<site_t> > includedCoordsOnTheSurface;

            for (site_t xCoord = 1; xCoord <= CubeSize; ++xCoord)
            {
              for (site_t yCoord = 1; yCoord <= CubeSize; ++yCoord)
              {
                for (site_t zCoord = 1; zCoord <= CubeSize; ++zCoord)
                {
                  if (xCoord == 1 || xCoord == CubeSize || yCoord == 1 || yCoord == CubeSize)
                  {
                    includedCoordsOnTheSurface.push_back(util::Vector3D<site_t>(xCoord,
                                                                                yCoord,
                                                                                zCoord));
                  }
                }
              }
            }

            TestExpectedIncludedSites(geometrySurfaceSelector, includedCoordsOnTheSurface);
          }

          void TestSurfacePointSelector()
          {
            TestOutOfGeometrySites(surfacePointSelector);

            // The only site we expect to be included is the cube corner (CubeSize-1, CubeSize-1, CubeSize-1) which is
            // exactly sqrt(3) times voxel size away from surfacePoint
            std::vector<util::Vector3D<site_t> > includedCoordsOnTheSurface;
            includedCoordsOnTheSurface.push_back(util::Vector3D<site_t>(CubeSize));

            TestExpectedIncludedSites(surfacePointSelector, includedCoordsOnTheSurface);
          }

          void TestSurfacePointSelectorMultipleHits()
          {
            TestOutOfGeometrySites(surfacePointSelectorMultipleHits);

            // surfacePointMultipleHits is one voxel away from the cube corner (CubeSize-1, CubeSize-1, CubeSize-1) in
            // the x direction, the geometry selector will pick three sites in a sqrt(3) times voxel size radius.
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

        private:
          void TestOutOfGeometrySites(hemelb::extraction::GeometrySelector* geometrySelector)
          {
            util::Vector3D<site_t> invalidLocation(0);
            CPPUNIT_ASSERT(!geometrySelector->Include(*dataSourceIterator, invalidLocation));

            invalidLocation = util::Vector3D<site_t>(CubeSize + 1);
            CPPUNIT_ASSERT(!geometrySelector->Include(*dataSourceIterator, invalidLocation));
          }

          void TestExpectedIncludedSites(hemelb::extraction::GeometrySelector* geometrySelector,
                                         std::vector<util::Vector3D<site_t> >& includedSites)
          {
            dataSourceIterator->Reset();
            while (dataSourceIterator->ReadNext())
            {
              bool expectedIncluded = std::count(includedSites.begin(),
                                                 includedSites.end(),
                                                 dataSourceIterator->GetPosition()) > 0;

              std::stringstream msg;

              msg << "Site at " << dataSourceIterator->GetPosition().x << ","
                  << dataSourceIterator->GetPosition().y << ","
                  << dataSourceIterator->GetPosition().z << " was ";

              if (!expectedIncluded)
              {
                msg << "not ";
              }

              msg << "expected to be included but actually was";

              if (!geometrySelector->Include(*dataSourceIterator, dataSourceIterator->GetPosition()))
              {
                msg << " not.";
              }
              else
              {
                msg << ".";
              }

              CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(),
                                           expectedIncluded,
                                           geometrySelector->Include(*dataSourceIterator,
                                                                     dataSourceIterator->GetPosition()));
            }
          }

          const double VoxelSize;
          const int CubeSize;
          const distribn_t CentreCoordinate;

          const util::Vector3D<distribn_t> planeNormal;
          const util::Vector3D<distribn_t> planePosition;
          const distribn_t planeRadius;
          const util::Vector3D<distribn_t> lineEndPoint1;
          const util::Vector3D<distribn_t> lineEndPoint2;
          const util::Vector3D<distribn_t> surfacePoint;
          const util::Vector3D<distribn_t> surfacePointMultipleHits;

          unittests::FourCubeLatticeData* latticeData;
          lb::SimulationState* simState;
          hemelb::lb::MacroscopicPropertyCache* propertyCache;
          hemelb::util::UnitConverter* unitConverter;
          hemelb::extraction::LbDataSourceIterator* dataSourceIterator;

          hemelb::extraction::PlaneGeometrySelector* planeGeometrySelector;
          hemelb::extraction::PlaneGeometrySelector* planeGeometrySelectorWithRadius;
          hemelb::extraction::StraightLineGeometrySelector* straightLineGeometrySelector;
          hemelb::extraction::WholeGeometrySelector* wholeGeometrySelector;
          hemelb::extraction::GeometrySurfaceSelector* geometrySurfaceSelector;
          hemelb::extraction::SurfacePointSelector* surfacePointSelector;
          hemelb::extraction::SurfacePointSelector* surfacePointSelectorMultipleHits;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION ( GeometrySelectorTests);

    }
  }
}

#endif /* HEMELB_UNITTESTS_EXTRACTION_GEOMETRYSELECTORTESTS_H */
