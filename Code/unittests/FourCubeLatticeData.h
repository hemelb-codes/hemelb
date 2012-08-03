// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H
#define HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H

#include <cstdlib>
#include "units.h"
#include "geometry/LatticeData.h"
#include "io/formats/geometry.h"

namespace hemelb
{
  namespace unittests
  {
    class TestSiteData : public geometry::SiteData
    {
      public:
        TestSiteData(geometry::SiteData& siteData) :
            geometry::SiteData(siteData)
        {

        }

        void SetHasBoundary(Direction direction)
        {
          unsigned newValue = geometry::SiteData::GetIntersectionData();
          newValue |= 1U << (direction - 1);
          boundaryIntersection = newValue;
        }
    };

    class FourCubeLatticeData : public geometry::LatticeData
    {
      public:

        /**
         * The create function makes a 4 x 4 x 4 cube of sites from (0,0,0) to (3,3,3).
         * The plane (x,y,0) is an inlet (boundary 0).
         * The plane (x,y,3) is an outlet (boundary 1).
         * The planes (0,y,z), (3,y,z), (x,0,z) and (x,3,z) are all walls.
         *
         * @return
         */
        static FourCubeLatticeData* Create(site_t sitesPerBlockUnit = 4, proc_t rankCount = 1)
        {
          hemelb::geometry::Geometry readResult(util::Vector3D<site_t>::Ones(),
                                                sitesPerBlockUnit,
                                                0.01,
                                                util::Vector3D<PhysicalLength>::Zero());

          hemelb::geometry::BlockReadResult& block = readResult.Blocks[0];
          block.Sites.resize(readResult.GetSitesPerBlock(), geometry::GeometrySite(false));

          site_t index = -1;
          for (site_t i = 0; i < sitesPerBlockUnit; ++i)
          {
            for (site_t j = 0; j < sitesPerBlockUnit; ++j)
            {
              for (site_t k = 0; k < sitesPerBlockUnit; ++k)
              {
                ++index;

                hemelb::geometry::GeometrySite& site = block.Sites[index];

                site.isFluid = true;
                site.targetProcessor = 0;

                for (Direction direction = 1; direction < lb::lattices::D3Q15::NUMVECTORS; ++direction)
                {
                  site_t neighI = i + lb::lattices::D3Q15::CX[direction];
                  site_t neighJ = j + lb::lattices::D3Q15::CY[direction];
                  site_t neighK = k + lb::lattices::D3Q15::CZ[direction];

                  hemelb::geometry::GeometrySiteLink link;

                  float randomDistance = (float(std::rand() % 10000) / 10000.0);

                  // The inlet is by the minimal z value.
                  if (neighK < 0)
                  {
                    link.ioletId = 0;
                    link.type = geometry::GeometrySiteLink::INLET_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }
                  // The outlet is by the maximal z value.
                  else if (neighK >= sitesPerBlockUnit)
                  {
                    link.ioletId = 0;
                    link.type = geometry::GeometrySiteLink::OUTLET_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }
                  // Walls are by extremes of x and y.
                  else if (neighI < 0 || neighJ < 0 || neighI >= sitesPerBlockUnit || neighJ >= sitesPerBlockUnit)
                  {
                    link.type = geometry::GeometrySiteLink::WALL_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }

                  site.links.push_back(link);
                }
              }
            }
          }

          FourCubeLatticeData* returnable = new FourCubeLatticeData(readResult);

          // First, fiddle with the fluid site count, for tests that require this set.
          returnable->fluidSitesOnEachProcessor.resize(rankCount);
          returnable->fluidSitesOnEachProcessor[0] = sitesPerBlockUnit * sitesPerBlockUnit * sitesPerBlockUnit;
          for (proc_t rank = 1; rank < rankCount; ++rank)
          {
            returnable->fluidSitesOnEachProcessor[rank] = rank * 1000;
          }

          return returnable;
        }

        /***
         Not used in setting up the four cube, but used in other tests to poke changes into the four cube for those tests.
         **/
        void SetHasBoundary(site_t site, Direction direction)
        {
          TestSiteData mutableSiteData(siteData[site]);
          mutableSiteData.SetHasBoundary(direction);
          siteData[site] = geometry::SiteData(mutableSiteData);
        }

        /***
         Not used in setting up the four cube, but used in other tests to poke changes into the four cube for those tests.
         **/
        void SetBoundaryDistance(site_t site, Direction direction, distribn_t distance)
        {
          distanceToWall[ (lb::lattices::D3Q15::NUMVECTORS - 1) * site + direction - 1] = distance;
        }

      protected:
        FourCubeLatticeData(hemelb::geometry::Geometry& readResult) :
            hemelb::geometry::LatticeData(lb::lattices::D3Q15::GetLatticeInfo(), readResult)
        {

        }

        FourCubeLatticeData() :
            hemelb::geometry::LatticeData(lb::lattices::D3Q15::GetLatticeInfo())
        {

        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H */
