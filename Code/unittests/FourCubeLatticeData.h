
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H
#define HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H

#include <cstdlib>
#include "units.h"
#include "geometry/LatticeData.h"
#include "io/formats/geometry.h"
#include "util/Vector3D.h"

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

        void SetHasWall(Direction direction)
        {
          unsigned newValue = geometry::SiteData::GetWallIntersectionData();
          newValue |= 1U << (direction - 1);
          wallIntersection = newValue;
        }

        void SetHasIolet(Direction direction)
        {
          unsigned newValue = geometry::SiteData::GetIoletIntersectionData();
          newValue |= 1U << (direction - 1);
          ioletIntersection = newValue;
        }

        void SetIoletId(int boundaryId)
        {
          ioletId = boundaryId;
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
        static FourCubeLatticeData* Create(const net::IOCommunicator& comm, site_t sitesPerBlockUnit = 6, proc_t rankCount = 1)
        {
          hemelb::geometry::Geometry readResult(util::Vector3D<site_t>::Ones(),
                                                sitesPerBlockUnit);
// VoxelSize                                               0.01,
// Origin                                               util::Vector3D<PhysicalDistance>::Zero());
          site_t sitesAlongCube = sitesPerBlockUnit - 2;
          site_t minInd = 1, maxInd = sitesAlongCube;

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

                if (i < minInd || i > maxInd || j < minInd || j > maxInd || k < minInd || k
                    > maxInd)
                  continue;

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
                  if (neighK < minInd)
                  {
                    link.ioletId = 0;
                    link.type = geometry::GeometrySiteLink::INLET_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }
                  // The outlet is by the maximal z value.
                  else if (neighK > maxInd)
                  {
                    link.ioletId = 0;
                    link.type = geometry::GeometrySiteLink::OUTLET_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }
                  // Walls are by extremes of x and y.
                  else if (neighI < minInd || neighJ < minInd || neighI > maxInd || neighJ > maxInd)
                  {
                    link.type = geometry::GeometrySiteLink::WALL_INTERSECTION;
                    link.distanceToIntersection = randomDistance;
                  }

                  site.links.push_back(link);

                }

                /*
                 * For sites at the intersection of two cube faces considered wall (i.e. perpendicular to the x or y
                 * axes), we arbitrarily choose the normal to lie along the y axis. The logic below must be consistent
                 * with Scripts/SimpleGeometryGenerationScripts/four_cube.py
                 */
                if (i == minInd)
                {
                  site.wallNormalAvailable = true;
                  site.wallNormal = util::Vector3D<float>(-1, 0, 0);
                }
                if (i == maxInd)
                {
                  site.wallNormalAvailable = true;
                  site.wallNormal = util::Vector3D<float>(1, 0, 0);
                }
                if (j == minInd)
                {
                  site.wallNormalAvailable = true;
                  site.wallNormal = util::Vector3D<float>(0, -1, 0);
                }
                if (j == maxInd)
                {
                  site.wallNormalAvailable = true;
                  site.wallNormal = util::Vector3D<float>(0, 1, 0);
                }
              }
            }
          }

          FourCubeLatticeData* returnable = new FourCubeLatticeData(readResult, comm);

          // First, fiddle with the fluid site count, for tests that require this set.
          returnable->fluidSitesOnEachProcessor.resize(rankCount);
          returnable->fluidSitesOnEachProcessor[0] = sitesAlongCube * sitesAlongCube
              * sitesAlongCube;
          for (proc_t rank = 1; rank < rankCount; ++rank)
          {
            returnable->fluidSitesOnEachProcessor[rank] = rank * 1000;
          }

          return returnable;
        }

        /***
         Not used in setting up the four cube, but used in other tests to poke changes into the four cube for those tests.
         **/
        void SetHasWall(site_t site, Direction direction)
        {
          TestSiteData mutableSiteData(siteData[site]);
          mutableSiteData.SetHasWall(direction);
          siteData[site] = geometry::SiteData(mutableSiteData);
        }

        /***
         Not used in setting up the four cube, but used in other tests to poke changes into the four cube for those tests.
         **/
        void SetHasIolet(site_t site, Direction direction)
        {
          TestSiteData mutableSiteData(siteData[site]);
          mutableSiteData.SetHasIolet(direction);
          siteData[site] = geometry::SiteData(mutableSiteData);
        }

        /***
         Not used in setting up the four cube, but used in other tests to poke changes into the four cube for those tests.
         **/
        void SetIoletId(site_t site, int id)
        {
          TestSiteData mutableSiteData(siteData[site]);
          mutableSiteData.SetIoletId(id);
          siteData[site] = geometry::SiteData(mutableSiteData);
        }

        /***
         Not used in setting up the four cube, but used in other tests to poke changes into the four cube for those tests.
         **/
        void SetBoundaryDistance(site_t site, Direction direction, distribn_t distance)
        {
          distanceToWall[ (lb::lattices::D3Q15::NUMVECTORS - 1) * site + direction - 1] = distance;
        }

        /***
         Not used in setting up the four cube, but used in other tests to poke changes into the four cube for those tests.
         **/
        void SetBoundaryNormal(site_t site, util::Vector3D<distribn_t> boundaryNormal)
        {
          wallNormalAtSite[site] = boundaryNormal;
        }

        /**
         * Used in unit tests for setting the fOld array, in a way that isn't possible in the main
         * part of the codebase.
         * @param site
         * @param fOldIn
         */
        template<class LatticeType>
        void SetFOld(site_t site, distribn_t* fOldIn)
        {
          for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
          {
            *GetFOld(site * LatticeType::NUMVECTORS + direction) = fOldIn[direction];
          }
        }

      protected:
        FourCubeLatticeData(hemelb::geometry::Geometry& readResult, const net::IOCommunicator& comms) :
          hemelb::geometry::LatticeData(lb::lattices::D3Q15::GetLatticeInfo(), readResult, comms)
        {

        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H */
