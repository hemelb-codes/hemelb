// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_FOURCUBELATTICEDATA_H
#define HEMELB_TESTS_HELPERS_FOURCUBELATTICEDATA_H

#include <cstdlib>
#include "units.h"
#include "geometry/LatticeData.h"
#include "io/formats/geometry.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace tests
  {
    class TestSiteData : public geometry::SiteData
    {
      public:
      TestSiteData(geometry::SiteData& siteData);

      void SetHasWall(Direction direction);
      void SetHasIolet(Direction direction);
      void SetIoletId(int boundaryId);
    };

    class FourCubeLatticeData : public geometry::LatticeData
    {
      public:

      // The create function makes a 4 x 4 x 4 cube of sites from (0,0,0) to (3,3,3).
      // The plane (x,y,0) is an inlet (boundary 0).
      // The plane (x,y,3) is an outlet (boundary 1).
      // The planes (0,y,z), (3,y,z), (x,0,z) and (x,3,z) are all walls.
      static FourCubeLatticeData* Create(const net::IOCommunicator& comm, site_t sitesPerBlockUnit = 6, proc_t rankCount = 1);

      // Not used in setting up the four cube, but used in other tests
      // to poke changes into the four cube for those tests.
      void SetHasWall(site_t site, Direction direction);
      void SetHasIolet(site_t site, Direction direction);
      void SetIoletId(site_t site, int id);
      void SetBoundaryDistance(site_t site, Direction direction, distribn_t distance);
      void SetBoundaryNormal(site_t site, util::Vector3D<distribn_t> boundaryNormal);

      // Used in unit tests for setting the fOld array, in a way that isn't possible in the main
      // part of the codebase.
      // @param site
      // @param fOldIn
      template<class LatticeType>
      void SetFOld(site_t site, distribn_t* fOldIn)
      {
	for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction) {
            *GetFOld(site * LatticeType::NUMVECTORS + direction) = fOldIn[direction];
          }
        }

      protected:
      FourCubeLatticeData(hemelb::geometry::Geometry& readResult, const net::IOCommunicator& comms);
    };
  }
}

#endif // ONCE
