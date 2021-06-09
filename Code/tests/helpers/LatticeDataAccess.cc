// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/LatticeDataAccess.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {

      LatticeDataAccess::LatticeDataAccess(geometry::LatticeData * const latDat) :
	latDat(latDat)
      {
      }

      distribn_t const * LatticeDataAccess::GetFNew(site_t index) const
      {
	return latDat->GetFNew(index);
      }

      void LatticeDataAccess::SetMinWallDistance(PhysicalDistance _mindist)
      {
        typedef std::vector<distribn_t>::iterator iterator;
        iterator i_first = latDat->distanceToWall.begin();
        iterator const i_end = latDat->distanceToWall.end();
        for (; i_first != i_end; ++i_first)
          if (*i_first > 0e0 and *i_first < _mindist)
            *i_first = _mindist;
      }

      void LatticeDataAccess::SetWallDistance(PhysicalDistance _mindist)
      {
        typedef std::vector<distribn_t>::iterator iterator;
        iterator i_first = latDat->distanceToWall.begin();
        iterator const i_end = latDat->distanceToWall.end();
        for (; i_first != i_end; ++i_first)
          if (*i_first > 0e0)
            *i_first = _mindist;
      }

      distribn_t const * GetFNew(geometry::LatticeData *latDat, site_t const &index)
      {
        return LatticeDataAccess(latDat).GetFNew(index);
      }

      void SetMinWallDistance(geometry::LatticeData * const latDat, PhysicalDistance _mindist)
      {
        LatticeDataAccess(latDat).SetMinWallDistance(_mindist);
      }
      void SetWallDistance(geometry::LatticeData * const latDat, PhysicalDistance _mindist)
      {
        LatticeDataAccess(latDat).SetWallDistance(_mindist);
      }

    }
  }
}
