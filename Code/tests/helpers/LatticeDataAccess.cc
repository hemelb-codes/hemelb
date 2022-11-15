// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/LatticeDataAccess.h"

namespace hemelb::tests::helpers
{

    LatticeDataAccess::LatticeDataAccess(geometry::FieldData * const latDat) :
            latDat(latDat)
    {
    }

    distribn_t const * LatticeDataAccess::GetFNew(site_t index) const
    {
        return latDat->GetFNew(index);
    }

    void LatticeDataAccess::SetMinWallDistance(PhysicalDistance _mindist)
    {
        for (auto& dist: domain->distanceToWall) {
            if (dist > 0e0)
                dist = _mindist;
        }
    }

    void LatticeDataAccess::SetWallDistance(PhysicalDistance _mindist)
    {
        for (auto& dist: domain->distanceToWall) {
            if (dist > 0e0)
                dist = _mindist;
        }
    }

    distribn_t const * GetFNew(geometry::FieldData& latDat, site_t const &index)
    {
        return LatticeDataAccess(&latDat).GetFNew(index);
    }

}
