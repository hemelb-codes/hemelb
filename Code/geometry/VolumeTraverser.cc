// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "constants.h"
#include "geometry/VolumeTraverser.h"

namespace hemelb::geometry
{

    Vec16 const& VolumeTraverser::GetCurrentLocation()
    {
        return mCurrentLocation;
    }

    site_t VolumeTraverser::GetCurrentIndex() const
    {
        return mCurrentNumber;
    }

    site_t VolumeTraverser::GetIndexFromLocation(Vec16 const& iLocation) const
    {
        return ( (iLocation.x() * GetYCount() + iLocation.y()) * GetZCount()) + iLocation.z();
    }

    bool VolumeTraverser::TraverseOne()
    {
        mCurrentNumber++;

        mCurrentLocation.z()++;
        if (mCurrentLocation.z() < GetZCount())
        {
            return true;
        }

        mCurrentLocation.z() = 0;
        mCurrentLocation.y()++;
        if (mCurrentLocation.y() < GetYCount())
        {
            return true;
        }

        mCurrentLocation.y() = 0;
        mCurrentLocation.x()++;

        if (mCurrentLocation.x() < GetXCount())
        {
            return true;
        }
        return false;
    }

    bool VolumeTraverser::CurrentLocationValid()
    {
        if (GetCurrentIndex() < 0)
        {
            return false;
        }

        if (!mCurrentLocation.IsInRange(
                Vec16::Zero(),
                Vec16(GetXCount() - 1, GetYCount() - 1, GetZCount() - 1)
        ))
        {
            return false;
        }

        return true;
    }
}
