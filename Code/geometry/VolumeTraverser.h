// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_VOLUMETRAVERSER_H
#define HEMELB_GEOMETRY_VOLUMETRAVERSER_H

#include <cstdint>

#include "util/Vector3D.h"

namespace hemelb::geometry
{
    /**
     * VolumeTraverser is used to sequentially traverse a 
     * 3D structure maintaining the index and Location 
     * within volume
     */
    class VolumeTraverser
    {
    public:
        virtual ~VolumeTraverser() = default;

        Vec16 const& GetCurrentLocation();

        site_t GetCurrentIndex() const;

        site_t GetIndexFromLocation(Vec16 const& iLocation) const;

        //Increments the index by one and update the location accordingly
        //Returns true if successful or false if the whole volume has been
        //traversed
        bool TraverseOne();

        bool CurrentLocationValid();

        //Virtual methods which must be defined for correct traversal
        virtual U16 GetXCount() const = 0;
        virtual U16 GetYCount() const = 0;
        virtual U16 GetZCount() const = 0;

    protected:
        VolumeTraverser() = default;

    private:
        Vec16 mCurrentLocation = Vec16::Zero();
        site_t mCurrentNumber = 0;
    };

}

#endif // HEMELB_GEOMETRY_VOLUMETRAVERSER_H
