// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_SITETRAVERSER_H
#define HEMELB_GEOMETRY_SITETRAVERSER_H

#include "geometry/VolumeTraverser.h"

namespace hemelb::geometry
{
    class Domain;

    /**
     * SiteTraverser is used to traverse the sites in a speficied block
     * within the lattice data
     */
    class SiteTraverser : public VolumeTraverser
    {
    public:
        SiteTraverser(const Domain & iLatticeDat);

        ~SiteTraverser() override = default;

        U16 GetXCount() const override;

        U16 GetYCount() const override;

        U16 GetZCount() const override;

    private:
        //Returns the block size in number of sites
        U16 GetBlockSize() const;

        const Domain & mLatticeData;
    };

}

#endif // HEMELB_GEOMETRY_SITETRAVERSER_H
