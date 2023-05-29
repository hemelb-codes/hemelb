// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TRAITS_H
#define HEMELB_TRAITS_H

#include "lb/Lattices.h"
#include "lb/Kernels.h"
#include "lb/Streamers.h"
#include "lb/Collisions.h"
#include "redblood/stencil.h"

namespace hemelb
{
    template <
        typename LATTICE = lb::DefaultLattice,
        template<lb::lattice_type> class KERNEL = lb::DefaultKernel,
        template<class> class COLLISION = lb::Normal,
        template<class> class STREAMER = lb::DefaultStreamer,
        template<class> class WALL_BOUNDARY = lb::DefaultWallStreamer,
        template<class> class INLET_BOUNDARY = lb::DefaultInletStreamer,
        template<class> class OUTLET_BOUNDARY = lb::DefaultOutletStreamer,
        typename STENCIL = redblood::stencil::DefaultStencil
    >
    struct Traits
    {
        using Lattice = LATTICE;
        using Kernel = KERNEL<Lattice>;
        using Collision = COLLISION<Kernel>;
        using Streamer = STREAMER<Collision>;
        using WallBoundary = WALL_BOUNDARY<Collision>;
        using InletBoundary = INLET_BOUNDARY<Collision>;
        using OutletBoundary = OUTLET_BOUNDARY<Collision>;
        using WallInletBoundary = typename lb::CombineWallAndIoletStreamers<WallBoundary, InletBoundary>::type;
        using WallOutletBoundary = typename lb::CombineWallAndIoletStreamers<WallBoundary, OutletBoundary>::type;
        using Stencil = STENCIL;
    };
}

#endif
