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

        //! Fully reinstantiate, where defaults are current choices.
        //! This is a convenience function to limit the amount of explicit resintantiation
        template <
            class NEW_LATTICE,
            template<class> class NEW_KERNEL = KERNEL,
            template<class> class NEW_COLLISION = COLLISION,
            template<class> class NEW_STREAMER = STREAMER,
            template<class> class NEW_WALL_BOUNDARY = WALL_BOUNDARY,
            template<class> class NEW_INLET_BOUNDARY = INLET_BOUNDARY,
            template<class> class NEW_OUTLET_BOUNDARY = OUTLET_BOUNDARY,
            typename NEW_STENCIL = STENCIL
        >
        struct Reinstantiate
        {
            using Type = Traits<
                    NEW_LATTICE,
                    NEW_KERNEL,
                    NEW_COLLISION,
                    NEW_STREAMER,
                    NEW_WALL_BOUNDARY,
                    NEW_INLET_BOUNDARY,
                    NEW_OUTLET_BOUNDARY,
                    NEW_STENCIL
            >;
        };

        //! Changes only kernel type
        //! This is a convenience function to limit the amount of explicit resintantiation
        template<template<class> class NEW_KERNEL>
        struct ChangeKernel
        {
            using Type = typename Reinstantiate<LATTICE, NEW_KERNEL>::Type;
        };

        //! Changes only stencil type
        //! This is a convenience function to limit the amount of explicit resintantiation
        template<class NEW_STENCIL>
        struct ChangeStencil
        {
            using Type = typename Reinstantiate<
                    LATTICE, KERNEL, COLLISION, STREAMER,
                    WALL_BOUNDARY, INLET_BOUNDARY, OUTLET_BOUNDARY,
                    NEW_STENCIL
            >::Type;
        };
    };
}

#endif
