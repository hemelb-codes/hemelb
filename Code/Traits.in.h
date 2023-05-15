// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TRAITS_IN_H
#define HEMELB_TRAITS_IN_H

#include "lb/streamers/SimpleCollideAndStream.h"
#include "lb/collisions/Normal.h"
#include "lb/BuildSystemInterface.h"
#include "redblood/stencil.h"

namespace hemelb
{
    template <
        typename LATTICE = lb::@HEMELB_LATTICE@,
        template<class> class KERNEL = lb::@HEMELB_KERNEL@,
        template<class> class COLLISION = lb::collisions::Normal,
        template<class> class STREAMER = lb::streamers::SimpleCollideAndStream,
        template<class> class WALL_BOUNDARY = lb::@HEMELB_WALL_BOUNDARY@,
        template<class> class INLET_BOUNDARY = lb::@HEMELB_INLET_BOUNDARY@,
        template<class> class OUTLET_BOUNDARY = lb::@HEMELB_OUTLET_BOUNDARY@,
        template<class> class WALL_INLET_BOUNDARY = lb::@HEMELB_WALL_INLET_BOUNDARY@,
        template<class> class WALL_OUTLET_BOUNDARY = lb::@HEMELB_WALL_OUTLET_BOUNDARY@,
        typename STENCIL = redblood::stencil::@HEMELB_STENCIL@
    >
    struct Traits
    {
        using Lattice = LATTICE;
        using Kernel = typename KERNEL<Lattice>::Type;
        using Collision = COLLISION<Kernel>;
        using Streamer = STREAMER<Collision>;
        using WallBoundary = typename WALL_BOUNDARY<Collision>::Type;
        using InletBoundary = typename INLET_BOUNDARY<Collision>::Type;
        using OutletBoundary = typename OUTLET_BOUNDARY<Collision>::Type;
        using WallInletBoundary = typename WALL_INLET_BOUNDARY<Collision>::Type;
        using WallOutletBoundary = typename WALL_OUTLET_BOUNDARY<Collision>::Type;
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
            template<class> class NEW_WALL_INLET_BOUNDARY = WALL_INLET_BOUNDARY,
            template<class> class NEW_WALL_OUTLET_BOUNDARY = WALL_OUTLET_BOUNDARY,
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
                    NEW_WALL_INLET_BOUNDARY,
                    NEW_WALL_OUTLET_BOUNDARY,
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
                    WALL_INLET_BOUNDARY, WALL_OUTLET_BOUNDARY,
                    NEW_STENCIL
            >::Type;
        };
    };
}

#endif
