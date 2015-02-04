//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_TRAITS_H
#define HEMELB_TRAITS_H

#include "lb/streamers/SimpleCollideAndStream.h"
#include "lb/collisions/Normal.h"
#include "lb/BuildSystemInterface.h"

namespace hemelb
{
  template
  <
    typename LATTICE = lb::lattices::HEMELB_LATTICE,
    template<class> class KERNEL = lb::HEMELB_KERNEL,
    template<class> class STREAMER = lb::streamers::SimpleCollideAndStream,
    template<class> class COLLISION = lb::collisions::Normal,
    template<class> class WALL_BOUNDARY = lb::HEMELB_WALL_BOUNDARY,
    template<class> class INLET_BOUNDARY = lb::HEMELB_INLET_BOUNDARY,
    template<class> class OUTLET_BOUNDARY = lb::HEMELB_OUTLET_BOUNDARY,
    template<class> class WALL_INLET_BOUNDARY = lb::HEMELB_WALL_INLET_BOUNDARY,
    template<class> class WALL_OUTLET_BOUNDARY = lb::HEMELB_WALL_OUTLET_BOUNDARY
  > struct Traits
  {
    typedef LATTICE Lattice;
    typedef typename KERNEL<Lattice>::Type Kernel;
    typedef COLLISION<Kernel> Collision;
    typedef STREAMER<Collision> Streamer;
    typedef typename WALL_BOUNDARY<Collision>::Type WallBoundary;
    typedef typename INLET_BOUNDARY<Collision>::Type InletBoundary;
    typedef typename OUTLET_BOUNDARY<Collision>::Type OutletBoundary;
    typedef typename WALL_INLET_BOUNDARY<Collision>::Type WallInletBoundary;
    typedef typename WALL_OUTLET_BOUNDARY<Collision>::Type WallOutletBoundary;
  };
}
#endif /* HEMELB_TRAITS */
