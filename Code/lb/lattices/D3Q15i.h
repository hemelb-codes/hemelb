// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_LATTICES_D3Q15I_H
#define HEMELB_LB_LATTICES_D3Q15I_H

#include "lb/lattices/Lattice.h"
#include "lb/lattices/D3Q15.h"

namespace hemelb::lb {
    // Use inheritance rather than an alias to get a nicer name when compiling/debugging.
    struct D3Q15i : Lattice<
            15,
            D3Q15::VECTORS,
            D3Q15::EQMWEIGHTS,
            false
    > {};
}
#endif
