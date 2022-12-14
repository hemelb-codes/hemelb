// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LATTICES_D3Q19_H
#define HEMELB_LB_LATTICES_D3Q19_H

#include "lb/lattices/Lattice.h"

namespace hemelb::lb::lattices
{
    using D3Q19 = Lattice<
            19,
            std::array<util::Vector3D<int>, 19>{
                    {
                            {0, 0, 0},
                            {1, 0, 0},
                            {-1, 0, 0},
                            {0, 1, 0},
                            {0, -1, 0},
                            {0, 0, 1},
                            {0, 0, -1},
                            {1, 1, 0},
                            {-1, -1, 0},
                            {1, -1, 0},
                            {-1, 1, 0},
                            {1, 0, 1},
                            {-1, 0, -1},
                            {1, 0, -1},
                            {-1, 0, 1},
                            {0, 1, 1},
                            {0, -1, -1},
                            {0, 1, -1},
                            {0, -1, 1},
                    }
            },
            std::array<distribn_t, 19>{
                    1.0 / 3.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0
                                                                   / 18.0,
                    1.0 / 18.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0
                                                                    / 36.0,
                    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0
                                                                    / 36.0,
                    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
            }
    >;

}

#endif
