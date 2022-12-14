// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LATTICES_D3Q27_H
#define HEMELB_LB_LATTICES_D3Q27_H

#include "lb/lattices/Lattice.h"

namespace hemelb::lb::lattices
{
    using D3Q27 = Lattice<
            27,
            std::array<util::Vector3D<int>, 27>{
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
                            {1, 1, 1},
                            {-1, -1, -1},
                            {1, 1, -1},
                            {-1, -1, 1},
                            {1, -1, 1},
                            {-1, 1, -1},
                            {1, -1, -1},
                            {-1, 1, 1}
                    }
            },
            std::array<distribn_t , 27>{
                    8.0 / 27.0,
                    2.0 / 27.0,
                    2.0 / 27.0,
                    2.0 / 27.0,
                    2.0 / 27.0,
                    2.0 / 27.0,
                    2.0 / 27.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 54.0,
                    1.0 / 216.0,
                    1.0 / 216.0,
                    1.0 / 216.0,
                    1.0 / 216.0,
                    1.0 / 216.0,
                    1.0 / 216.0,
                    1.0 / 216.0,
                    1.0 / 216.0
            }
    >;

}

#endif
