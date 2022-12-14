// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LATTICES_D3Q15_H
#define HEMELB_LB_LATTICES_D3Q15_H

#include "lb/lattices/Lattice.h"

namespace hemelb::lb::lattices
{

    using D3Q15 = Lattice<15,
            std::array<util::Vector3D<int>, 15>{{
                                                        { 0, 0, 0},

                                                        { 1, 0, 0},
                                                        {-1, 0, 0},
                                                        { 0, 1, 0},
                                                        { 0,-1, 0},
                                                        { 0, 0, 1},
                                                        { 0, 0,-1},

                                                        { 1, 1, 1},
                                                        {-1,-1,-1},
                                                        { 1, 1,-1},
                                                        {-1,-1, 1},
                                                        { 1,-1, 1},
                                                        {-1, 1,-1},
                                                        { 1,-1,-1},
                                                        {-1, 1, 1}
                                                }},
            std::array<distribn_t, 15>{
                    2.0 / 9.0,

                    1.0 / 9.0,
                    1.0 / 9.0,
                    1.0 / 9.0,
                    1.0 / 9.0,
                    1.0 / 9.0,
                    1.0 / 9.0,

                    1.0 / 72.0,
                    1.0 / 72.0,
                    1.0 / 72.0,
                    1.0 / 72.0,
                    1.0 / 72.0,
                    1.0 / 72.0,
                    1.0 / 72.0,
                    1.0 / 72.0
            }>;
}
#endif
