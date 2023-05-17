// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cassert>
#include "lb/kernels/momentBasis/DHumieresD3Q15MRTBasis.h"

namespace hemelb::lb
{
    auto DHumieresD3Q15MRTBasis::SetUpCollisionMatrix(distribn_t tau) -> std::array<distribn_t, NUMMOMENTS>
    {
        // Relaxation values taken from d'Humieres 2002, except for the kinematic viscosity where the usual tau formula is used.
        return {
                1.6, // e (s1)
                1.2, // epsilon (s2)
                1.6, // q_x (s4)
                1.6, // q_y (s4)
                1.6, // q_z (s4)
                1.0 / tau, // 3p_xx (s9)
                1.0 / tau, // p_ww (s9)
                1.0 / tau, // p_xy (s11 = s9)
                1.0 / tau, // p_yz (s11 = s9)
                1.0 / tau, // p_zx (s11 = s9)
                1.2 // m_xyz (s14)
        };
    }
}
