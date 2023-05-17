// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/kernels/momentBasis/DHumieresD3Q19MRTBasis.h"

namespace hemelb::lb
{

    auto DHumieresD3Q19MRTBasis::SetUpCollisionMatrix(distribn_t tau) -> std::array <distribn_t, NUMMOMENTS>
    {
        // Relaxation values taken from d'Humieres 2002, except for the kinematic viscosity where the usual tau formula is used.
        return {
                1.19, // e (s1)

                1.4, // epsilon (s2)

                1.2, // q_x (s4)
                1.2, // q_y (s4)
                1.2, // q_z (s4)

                1.0 / tau, // 3p_xx (s9)
                1.4, // 3pi_xx s10
                1.0 / tau, // 3p_ww (s9)
                1.4, // 3pi_ww s10

                1.0 / tau, // p_xy (s13 = s9)
                1.0 / tau, // p_yz (s13 = s9)
                1.0 / tau, // p_xz (s13 = s9)

                1.98, // m_x (s16)
                1.98, // m_y (s16)
                1.98 // m_z (s16)
        };
    }
}
