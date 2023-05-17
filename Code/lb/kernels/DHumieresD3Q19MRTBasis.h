// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_DHUMIERESD3Q19MRTBASIS_H
#define HEMELB_LB_KERNELS_DHUMIERESD3Q19MRTBASIS_H

#include "lb/lattices/D3Q19.h"
#include "lb/kernels/basis_helpers.h"

namespace hemelb::lb
{
    /**
     *  Class implementing the Multiple Relaxation Time (MRT) moment basis presented in in d'Humieres et al. (2002)
     *  "Multiple-relaxation-time lattice Boltzmann models in three dimensions" for the D3Q19 lattice
     */
    class DHumieresD3Q19MRTBasis
    {
    public:
        /**
         * The lattice that suits this basis.
         */
        using Lattice = D3Q19;

        /**
         * Moments can be separated into two groups:
         * a) conserved (density and momentum) and
         * b) non-conserved (stress, kinetic/ghost modes).
         */
        static constexpr unsigned NUMMOMENTS = 15;

        using Traits = moment_traits<distribn_t, NUMMOMENTS, Lattice::NUMVECTORS>;
        using MatrixType = typename Traits::MatrixType;

        /**
         * Matrix used to convert from the velocities space to the reduced moment space containing only kinetic moments.
         *
         * Kinetic moments defined in d'Humieres et al. 2002. To get the matrix below, swap
         * columns 9 and 11, 13 and 15, 17 and 19 to match HemeLB's lattice velocity ordering.
         *
         * See publication for  meaning of e, epsilon, etc.
         */
        static constexpr MatrixType REDUCED_MOMENT_BASIS =
                {{ { -30, -11, -11, -11, -11, -11, -11, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 }, //e
                   { 12, -4, -4, -4, -4, -4, -4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, // epsilon
                   { 0, -4, 4, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0 }, // q_x
                   { 0, 0, 0, -4, 4, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1 }, //q_y
                   { 0, 0, 0, 0, 0, -4, 4, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1 }, // q_z
                   { 0, 2, 2, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2 }, // 3p_xx
                   { 0, -4, -4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2 }, // 3pi_xx
                   { 0, 0, 0, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0 }, // p_ww
                   { 0, 0, 0, -2, -2, 2, 2, 1, 1, 1, 1, -1, -1, -1, -1, 0, 0, 0, 0 }, // pi_ww
                   { 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0 }, // p_xy
                   { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1 }, // p_yz
                   { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0 }, // p_xz
                   { 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0 }, // m_x
                   { 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1 }, // m_y
                   { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1 } //m_z
                 }};

        /** Diagonal matrix REDUCED_MOMENT_BASIS * REDUCED_MOMENT_BASIS'. */
        static constexpr auto BASIS_TIMES_BASIS_TRANSPOSED = Traits::DiagSelfProduct(REDUCED_MOMENT_BASIS);

        /**
         * Sets up the MRT collision matrix \hat{S}
         *
         * @return collisionMatrix MRT collision matrix, diagonal
         * @param tau LB relaxation time used to relax some of the moments
         */
        static std::array<distribn_t, NUMMOMENTS>
        SetUpCollisionMatrix(distribn_t tau);
    };

}
#endif //HEMELB_LB_KERNELS_DHUMIERESD3Q19MRTBASIS_H
