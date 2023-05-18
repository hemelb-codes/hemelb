// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_MRT_H
#define HEMELB_LB_KERNELS_MRT_H

#include "lb/SimulationState.h"
#include <cassert>
#include <cmath>

namespace hemelb::lb
{
    // Forward declaration needed by the struct
    template<moment_basis> class MRT;

    template<class MomentBasis>
    struct HydroVars<MRT<MomentBasis> > : public HydroVarsBase<typename MomentBasis::Lattice>
    {
        using HydroVarsBase<typename MomentBasis::Lattice>::HydroVarsBase;
        /** Equilibrium velocity distribution in the momentum space. */
        std::array<distribn_t, MomentBasis::NUMMOMENTS> m_neq;
    };

    /**
     * This class implements the MRT collision operator.
     *
     *  \Omega(f) = - M^{-1} * \hat{S} * M (f - f_{eq})
     *            = - M^T * (M * M^T)^{-1} * \hat{S} * m_{neq}
     *
     *  where f is the velocity distribution function, \Omega() is the collision operator,
     *  M is the momentum space basis, m=Mf is the vector of momentums, \hat{S} is the
     *  collision matrix, and {m,f}_{eq} and {m,f}_{neq} are the equilibrium and non-equilibrium
     *  versions of {m,f}.
     *
     *  (M * M^T)^{-1} and \hat{S} are diagonal matrices.
     */
    template<moment_basis M>
    class MRT
    {
    public:
        using MomentType = M;
        using LatticeType = typename MomentType::Lattice;
        using VarsType = HydroVars<MRT>;
        using MatrixType = typename MomentType::MatrixType;

        static constexpr std::size_t NUMMOMENTS = MomentType::NUMMOMENTS;
        static constexpr std::size_t NUMVECTORS = LatticeType::NUMVECTORS;

        MRT(InitParams& initParams) :
                collisionMatrixDiagonals(MomentType::SetUpCollisionMatrix(initParams.lbmParams->GetTau()))
        {
        }

        void CalculateDensityMomentumFeq(VarsType& hydroVars, site_t index)
        {
            LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum,
                                                     hydroVars.velocity,
                                                     hydroVars.f_eq);

            for (unsigned int ii = 0; ii < NUMVECTORS; ++ii)
            {
                hydroVars.f_neq[ii] = hydroVars.f[ii] - hydroVars.f_eq[ii];
            }

            /** @todo #222 consider computing m_neq directly in the momentum space. See d'Humieres 2002. */
            ProjectVelsIntoMomentSpace(hydroVars.f_neq, hydroVars.m_neq);
        }

        void CalculateFeq(VarsType& hydroVars, site_t index)
        {
            LatticeType::CalculateFeq(hydroVars.density,
                                               hydroVars.momentum,
                                               hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }

            /** @todo #222 consider computing m_neq directly in the moment space. See d'Humieres 2002. */
            MomentType::ProjectVelsIntoMomentSpace(hydroVars.f_neq.f, hydroVars.m_neq);
          }

        void Collide(const LbmParameters* const lbmParams, VarsType& hydroVars)
        {
            for (Direction direction = 0; direction < NUMVECTORS; ++direction)
            {
              /** @todo #222 many optimisations possible (and necessary!).
               *  - Compute the loop below as a matrix product in DoCalculate*, alternatively we could consider reimplementing DoCollide to work with whole arrays (consider libraries boost::ublas or Armadillo)
               */
              distribn_t collision = 0.;
              for (unsigned momentIndex = 0; momentIndex < NUMMOMENTS;
                  momentIndex++)
              {
                collision += collisionMatrixDiagonals[momentIndex]
                             * normalisedReducedMomentBasis[momentIndex][direction]
                             * hydroVars.m_neq[momentIndex];
              }
              hydroVars.SetFPostCollision(direction, hydroVars.f[direction] - collision);
            }
          }

        /**
         *  This method is used in unit testing in order to make an MRT kernel behave as LBGK, regardless of the
         *  moment basis, by setting all the relaxation parameters to be the same.
         */
        void SetMrtRelaxationParameters(ConstDistSpan<NUMMOMENTS> newRelaxationParameters)
        {
            std::copy(newRelaxationParameters.begin(), newRelaxationParameters.end(), collisionMatrixDiagonals.begin());
        }

    private:
        /**
         * Projects a velocity distributions vector into the (reduced) MRT moment space.
         *
         * @param velDistributions velocity distributions vector
         * @param moments equivalent vector in the moment space
         */
        static void ProjectVelsIntoMomentSpace(ConstDistSpan<NUMVECTORS> velDistributions,
                                               MutDistSpan<NUMMOMENTS> moments)
        {
            for (unsigned momentIndex = 0; momentIndex < NUMMOMENTS; momentIndex++)
            {
                moments[momentIndex] = 0.;
                for (Direction velocityIndex = 0; velocityIndex < NUMVECTORS; velocityIndex++)
                {
                    moments[momentIndex] +=
                            MomentType::REDUCED_MOMENT_BASIS[momentIndex][velocityIndex] * velDistributions[velocityIndex];
                }
            }
        }

        /** MRT collision matrix (\hat{S}, diagonal). It corresponds to the inverse of the relaxation time for each mode. */
        std::array<distribn_t, NUMMOMENTS> collisionMatrixDiagonals;
        static constexpr MatrixType Normalise() {
            // Pre-compute the reduced moment basis divided by the basis times basis transposed.
            MatrixType ans;
            for (Direction direction = 0; direction < NUMVECTORS; ++direction)
            {
                for (unsigned momentIndex = 0; momentIndex < NUMMOMENTS; ++momentIndex)
                {
                    ans[momentIndex][direction] = MomentType::REDUCED_MOMENT_BASIS[momentIndex][direction]
                            / MomentType::BASIS_TIMES_BASIS_TRANSPOSED[momentIndex];
                }
            }
            return ans;
        }
        static constexpr MatrixType normalisedReducedMomentBasis = Normalise();
      };
}

#endif
