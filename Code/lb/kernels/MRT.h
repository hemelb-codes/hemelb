#ifndef HEMELB_LB_STREAMERS_MRT_H
#define HEMELB_LB_STREAMERS_MRT_H

#include "lb/kernels/BaseKernel.h"
#include "lb/SimulationState.h"
#include <cassert>
#include <cmath>

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      // Forward declaration needed by the struct
      template<class LatticeType> class MRT;

      template<class LatticeType>
      struct HydroVars<MRT<LatticeType> > : public HydroVarsBase<LatticeType>
      {
        public:
          HydroVars(const distribn_t* const f) :
              HydroVarsBase<LatticeType>(f)
          {
          }

          /** Equilibrium velocity distribution in the momentum space. */
          distribn_t m_neq[LatticeType::NUM_KINETIC_MOMENTS];
      };

      /**
       * This class implements the Multiple Relaxation Time (MRT) collision operator.
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
      template<class LatticeType>
      class MRT : public BaseKernel<MRT<LatticeType>, LatticeType>
      {
        public:

          MRT(InitParams& initParams) :
              collisionMatrix(initParams.lbmParams->GetMrtRelaxationParameters())
          {
          }

          inline void DoCalculateDensityVelocityFeq(HydroVars<MRT>& hydroVars, site_t index)
          {
            LatticeType::CalculateDensityVelocityFEq(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.v_x,
                                                     hydroVars.v_y,
                                                     hydroVars.v_z,
                                                     hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }

            /** @todo #61 consider computing m_neq directly in the momentum space. See d'Humieres 2002. */
            LatticeType::ProjectVelsIntoMomentSpace(hydroVars.f_neq.f, hydroVars.m_neq);
          }

          inline void DoCalculateFeq(HydroVars<MRT>& hydroVars, site_t index)
          {
            LatticeType::CalculateFeq(hydroVars.density, hydroVars.v_x, hydroVars.v_y, hydroVars.v_z, hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }

            /** @todo #61 consider computing m_neq directly in the momentum space. See d'Humieres 2002. */
            LatticeType::ProjectVelsIntoMomentSpace(hydroVars.f_neq.f, hydroVars.m_neq);
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<MRT>& hydroVars)
          {
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
              /** @todo #61 many optimisations possible (and necessary!).
               *  - Store the product of REDUCED_MOMENT_BASIS and 1/BASIS_TIMES_BASIS_TRANSPOSED instead of REDUCED_MOMENT_BASIS
               *  - Compute the loop below as a matrix product in DoCalculate*, alternatively we could consider reimplementing DoCollide to work with whole arrays (consider libraries boost::ublas or Armadillo)
               */
              distribn_t collision = 0.;
              for (unsigned momentIndex = 0; momentIndex < LatticeType::NUM_KINETIC_MOMENTS; momentIndex++)
              {
                collision += (collisionMatrix[momentIndex] / LatticeType::BASIS_TIMES_BASIS_TRANSPOSED[momentIndex])
                    * LatticeType::REDUCED_MOMENT_BASIS[momentIndex][direction] * hydroVars.m_neq[momentIndex];
              }
              hydroVars.GetFPostCollision()[direction] = hydroVars.f[direction] - collision;
            }
          }

          inline void DoReset(InitParams* initParams)
          {
          }

        private:
          /** MRT collision matrix (\hat{S}, diagonal). It corresponds to the inverse of the relaxation time for each mode. */
          const std::vector<distribn_t>& collisionMatrix;

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_MRT_H */
