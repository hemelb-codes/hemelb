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
      class MRT;

      template<>
      struct HydroVars<MRT> : public HydroVarsBase
      {
        public:
          HydroVars(const distribn_t* const f) :
              HydroVarsBase(f)
          {
          }

          distribn_t m_neq[D3Q15::NUM_KINETIC_MOMENTS];
      };

      class MRT : public BaseKernel<MRT>
      {
        public:

          MRT(InitParams& initParams):
            collisionMatrix(initParams.lbmParams->GetMrtRelaxationParameters())
          {
          }

          void DoCalculateDensityVelocityFeq(HydroVars<MRT>& hydroVars, site_t index)
          {
            D3Q15::CalculateDensityVelocityFEq(hydroVars.f,
                                               hydroVars.density,
                                               hydroVars.v_x,
                                               hydroVars.v_y,
                                               hydroVars.v_z,
                                               hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }

            /** @todo #61 consider computing m_neq directly in the momentum space. See d'Humieres 2002. */
            D3Q15::ProjectVelsIntoMomentSpace(hydroVars.f_neq.f, hydroVars.m_neq);
          }

          void DoCalculateFeq(HydroVars<MRT>& hydroVars, site_t index)
          {
            D3Q15::CalculateFeq(hydroVars.density,
                                hydroVars.v_x,
                                hydroVars.v_y,
                                hydroVars.v_z,
                                hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }

            /** @todo #61 consider computing m_neq directly in the momentum space. See d'Humieres 2002. */
            D3Q15::ProjectVelsIntoMomentSpace(hydroVars.f_neq.f, hydroVars.m_neq);
          }

          distribn_t DoCollide(const LbmParameters* const lbmParams,
                               HydroVars<MRT>& hydroVars,
                               unsigned int direction)
          {
            /** @todo #61 many optimisations possible (and necessary!).
             *  - Store the product of REDUCED_MOMENT_BASIS and 1/BASIS_TIMES_BASIS_TRANSPOSED instead of REDUCED_MOMENT_BASIS
             *  - Compute the loop below as a matrix product in DoCalculate*, alternatively we could consider reimplementing DoCollide to work with whole arrays (consider libraries boost::ublas or Armadillo)
             */
            distribn_t collision = 0.;
            for (unsigned momentIndex = 0; momentIndex < D3Q15::NUM_KINETIC_MOMENTS; momentIndex++)
            {
              collision += (collisionMatrix[momentIndex]
                  / D3Q15::BASIS_TIMES_BASIS_TRANSPOSED[momentIndex])
                  * D3Q15::REDUCED_MOMENT_BASIS[momentIndex][direction]
                  * hydroVars.m_neq[momentIndex];
            }

            return hydroVars.f[direction] - collision;
          }

          void DoReset(InitParams* initParams)
          {
          }

        private:
          const std::vector<distribn_t>& collisionMatrix;

      };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_MRT_H */
