#ifndef HEMELB_LB_KERNELS_LBGK_H
#define HEMELB_LB_KERNELS_LBGK_H

#include <cstdlib>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      /**
       * LBGK: This class implements the LBGK single-relaxation time kernel.
       */
      class LBGK : public BaseKernel<LBGK>
      {
        public:
          LBGK(InitParams& initParams)
          {
          }

          inline void DoCalculateDensityVelocityFeq(HydroVars<LBGK>& hydroVars, site_t index)
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
          }

          inline void DoCalculateFeq(HydroVars<LBGK>& hydroVars, site_t index)
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
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<LBGK>& hydroVars)
          {
            for (Direction direction = 0; direction < D3Q15::NUMVECTORS; ++direction)
            {
              hydroVars.GetFPostCollision()[direction] = hydroVars.f[direction]
                  + hydroVars.f_neq.f[direction] * lbmParams->GetOmega();
            }
          }

          inline void DoReset(InitParams* init)
          {

          }
      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_LBGK_H */
