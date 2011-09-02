#ifndef HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUM_H
#define HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUM_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      /**
       * Collision operator that maintains the density, set velocity to 0 and fully relaxes
       * to equilibrium.
       */
      template<typename KernelType>
      class ZeroVelocityEquilibrium : public BaseCollision<ZeroVelocityEquilibrium<KernelType>,
          KernelType>
      {
        public:
          typedef KernelType CKernel;

          ZeroVelocityEquilibrium(kernels::InitParams& initParams) :
              BaseCollision<ZeroVelocityEquilibrium<KernelType>, KernelType>(), kernel(initParams)
          {

          }

          void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                       const site_t index)
          {
            hydroVars.density = 0.0;

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ++ii)
            {
              hydroVars.density += hydroVars.f[ii];
            }

            hydroVars.v_x = 0.0;
            hydroVars.v_y = 0.0;
            hydroVars.v_z = 0.0;

            kernel.CalculateFeq(hydroVars, index);
          }

          distribn_t DoCollide(const LbmParameters* lbmParams,
                               unsigned int directionIndex,
                               kernels::HydroVars<KernelType>& iHydroVars)
          {
            return iHydroVars.f_eq[directionIndex];
          }

          void DoReset(kernels::InitParams* init)
          {
            kernel.Reset(init);
          }

        private:
          KernelType kernel;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUM_H */
