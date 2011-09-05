#ifndef HEMELB_LB_COLLISIONS_NORMAL_H
#define HEMELB_LB_COLLISIONS_NORMAL_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      /**
       * Normal collisions - use a formula to relax the distribution towards equilibrium.
       */
      template<typename KernelType>
      class Normal : public BaseCollision<Normal<KernelType>, KernelType>
      {
        public:
          typedef KernelType CKernel;

          Normal(kernels::InitParams& initParams) :
              kernel(initParams)
          {
          }

          void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                       const site_t index)
          {
            kernel.CalculateDensityVelocityFeq(hydroVars, index);
          }

          distribn_t DoCollide(const LbmParameters* lbmParams,
                               unsigned int directionIndex,
                               kernels::HydroVars<KernelType>& iHydroVars)
          {
            return kernel.Collide(lbmParams, iHydroVars, directionIndex);
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

#endif /* HEMELB_LB_COLLISIONS_NORMAL_H */
