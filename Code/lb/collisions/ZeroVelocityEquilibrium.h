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
      class ZeroVelocityEquilibrium : public BaseCollision<ZeroVelocityEquilibrium<KernelType>, KernelType>
      {
        public:
          typedef KernelType CKernel;

          ZeroVelocityEquilibrium(kernels::InitParams& initParams) :
              BaseCollision<ZeroVelocityEquilibrium<KernelType>, KernelType>(), kernel(initParams)
          {

          }

          inline void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars, const geometry::Site& site)
          {
            hydroVars.density = 0.0;

            for (unsigned int ii = 0; ii < CKernel::LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.density += hydroVars.f[ii];
            }

            hydroVars.momentum = util::Vector3D<distribn_t>::Zero();

            kernel.CalculateFeq(hydroVars, site.GetIndex());
          }

          inline void DoCollide(const LbmParameters* lbmParams, kernels::HydroVars<KernelType>& iHydroVars)
          {
            for (Direction direction = 0; direction < CKernel::LatticeType::NUMVECTORS; ++direction)
            {
              iHydroVars.SetFPostCollision(direction, iHydroVars.GetFEq()[direction]);
            }
          }

          inline void DoReset(kernels::InitParams* init)
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
