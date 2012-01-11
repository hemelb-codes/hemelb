#ifndef HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H
#define HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      template<typename KernelType>
      class NonZeroVelocityEquilibriumFixedDensity : public BaseCollision<
          NonZeroVelocityEquilibriumFixedDensity<KernelType>, KernelType>
      {
        public:
          typedef KernelType CKernel;

          NonZeroVelocityEquilibriumFixedDensity(kernels::InitParams& initParams) :
              kernel(initParams), boundaryObject(initParams.boundaryObject), latDat(initParams.latDat)
          {

          }

          inline void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                              const site_t index)
          {
            D3Q15::CalculateDensityAndVelocity(hydroVars.f,
                                               hydroVars.density,
                                               hydroVars.v_x,
                                               hydroVars.v_y,
                                               hydroVars.v_z);

            // Externally impose a density.
            hydroVars.density = boundaryObject->GetBoundaryDensity(latDat->GetBoundaryId(index));

            kernel.CalculateFeq(hydroVars, index);
          }

          inline void DoCollide(const LbmParameters* lbmParams,
                                kernels::HydroVars<KernelType>& iHydroVars)
          {
            for (Direction direction = 0; direction < D3Q15::NUMVECTORS; ++direction)
            {
              iHydroVars.GetFPostCollision()[direction] = iHydroVars.GetFEq()[direction];
            }
          }

          inline void DoReset(kernels::InitParams* init)
          {
            kernel.Reset(init);
          }

        private:
          KernelType kernel;
          boundaries::BoundaryValues* boundaryObject;
          const geometry::LatticeData* latDat;
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H */
