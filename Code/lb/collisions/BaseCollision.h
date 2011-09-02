#ifndef HEMELB_LB_COLLISIONS_BASECOLLISION_H
#define HEMELB_LB_COLLISIONS_BASECOLLISION_H

#include "constants.h"
#include "lb/LbmParameters.h"
#include "D3Q15.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      /**
       * BaseCollision: inheritable base class for the collision operator. The public interface
       * here defines the complete interface usable by streaming operators:
       *  - CalculatePreCollision(CHydroVars&, site_t)
       *  - Collide(const LbmParameters*, unsigned int, CHydroVars&)
       *  - Reset(InitParams*)
       *
       * The following must be implemented must be collisions (which derive from this class
       * using the CRTP).
       *  - CKernel, the type of the kernel to be used.
       *  - DoCalculatePreCollision(CHydroVars&, site_t)
       *  - DoCollide(const LbmParameters*, unsigned int, CHydroVars&) returns distribn_t
       *  - DoReset(InitParams*)
       */
      template<typename CollisionImpl, typename KernelImpl>
      class BaseCollision
      {
        public:
          typedef KernelImpl CKernel;

          void CalculatePreCollision(kernels::HydroVars<KernelImpl>& hydroVars, site_t index)
          {
            static_cast<CollisionImpl*>(this)->DoCalculatePreCollision(hydroVars, index);
          }

          distribn_t Collide(const LbmParameters* lbmParams,
                             unsigned int directionIndex,
                             kernels::HydroVars<KernelImpl>& hydroVars)
          {
            return static_cast<CollisionImpl*>(this)->DoCollide(lbmParams,
                                                                directionIndex,
                                                                hydroVars);
          }

          void Reset(kernels::InitParams* init)
          {
            static_cast<CollisionImpl*>(this)->DoReset(init);
          }
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_BASECOLLISION_H */
