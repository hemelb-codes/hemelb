#ifndef HEMELB_LB_KERNELS_BASEKERNEL_H
#define HEMELB_LB_KERNELS_BASEKERNEL_H

#include <cstdlib>
#include "constants.h"
#include "lb/boundaries/BoundaryValues.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"

namespace hemelb
{
  namespace geometry{
    namespace neighbouring{
      class NeighbouringDataManager;
    }
  }
  namespace lb
  {
    namespace kernels
    {

      /**
       * HydroVars: struct for storing all of the hydrodynamics variables for doing the
       * colliding and streaming of Lattice-Boltzmann.
       *
       * This is specialised to each kernel implementation using a traits-type mechanism.
       * (Google for 'CRTP traits'). Each kernel implementation can also specialise this template.
       */

      template<class LatticeType>
      struct FVector
      {
        public:
          distribn_t f[LatticeType::NUMVECTORS];

          inline distribn_t& operator[](unsigned index)
          {
            return f[index];
          }

          inline const distribn_t& operator[](unsigned index) const
          {
            return f[index];
          }
      };

      template<class LatticeType>
      struct HydroVarsBase
      {
          template<class LatticeImpl> friend class Entropic;
          template<class LatticeImpl> friend class EntropicAnsumali;
          template<class LatticeImpl> friend class EntropicChik;
          template<class LatticeImpl> friend class LBGK;
          template<class rheologyModel, class LatticeImpl> friend class LBGKNN;
          template<class LatticeImpl> friend class MRT;

        protected:
          HydroVarsBase(const distribn_t* const f) :
              f(f)
          {
          }

        public:
          distribn_t density, v_x, v_y, v_z, tau;
          const distribn_t* const f;

          inline const FVector<LatticeType>& GetFEq()
          {
            return f_eq;
          }

          inline const FVector<LatticeType>& GetFNeq()
          {
            return f_neq;
          }

          // This is necessary as some of the streamers need the post-collision distribution.
          // It is calculated by collisions and kernels.
          inline void SetFPostCollision(Direction direction, distribn_t value)
          {
            fPostCollision[direction] = value;
          }

          inline const FVector<LatticeType>& GetFPostCollision()
          {
            return fPostCollision;
          }

        protected:
          FVector<LatticeType> f_eq, f_neq, fPostCollision;
      };

      template<typename KernelImpl>
      struct HydroVars : HydroVarsBase<typename KernelImpl::LatticeType>
      {
        public:
          HydroVars(const distribn_t* const f) :
              HydroVarsBase<typename KernelImpl::LatticeType>(f)
          {

          }
      };

      /**
       * InitParams: struct for passing variables into streaming, collision and kernel operators
       * to initialise them.
       *
       * When a newly-developed kernel, collider or streamer requires extra parameters to be
       * passed in for initialisation, it's annoying to have to change the constructors in
       * multiple places to make them all consistent (so that higher-up code can seamlessly
       * construct one kind or another).
       *
       * Instead, new parameters can be added to this single object, which should be the only
       * constructor argument used by any kernel / collision / streaming implementation.
       */
      struct InitParams
      {
        public:

          // Assume the first site to be used in the kernel is the first site in the core, unless otherwise specified
          InitParams()
          {
          }

          // The number of sites using this kernel instance.
          site_t siteCount;

          // The array with the imposed density at each boundary.
          boundaries::BoundaryValues* boundaryObject;

          // The lattice data object. Currently only used for accessing the boundary id
          // of each site next to an inlet or an outlet.
          const geometry::LatticeData* latDat;

          // The LB parameters object. Currently only used in LBGKNN to access the current
          // time step.
          const LbmParameters* lbmParams;
          geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager;
      };

      /**
       * BaseKernel: inheritable base class for the kernel. The public interface here define the
       * complete interface usable by collision operators:
       *  - Constructor(InitParams&)
       *  - KHydroVars, the type name for the kernel's hydrodynamic variable object.
       *  - LatticeType, the type of lattice being used (D3Q15, D3Q19 etc)
       *  - CalculateDensityVelocityFeq(KHydroVars&, site_t) for calculating
       *      the density, velocity and equilibrium distribution
       *  - Collide(const LbmParameters*, KHydroVars& hydroVars, unsigned int directionIndex)
       *  - Reset(InitParams*)
       *
       * The following must be implemented must be kernels (which derive from this class
       * using the CRTP).
       *  - Constructor(InitParams&)
       *  - DoCalculateDensityVelocityFeq(KHydroVars&, site_t)
       *  - DoCollide(const LbmParameters*, KHydroVars&, unsigned int) returns distibn_t
       *  - DoReset(InitParams*)
       */
      template<typename KernelImpl, typename LatticeImpl>
      class BaseKernel
      {
        public:
          typedef HydroVars<KernelImpl> KHydroVars;
          typedef LatticeImpl LatticeType;

          inline void CalculateDensityVelocityFeq(KHydroVars& hydroVars, site_t index)
          {
            static_cast<KernelImpl*>(this)->DoCalculateDensityVelocityFeq(hydroVars, index);
          }

          inline void CalculateFeq(KHydroVars& hydroVars, site_t index)
          {
            static_cast<KernelImpl*>(this)->DoCalculateFeq(hydroVars, index);
          }

          inline void Collide(const LbmParameters* lbmParams, KHydroVars& hydroVars)
          {
            static_cast<KernelImpl*>(this)->DoCollide(lbmParams, hydroVars);
          }

          inline void Reset(InitParams* init)
          {
            static_cast<KernelImpl*>(this)->DoReset(init);
          }
      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_BASEKERNEL_H */
