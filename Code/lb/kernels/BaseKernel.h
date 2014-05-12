// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_KERNELS_BASEKERNEL_H
#define HEMELB_LB_KERNELS_BASEKERNEL_H

#include <cstdlib>
#include "constants.h"
#include "lb/iolets/BoundaryValues.h"
#include "lb/kernels/rheologyModels/RheologyModels.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "util/Vector3D.h"

namespace hemelb
{
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
          HydroVarsBase(const distribn_t* const f, const util::Vector3D<distribn_t>* const force) :
            f(f), force(force)
          {
          }
          HydroVarsBase(const distribn_t* const f) :
                      f(f), force(NULL)
                    {
                    }

        public:
          distribn_t density, tau;
          util::Vector3D<distribn_t> momentum;
          util::Vector3D<distribn_t> velocity;

          const distribn_t* const f;
          // This is pointing to the vector-field of external forces
          // defined in LatticeData::forceAtSite
          const util::Vector3D<distribn_t>* const force;

          // Guo lattice distribution of external force contributions
          // as calculated in lattice::CalculateForceDistribution.
          inline const FVector<LatticeType>& GetForceDist() const
          {
            return forceDist;
          }
          inline void SetForceDist(Direction i, distribn_t val)
          {
             forceDist[i] = val;
          }

          inline const FVector<LatticeType>& GetFEq() const
          {
            return f_eq;
          }

          inline void SetFEq(Direction i, distribn_t val)
          {
            f_eq[i] = val;
          }

          inline distribn_t* GetFEqPtr()
          {
            return f_eq.f;
          }

          inline const FVector<LatticeType>& GetFNeq() const
          {
            return f_neq;
          }
          // This is necessary as some of the streamers need the post-collision distribution.
          // It is calculated by collisions and kernels.
          inline void SetFNeq(Direction direction, distribn_t value)
          {
            f_neq[direction] = value;
          }

          inline distribn_t* GetFNeqPtr()
          {
            return f_neq.f;
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
          FVector<LatticeType> f_eq, f_neq, fPostCollision, forceDist;
      };

      template<typename KernelImpl>
      struct HydroVars : HydroVarsBase<typename KernelImpl::LatticeType>
      {
        public:
          HydroVars(const distribn_t* const f, const util::Vector3D<distribn_t>* const force) :
            HydroVarsBase<typename KernelImpl::LatticeType> (f,force)
          {

          }

          HydroVars(const distribn_t* const f) :
                      HydroVarsBase<typename KernelImpl::LatticeType> (f)
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

          // Each streamer is responsible for updating certain types of sites. These are arranged such they are largely
          // contiguous in memory (the local contiguous site id). This data structure refers to which of those are handled
          // by the current streamer. These are given as a collection of contiguous site ids, running from e.g.
          // siteRanges[0].first to siteRanges[0].second-1 (inclusive).
          std::vector<std::pair<site_t, site_t> > siteRanges;

          // The array with the imposed density at each boundary.
          iolets::BoundaryValues* boundaryObject;

          // The lattice data object. Currently only used for accessing the boundary id
          // of each site next to an inlet or an outlet.
          const geometry::LatticeData* latDat;

          // The LB parameters object. Currently only used in LBGKNN to access the current
          // time step.
          const LbmParameters* lbmParams;

          // The neighbouring data manager, for kernels / collisions / streamers that
          // require data from other cores.
          geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager;
      };

      /**
       * BaseKernel: inheritable base class for the kernel. The public interface here define the
       * complete interface usable by collision operators:
       *  - Constructor(InitParams&)
       *  - KHydroVars, the type name for the kernel's hydrodynamic variable object.
       *  - LatticeType, the type of lattice being used (D3Q15, D3Q19 etc)
       *  - CalculateDensityMomentumFeq(KHydroVars&, site_t) for calculating
       *      the density, momentum and equilibrium distribution
       *  - Collide(const LbmParameters*, KHydroVars& hydroVars, unsigned int directionIndex)
       *  - Reset(InitParams*)
       *
       * The following must be implemented must be kernels (which derive from this class
       * using the CRTP).
       *  - Constructor(InitParams&)
       *  - DoCalculateDensityMomentumFeq(KHydroVars&, site_t)
       *  - DoCollide(const LbmParameters*, KHydroVars&, unsigned int) returns distibn_t
       *  - DoReset(InitParams*)
       */
      template<typename KernelImpl, typename LatticeImpl>
      class BaseKernel
      {
        public:
          typedef HydroVars<KernelImpl> KHydroVars;
          typedef LatticeImpl LatticeType;

          inline void CalculateDensityMomentumFeq(KHydroVars& hydroVars, site_t index)
          {
            static_cast<KernelImpl*> (this)->DoCalculateDensityMomentumFeq(hydroVars, index);
          }

          inline void CalculateFeq(KHydroVars& hydroVars, site_t index)
          {
            static_cast<KernelImpl*> (this)->DoCalculateFeq(hydroVars, index);
          }

          inline void Collide(const LbmParameters* lbmParams, KHydroVars& hydroVars)
          {
            static_cast<KernelImpl*> (this)->DoCollide(lbmParams, hydroVars);
          }

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_BASEKERNEL_H */
