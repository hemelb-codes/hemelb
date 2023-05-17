// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_BASEKERNEL_H
#define HEMELB_LB_KERNELS_BASEKERNEL_H

#include <cstdlib>
#include "constants.h"
#include "lb/concepts.h"
#include "lb/lattices/Lattice.h"
#include "lb/iolets/BoundaryValues.h"
#include "lb/kernels/RheologyModels.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "geometry/SiteData.h"
#include "util/Vector3D.h"

namespace hemelb::lb
{

    template <class LatticeType>
    using FVector = std::array<distribn_t, LatticeType::NUMVECTORS>;

    /**
     * HydroVars: struct for storing all of the hydrodynamics variables for doing the
     * colliding and streaming of Lattice-Boltzmann.
     *
     * This base is needed by all kernels. The primary template class should be used
     * and can be specialised for those kernels that need extra data/behaviour.
     */
    template<lattice_type LatticeType>
    struct HydroVarsBase
    {
        using const_span = typename LatticeType::const_span;
        using mut_span = typename LatticeType::mut_span;

        template<lattice_type> friend class EntropicBase;
        template<lattice_type> friend class EntropicAnsumali;
        template<lattice_type> friend class EntropicChik;
        template<lattice_type> friend class LBGK;
        template<class rheologyModel, lattice_type> friend class LBGKNN;
        template<moment_basis> friend class MRT;
        template<lattice_type> friend class TRT;

        HydroVarsBase(const_span s) : f(s) {
        }

        HydroVarsBase(const distribn_t* const f) :
                f(f, LatticeType::NUMVECTORS)
        {
        }
        template<class DataSource>
        HydroVarsBase(geometry::Site<DataSource> const &_site) :
                f(_site.template GetFOld<LatticeType>())
        {
        }

        distribn_t density, tau;
        util::Vector3D<distribn_t> momentum;
        util::Vector3D<distribn_t> velocity;

        const_span f;

        mut_span GetFEq()
        {
            return f_eq;
        }
        const_span GetFEq() const
        {
            return f_eq;
        }

        void SetFEq(Direction i, distribn_t val)
        {
            f_eq[i] = val;
        }

        distribn_t* GetFEqPtr()
        {
            return f_eq.f;
        }

        const FVector<LatticeType>& GetFNeq() const
        {
            return f_neq;
        }

        // This is necessary as some of the streamers need the post-collision distribution.
        // It is calculated by collisions and kernels.
        void SetFNeq(Direction direction, distribn_t value)
        {
            f_neq[direction] = value;
        }

        distribn_t* GetFNeqPtr()
        {
            return f_neq.data();
        }

        // This is necessary as some of the streamers need the post-collision distribution.
        // It is calculated by collisions and kernels.
        void SetFPostCollision(Direction direction, distribn_t value)
        {
            fPostCollision[direction] = value;
        }

        const FVector<LatticeType>& GetFPostCollision()
        {
            return fPostCollision;
        }

    protected:
        FVector<LatticeType> f_eq, f_neq, fPostCollision;
    };

    /// This exists to be specialised for some kernels that need extra members.
    template <typename KernelImpl>
    struct HydroVars : HydroVarsBase<typename KernelImpl::LatticeType>
    {
        using HydroVarsBase<typename KernelImpl::LatticeType>::HydroVarsBase;
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
          InitParams() = default;

          // The number of sites using this kernel instance.
          site_t siteCount;

          // Each streamer is responsible for updating certain types of sites. These are arranged such they are largely
          // contiguous in memory (the local contiguous site id). This data structure refers to which of those are handled
          // by the current streamer. These are given as a collection of contiguous site ids, running from e.g.
          // siteRanges[0].first to siteRanges[0].second-1 (inclusive).
          std::vector<std::pair<site_t, site_t> > siteRanges;

          // The array with the imposed density at each boundary.
          BoundaryValues* boundaryObject;

          // The lattice data object. Currently only used for accessing the boundary id
          // of each site next to an inlet or an outlet.
          const geometry::Domain* latDat;

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
    template <typename KernelImpl, lattice_type LatticeImpl>
    class BaseKernel
    {
    public:
        using KHydroVars = HydroVars<KernelImpl>;
        using LatticeType = LatticeImpl;

        inline void CalculateDensityMomentumFeq(KHydroVars& hydroVars, site_t index)
        {
            static_cast<KernelImpl*>(this)->DoCalculateDensityMomentumFeq(hydroVars, index);
        }

        inline void CalculateFeq(KHydroVars& hydroVars, site_t index)
        {
            static_cast<KernelImpl*>(this)->DoCalculateFeq(hydroVars, index);
        }

        inline void Collide(const LbmParameters* lbmParams, KHydroVars& hydroVars)
        {
            static_cast<KernelImpl*>(this)->DoCollide(lbmParams, hydroVars);
        }

    };

}

#endif /* HEMELB_LB_KERNELS_BASEKERNEL_H */
