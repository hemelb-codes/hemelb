// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_HYDROVARS_H
#define HEMELB_LB_HYDROVARS_H

#include <array>

#include "units.h"
#include "lb/concepts.h"
#include "geometry/Site.h"
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

}

#endif
