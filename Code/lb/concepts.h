// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_CONCEPTS_H
#define HEMELB_LB_CONCEPTS_H

#include <concepts>
#include "lb/LbmParameters.h"
#include "lb/lattices/Lattice.h"

namespace hemelb::geometry {
        class FieldData;
}
namespace hemelb::lb {
    class MacroscopicPropertyCache;

    namespace detail {
        template<typename Base, typename Derived>
        concept base_of = std::derived_from<Derived, Base>;
    }

    // A lattice is an instantiation of Lattice (or a type derived from such)
    template<typename L>
    concept lattice_type = requires(L l) {
        { Lattice{l}} -> detail::base_of<L>;
    };

    template <typename K>
    concept kernel_type = requires {
        typename K::LatticeType;
        requires lattice_type<typename K::LatticeType>;
        typename K::VarsType;
        requires requires (K k, typename K::VarsType v, site_t idx, LbmParameters const* lbmParams) {
            { k.CalculateDensityMomentumFeq(v, idx) };
            { k.Collide(lbmParams, v)};
        };
    };

    // For MRT kernels
    template<typename MB>
    concept moment_basis = requires {
        // A lattice
        requires lattice_type<typename MB::Lattice>;
        // That the number of kinetic/non-hydrodynamic moments be consistent with the Lattice
        // (as there are 4 conserved moments in 3D).
        typename std::enable_if<MB::NUMMOMENTS + 4 == MB::Lattice::NUMVECTORS>::type;
        typename MB::MatrixType;
        // That there be a static member to compute the (diagonal) collision matrix
        { MB::SetUpCollisionMatrix(0.0) } -> std::same_as<std::array<distribn_t, MB::NUMMOMENTS>>;
    };

    template<typename S>
    concept streamer = requires(S s,
            InitParams& i,
            site_t site_idx, LbmParameters const* lbmParameters,
            geometry::FieldData& dom, MacroscopicPropertyCache& cache
            ) {
        // public typedef of collision
        typename S::CollisionType;
        // Constructor accepting InitParams&
        { S{i} };
        // Main work functions!
        { s.StreamAndCollide(site_idx, site_idx, lbmParameters, dom, cache) };
        { s.PostStep(site_idx, site_idx, lbmParameters, dom, cache) };
    };
}
#endif