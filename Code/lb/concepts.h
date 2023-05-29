// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_CONCEPTS_H
#define HEMELB_LB_CONCEPTS_H

#include <concepts>
#include "lb/LbmParameters.h"
#include "lb/lattices/Lattice.h"
#include "util/concepts.h"

namespace hemelb::geometry {
    template <typename>
    class Site;
    class Domain;
    class FieldData;
}
namespace hemelb::lb {
    class MacroscopicPropertyCache;

    // Would like the template parameter to be kernel_type, but this would be circular.
    template <typename>
    struct HydroVars;

    // A lattice is an instantiation of Lattice (or a type derived from such)
    template<typename L>
    concept lattice_type = requires(L l) {
        { Lattice{l}} -> util::base_of<L>;
    };

    template <typename K>
    concept kernel_type =
    lattice_type<typename K::LatticeType> &&
    requires (InitParams& i, K k, typename K::VarsType v, site_t idx, LbmParameters const* lbmParams) {
        typename K::VarsType;
        { K(i) };
        { k.CalculateDensityMomentumFeq(v, idx) };
        { k.Collide(lbmParams, v)};
    };

    // For MRT kernels
    template<typename MB>
    concept moment_basis =
    lattice_type<typename MB::Lattice> && // A lattice
    requires {
        // That the number of kinetic/non-hydrodynamic moments be consistent with the Lattice
        // (as there are 4 conserved moments in 3D).
        typename std::enable_if<MB::NUMMOMENTS + 4 == MB::Lattice::NUMVECTORS>::type;
        typename MB::MatrixType;
        // That there be a static member to compute the (diagonal) collision matrix
        { MB::SetUpCollisionMatrix(0.0) } -> std::same_as<std::array<distribn_t, MB::NUMMOMENTS>>;
    };

    // Collision - basically a wrapper around a kernel
    template <typename C>
    concept collision_type =
    kernel_type<typename C::KernelType> &&
    lattice_type<typename C::LatticeType> &&
    requires (
            C c,
            InitParams& i,
            typename C::VarsType v,
            geometry::Site<geometry::Domain> const site,
            LbmParameters const* lbmParams
    ) {
        { C(i) };
        { c.CalculatePreCollision(v, site) };
        { c.Collide(lbmParams, v) };
    };

    /// Concept for doing collide-and-stream on one link (lattice vector)
    template <typename T>
    concept link_streamer =
    collision_type<typename T::CollisionType> &&
    lattice_type<typename T::LatticeType> &&
    requires(
            T& linkStreamer,
            LbmParameters const* lbmParams,
            geometry::FieldData& data,
            geometry::Site<geometry::FieldData> const& site,
            typename T::CollisionType::VarsType& hydroVars,
            Direction d
    ) {
        { linkStreamer.StreamLink(lbmParams, data, site, hydroVars, d) };
        { linkStreamer.PostStepLink(data, site, d) };
    };
    template<typename S>
    concept streamer =
    collision_type<typename S::CollisionType> && // collision
    requires(S s,
             InitParams& i,
             site_t site_idx, LbmParameters const* lbmParameters,
             geometry::FieldData& dom, MacroscopicPropertyCache& cache
    ) {
        // Constructor accepting InitParams&
        { S{i} };
        // Main work functions!
        { s.StreamAndCollide(site_idx, site_idx, lbmParameters, dom, cache) };
        { s.PostStep(site_idx, site_idx, lbmParameters, dom, cache) };
    };
}
#endif