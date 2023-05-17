// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_FIELDDATA_H
#define HEMELB_GEOMETRY_FIELDDATA_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <vector>

#include "units.h"
#include "geometry/Domain.h"
#include "geometry/Site.h"
#include "geometry/neighbouring/NeighbouringDomain.h"
#include "util/Vector3D.h"

namespace hemelb::net { class Net; }
namespace hemelb::tests::helpers { class LatticeDataAccess; }

namespace hemelb::geometry {

    // Hold field data across a geometry, as described by a domain_type.
    class FieldData {
    public:
        friend class tests::helpers::LatticeDataAccess;

        using domain_type = Domain;
    protected:
        std::shared_ptr <domain_type> m_domain;
        // For now just, list our fields.
        std::vector <distribn_t> m_currentDistributions; //! The distribution values at the start of the current time step.
        std::vector <distribn_t> m_nextDistributions; //! The distribution values for the next time step.
        std::vector <LatticeForceVector> m_force; //! Holds the force vector at a fluid site

        std::unique_ptr <neighbouring::NeighbouringFieldData> m_neighbouringFields;

        static std::size_t CalcDistSize(Domain const &d);

    public:
        FieldData() = default;

        FieldData(std::shared_ptr <domain_type> d);

        // Access the domain
        inline domain_type const &GetDomain() const {
            return *m_domain.get();
        }

        inline domain_type &GetDomain() {
            return *m_domain.get();
        }

        inline neighbouring::NeighbouringFieldData &GetNeighbouringData() {
            return *m_neighbouringFields;
        }

        inline neighbouring::NeighbouringFieldData const &GetNeighbouringData() const {
            return *m_neighbouringFields;
        }

        /**
         * Get a site object for the given index.
         * @param localIndex
         * @return
         */
        inline Site <FieldData> GetSite(site_t localIndex) {
            return Site<FieldData>(localIndex, *this);
        }

        /**
         * Get a site object for the given position.
         */
        inline Site <FieldData> GetSite(LatticeVector const &pos) {
            return GetSite(m_domain->GetContiguousSiteId(pos));
        }

        /**
         * Get a site object for the given position.
         */
        inline Site<const FieldData> GetSite(LatticeVector const &pos) const {
            return GetSite(m_domain->GetContiguousSiteId(pos));
        }

        /**
         * Get a site object for the given position.
         */
        inline Site <FieldData> GetSite(site_t x, site_t y, site_t z) {
            return GetSite(LatticeVector(x, y, z));
        }

        /**
         * Get a site object for the given position.
         */
        inline Site<const FieldData> GetSite(site_t x, site_t y, site_t z) const {
            return GetSite(LatticeVector(x, y, z));
        }

        /**
         * Get a const site object for the given index.
         * @param localIndex
         * @return
         */
        inline Site<const FieldData> GetSite(site_t localIndex) const {
            return Site<const FieldData>(localIndex, *this);
        }

        /**
         * Get a pointer to the fOld array starting at the requested index
         * @param distributionIndex
         * @return
         */
        // Method should remain protected, intent is to access this information via Site
        inline distribn_t *GetFOld(site_t distributionIndex) {
            return &m_currentDistributions[distributionIndex];
        }

        /**
         * Get a pointer to the fOld array starting at the requested index. This version
         * of the function allows us to access the fOld array in a const way from a const
         * domain_type
         * @param distributionIndex
         * @return
         */
        // Method should remain protected, intent is to access this information via Site
        inline const distribn_t *GetFOld(site_t distributionIndex) const {
            return &m_currentDistributions[distributionIndex];
        }

        /**
         * Get a pointer into the fNew array at the given index
         * @param distributionIndex
         * @return
         */
        inline distribn_t *GetFNew(site_t distributionIndex) {
            return &m_nextDistributions[distributionIndex];
        }

        template <typename LatticeType>
        auto GetFNew(site_t site_idx) {
            constexpr auto Q = LatticeType::NUMVECTORS;
            return MutDistSpan<Q>{&m_nextDistributions[site_idx * Q], Q};
        }

        /**
         * Get a pointer into the fNew array at the given index. This version of the above lets us
         * use a const version of a domain_type to get a const *.
         * @param distributionIndex
         * @return
         */
        inline const distribn_t *GetFNew(site_t distributionIndex) const {
            return &m_nextDistributions[distributionIndex];
        }

        template <typename LatticeType>
        auto GetFNew(site_t site_idx) const {
            constexpr auto Q = LatticeType::NUMVECTORS;
            return ConstDistSpan<Q>{&m_nextDistributions[site_idx * Q], Q};
        }

        //! Swap the fOld and fNew arrays around.
        inline void SwapOldAndNew() {
            m_currentDistributions.swap(m_nextDistributions);
        }

        //! Reset forces to some constant value
        //! Could be zero. Or could be something like a constant term for gravity.
        inline void ResetForces(LatticeForceVector const &force = LatticeForceVector(0, 0, 0)) {
            std::fill(m_force.begin(), m_force.end(), force);
        }

        // Method should remain protected, intent is to access this information via Site
        inline LatticeForceVector const &GetForceAtSite(site_t iSiteIndex) const {
            return m_force[iSiteIndex];
        }
        /**
         * Set the force vector at the given site
         * @param iSiteIndex
         * @param force
         * @return
         */
        // Method should remain protected, intent is to set this information via Site
        inline void SetForceAtSite(site_t iSiteIndex, LatticeForceVector const &force) {
            assert(iSiteIndex >= site_t(0));
            assert(m_force.size() > std::size_t(iSiteIndex));
            m_force[iSiteIndex] = force;
        }

        inline void AddToForceAtSite(site_t iSiteIndex, LatticeForceVector const &force) {
            assert(iSiteIndex >= site_t(0));
            assert(m_force.size() > std::size_t(iSiteIndex));
            m_force[iSiteIndex] += force;
        }
        /**
         * Set a vertical force vector at the given site
         * @param iSiteIndex
         * @param force
         * @return
         */
        // Method should remain protected, intent is to set this information via Site
        inline void SetForceAtSite(site_t iSiteIndex, LatticeForce force) {
            assert(iSiteIndex >= site_t(0));
            assert(m_force.size() > std::size_t(iSiteIndex));
            m_force[iSiteIndex] = util::Vector3D<distribn_t>(0.0, 0.0, force);
        }

        void SendAndReceive(net::Net *net);

        void CopyReceived();

    };
}
#endif // once