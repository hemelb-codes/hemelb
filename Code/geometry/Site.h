// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_SITE_H
#define HEMELB_GEOMETRY_SITE_H

#include <span>

#include "units.h"
#include "geometry/SiteData.h"
#include "util/Vector3D.h"

namespace hemelb::geometry
{
    namespace detail {
        // Implement a simple trait for determining if a type has an inner type "domain_type".
        template <typename T, typename = void>
        struct has_domain_type : std::false_type {};
        // Specialisation for true
        template <typename T>
        struct has_domain_type<T, std::void_t<typename T::domain_type>> : std::true_type {};
        // Helper
        template <typename T>
        inline constexpr bool has_domain_v = has_domain_type<T>::value;

        template <typename SOURCE, typename X>
        using conditionally_add_const_t = typename std::conditional<
                std::is_const_v<SOURCE>,
                std::add_const_t<X>,
                std::remove_const_t<X>
        >::type;

        template <typename DS, typename = void>
        struct SiteDataSourceTraits {
            using data_source_type = DS;
            static constexpr bool data_source_has_domain = false;

            using domain_type = DS;
            using field_type = void*;

            static domain_type& GetDomain(DS& ds) {
                return ds;
            }
            static field_type* GetField(DS& ds) {
                return nullptr;
            }
        };

        template <typename DS>
        struct SiteDataSourceTraits<DS, std::void_t<typename DS::domain_type>> {
            using data_source_type = DS;
            static constexpr bool data_source_has_domain = true;

            using domain_type = conditionally_add_const_t<DS, typename DS::domain_type>;
            using field_type = DS;

            static domain_type& GetDomain(DS &ds) {
                return ds.GetDomain();
            }

            static field_type *GetField(DS &ds) {
                return &ds;
            }
        };
    }

    // There are two types of data about a site: the geometrical
    // information and the field information. The former comes from
    // a domain object and the latter from a field data.
    //
    // Template parameter can be either of these, but obviously it
    // will be an error to try and access field data if instantiated
    // on a domain_type-like type.
    template<class DataSource>
    class Site
    {
    public:
        using traits = detail::SiteDataSourceTraits<DataSource>;
        using domain_type = typename traits::domain_type;
        using field_type = typename traits::field_type;

    private:
        // Members
        site_t index;
        field_type* m_fieldData;
        domain_type* m_domain;

    public:
        Site(site_t localContiguousIndex, DataSource &latticeData) :
                index{localContiguousIndex}, m_fieldData{traits::GetField(latticeData)}, m_domain{&traits::GetDomain(latticeData)}
        {
        }

        // Allow conversion from site<field> -> site<domain>
        operator Site<domain_type>() const {
                return Site<domain_type>{index, *m_domain};
        }

        bool IsWall() const
        {
          return GetSiteData().IsWall();
        }

        bool IsSolid() const
        {
          return GetSiteData().IsSolid();
        }

        unsigned GetCollisionType() const
        {
          return GetSiteData().GetCollisionType();
        }

        SiteType GetSiteType() const
        {
          return GetSiteData().GetSiteType();
        }

        int GetIoletId() const
        {
          return GetSiteData().GetIoletId();
        }

        bool HasWall(Direction direction) const
        {
          return GetSiteData().HasWall(direction);
        }

        bool HasIolet(Direction direction) const
        {
          return GetSiteData().HasIolet(direction);
        }

        template<typename LatticeType>
        distribn_t GetWallDistance(Direction direction) const
        {
          return m_domain->template GetCutDistance<LatticeType>(index, direction);
        }

        auto* GetWallDistances()
        {
          return m_domain->GetCutDistances(index);
        }

        const distribn_t* GetWallDistances() const
        {
          return m_domain->GetCutDistances(index);
        }

        util::Vector3D<distribn_t> const& GetWallNormal() const
        {
          return m_domain->GetNormalToWall(index);
        }

        // Return Vector3D<distribn_t>& may be const qualified if domain_type is const
        auto& GetWallNormal()
        {
          return m_domain->GetNormalToWall(index);
        }
        const LatticeForceVector& GetForce() const
        {
          return m_fieldData->GetForceAtSite(index);
        }
        void SetForce(LatticeForceVector const &_force)
        {
          return m_fieldData->SetForceAtSite(index, _force);
        }
        void AddToForce(LatticeForceVector const &_force)
        {
          return m_fieldData->AddToForceAtSite(index, _force);
        }
        site_t GetIndex() const
        {
          return index;
        }

        /**
         * This returns the index of the distribution to stream to. If streaming would take the
         * distribution out of the geometry, we instead stream to the 'rubbish site', an extra
         * position in the array that doesn't correspond to any site in the geometry.
         *
         * @param direction
         * @return
         */
        template<typename LatticeType>
        site_t GetStreamedIndex(Direction direction) const
        {
          return m_domain->template GetStreamedIndex<LatticeType>(index, direction);
        }

        template<typename LatticeType>
        auto GetFOld() const
        {
            return ConstDistSpan<LatticeType::NUMVECTORS>{m_fieldData->GetFOld(index * LatticeType::NUMVECTORS), LatticeType::NUMVECTORS};
        }

        // Non-templated version of GetFOld, for when you haven't got a lattice type handy
        const distribn_t* GetFOld(int numvectors) const
        {
          return m_fieldData->GetFOld(index * numvectors);
        }

        // Note that for const qualified field_type, dists are const too.
        template<typename LatticeType>
        auto GetFOld()
        {
            auto ptr = m_fieldData->GetFOld(index * LatticeType::NUMVECTORS);
            // To correctly return the Const/Mut span
            return std::span<
                    typename std::pointer_traits<decltype(ptr)>::element_type,
                    LatticeType::NUMVECTORS
            >{ptr, LatticeType::NUMVECTORS};
        }

        // Non-templated version of GetFOld, for when you haven't got a lattice type handy
        auto GetFOld(int numvectors)
        {
            return m_fieldData->GetFOld(index * numvectors);
        }

        const SiteData& GetSiteData() const
        {
          return m_domain->GetSiteData(index);
        }

        // Note use auto to deduce correct const/non-const of SiteDate return type.
        auto& GetSiteData()
        {
          return m_domain->GetSiteData(index);
        }

        const LatticeVector& GetGlobalSiteCoords() const
        {
          return m_domain->GetGlobalSiteCoords(index);
        }
    };
}

#endif
