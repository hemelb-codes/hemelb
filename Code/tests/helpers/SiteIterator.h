// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_SITEITERATOR_H
#define HEMELB_TESTS_HELPERS_SITEITERATOR_H

#include <memory>
#include <iterator>
#include "units.h"
#include "geometry/Domain.h"
#include "geometry/Site.h"

namespace hemelb
{
  namespace tests
  {
    //! Enables ranged-for iterations over lattice sites
    class site_iterator
    {
      public:
        using value_type = hemelb::geometry::Site<geometry::FieldData const>;
        using reference = value_type const& ;
        using pointer = std::unique_ptr<value_type>;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        site_iterator(geometry::FieldData const& source, site_t current)
          : m_current(current), m_source(&source)
        {
        }

        value_type operator*() const
        {
          return m_source->GetSite(m_current);
        }
        site_iterator & operator++()
        {
          ++m_current;
          return *this;
        }
        site_iterator operator++(int)
        {
          auto const result = site_iterator(*m_source, m_current);
          ++m_current;
          return result;
        }
        friend bool operator!=(site_iterator const &lhs, site_iterator const &rhs)
        {
          assert(lhs.m_source == rhs.m_source);
          return lhs.m_current != rhs.m_current;
        }

      protected:
        //! Current site
        site_t m_current;
        //! Reference to data source
        geometry::FieldData const* m_source;
    };
  }
}

// TODO: remove UB!
namespace std
{
  inline hemelb::tests::site_iterator begin(hemelb::geometry::FieldData const& latDat)
  {
    return {latDat, hemelb::site_t(0)};
  }
  inline hemelb::tests::site_iterator end(hemelb::geometry::FieldData const& latDat)
  {
    return {latDat, latDat.GetDomain().GetTotalFluidSites()};
  }
}

#endif
