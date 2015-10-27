//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_SITE_ITERATOR_H
#define HEMELB_UNITTESTS_SITE_ITERATOR_H

#include <memory>
#include <iterator>
#include "units.h"
#include "geometry/LatticeData.h"
#include "geometry/Site.h"

namespace hemelb
{
  namespace unittests
  {
    //! Enables ranged-for iterations over lattice sites
    class site_iterator
    {
      public:
        typedef hemelb::geometry::Site<geometry::LatticeData const> value_type;
        typedef value_type const& reference;
        typedef std::unique_ptr<value_type> pointer;
        typedef std::ptrdiff_t difference_type;
        typedef std::input_iterator_tag iterator_category;

        site_iterator(geometry::LatticeData const& source, site_t current)
          : current(current), source(source)
        {
        }

        value_type operator*() const
        {
          return source.GetSite(current);
        }
        site_iterator & operator++()
        {
          ++current;
          return *this;
        }
        site_iterator operator++(int)
        {
          auto const result = site_iterator(source, current);
          ++current;
          return result;
        }
        bool operator!=(site_iterator const &other) const
        {
          assert(&source == &other.source);
          return current != other.current;
        }

      protected:
        //! Current site
        site_t current;
        //! Reference to data source
        geometry::LatticeData const & source;
    };
  }
}

namespace std
{
  inline hemelb::unittests::site_iterator begin(hemelb::geometry::LatticeData const& latDat)
  {
    return {latDat, hemelb::site_t(0)};
  }
  inline hemelb::unittests::site_iterator end(hemelb::geometry::LatticeData const& latDat)
  {
    return {latDat, latDat.GetTotalFluidSites()};
  }
}

#endif /* HEMELB_UNITTESTS_FOURCUBELATTICEDATA_H */
