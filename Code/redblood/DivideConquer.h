//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_DIVIDE_AND_CONQUER_H
#define HEMELB_REDBLOOD_DIVIDE_AND_CONQUER_H

#include <vector>
#include <map>
#include <cmath>
#include "units.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace redblood
  {
    namespace details
    {
      namespace
      {
        // Short-hand to get type of the base of Divide and Conquer class
        // Also aggregates key type and compare functor
        template <class T>
        struct DnCBase
        {
          typedef LatticeVector key_type;
          struct CompareKeys
          {
            bool operator()(key_type const &a, key_type const &b) const
            {
              if (a.x > b.x)
              {
                return false;
              }
              else if (a.x < b.x)
              {
                return true;
              }

              if (a.y > b.y)
              {
                return false;
              }
              else if (a.y < b.y)
              {
                return true;
              }

              return a.z < b.z;
            }
          };
          typedef std::multimap<key_type, T, CompareKeys> type;
        };
      }
    }

    //! \brief Multimap for divide and conquer algorithms
    //! \details Items at a position x are mapped into boxes of a given size.
    template <class T>
    class DivideConquer : public details::DnCBase<T>::type
    {
      typedef typename details::DnCBase<T>::type base_type;

      public:
      typedef typename base_type::key_type key_type;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::reference reference;
      typedef typename base_type::const_reference const_reference;
      typedef typename base_type::iterator iterator;
      typedef typename base_type::const_iterator const_iterator;
      typedef std::pair<iterator, iterator> range;
      typedef std::pair<const_iterator, const_iterator> const_range;

      //! Constructor sets size of cutoff
      DivideConquer(PhysicalDistance boxsize) : base_type(), boxsize(boxsize)
      {
      }
      //! Insert into divide and conquer container
      iterator insert(LatticePosition const &pos, T const &value)
      {
        return base_type::insert(value_type(DowngradeKey(pos), value));
      }
      //! Insert into divide and conquer container
      iterator insert(key_type const &pos, T const &value)
      {
        return base_type::insert(value_type(pos, value));
      }
      //! All objects in a single divide and conquer box
      range equal_range(LatticePosition const &pos)
      {
        return DivideConquer<T>::equal_range(DowngradeKey(pos));
      }
      //! All objects in a single divide and conquer box
      range equal_range(key_type const &pos)
      {
        return base_type::equal_range(pos);
      }
      //! All objects in a single divide and conquer box
      const_range equal_range(LatticePosition const &pos) const
      {
        return DivideConquer<T>::equal_range(DowngradeKey(pos));
      }
      //! All objects in a single divide and conquer box
      const_range equal_range(key_type const &pos) const
      {
        return base_type::equal_range(pos);
      }

      //! Length of each box
      PhysicalDistance GetBoxSize() const
      {
        return boxsize;
      }

      //! Converts from position to box index
      key_type DowngradeKey(LatticePosition const &pos) const
      {
        return key_type(static_cast<LatticeCoordinate>(std::floor(pos.x / boxsize)),
                        static_cast<LatticeCoordinate>(std::floor(pos.y / boxsize)),
                        static_cast<LatticeCoordinate>(std::floor(pos.z / boxsize)));
      }
      //! No conversion since in box index type already
      key_type DowngradeKey(key_type const &pos) const
      {
        return pos;
      }

      protected:
      PhysicalDistance const boxsize;
    };
  }
}  // hemelb::redblood

#endif
