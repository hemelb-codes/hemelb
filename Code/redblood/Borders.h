// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_BORDERS_H
#define HEMELB_REDBLOOD_BORDERS_H

#include <cassert>
#include "redblood/DivideConquer.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    //! Names for each border
    enum class Borders
      : size_t
      {
        NONE = 0,
      TOP = 1,
      BOTTOM = 2,
      NORTH = 4,
      SOUTH = 8,
      WEST = 16,
      EAST = 32,
      LAST = 64,
      CENTER = 128
    };
    //! Converts border enum to actual lattice direction
    template<class T> util::Vector3D<T> direction(Borders border);
    //! Converts border enum to actual lattice direction
    template<class T> util::Vector3D<T> direction(size_t border);
    //! \brief Encodes how close a vertex is to neighboring boxes
    //! \param[in] dnc Divide and conquer boxes
    //! \param[in] vertex Position to investigate
    //! \param[in] haloLenght the meaning of close
    template<class T>
    size_t figureNearness(DivideConquer<T> const &dnc, LatticePosition const &vertex,
                          LatticeDistance const &haloLength);

    //! \brief Iterates only over nearest borders
    class BorderBoxIterator
    {
      public:
        //! Type returned by dereference
        typedef LatticeVector value_type;
        typedef value_type const& const_reference;
        typedef value_type const* const_pointer;

        //! \brief Iterates over all boxes for which a bit is set
        //! \see Borders
        BorderBoxIterator(size_t const nearness);
        //! Current direction to look at
        const_reference operator*() const
        {
          return current;
        }
        //! Current direction to look at
        const_pointer operator->() const
        {
          return &current;
        }
        //! Whether this iterator is valid
        operator bool() const
        {
          return isValid;
        }
        //! Go to next position
        BorderBoxIterator& operator++();

      protected:
        //! Current direction;
        value_type current;
        //! Whether the object is valid
        bool isValid;
        //! Which box to do
        bool doCenter, doTop, doBottom, doNorth, doSouth, doWest, doEast;

        //! Returns true if this is a box we want to visit
        bool wannaSee(value_type const &current) const;
    };

    template<class T>
    util::Vector3D<T> direction(Borders border)
    {
      switch (border)
      {
        case Borders::TOP:
          return {1, 0, 0};
        case Borders::BOTTOM:
          return {-1, 0, 0};
        case Borders::NORTH:
          return {0, 1, 0};
        case Borders::SOUTH:
          return {0, -1, 0};
        case Borders::WEST:
          return {0, 0, -1};
        case Borders::EAST:
          return {0, 0, 1};
        default:
          return {0, 0, 0};
      };
      return {0, 0, 0};
    }
    template<class T> util::Vector3D<T> direction(size_t border)
    {
      Borders const b(static_cast<Borders>(border));
      assert(b == Borders::TOP or Borders(b) == Borders::BOTTOM or b == Borders::NORTH
          or b == Borders::SOUTH or b == Borders::EAST or b == Borders::WEST);
      return direction<T>(b);
    }

    template<class T>
    size_t figureNearness(DivideConquer<T> const &dnc, LatticeVector const &key,
                          LatticePosition const &vertex, LatticeDistance const &haloLength)
    {
      if (haloLength > dnc.GetBoxSize())
      {
        return static_cast<size_t>(Borders::CENTER) bitor static_cast<size_t>(Borders::NORTH)
            bitor static_cast<size_t>(Borders::SOUTH) bitor static_cast<size_t>(Borders::EAST)
            bitor static_cast<size_t>(Borders::WEST) bitor static_cast<size_t>(Borders::TOP)
            bitor static_cast<size_t>(Borders::BOTTOM);
      }
      else if (haloLength < 1e-12)
      {
        return static_cast<size_t>(Borders::CENTER);
      }

      size_t result = 0;
      typedef LatticePosition::value_type Distance;

      for (size_t d(1); d < (1 << 6); d <<= 1)
      {
        Borders const border(static_cast<Borders>(d));
        LatticePosition const translated(direction<Distance>(border) * haloLength);

        if (not (key == dnc.DowngradeKey(vertex + translated)))
        {
          result |= d;
        }
      }

      return result bitor static_cast<size_t>(Borders::CENTER);
    }
    template<class T>
    size_t figureNearness(DivideConquer<T> const &dnc, LatticePosition const &vertex,
                          LatticeDistance const &haloLength)
    {
      return figureNearness<T>(dnc, dnc.DowngradeKey(vertex), vertex, haloLength);
    }

  }
} // hemelb::redblood

#endif
