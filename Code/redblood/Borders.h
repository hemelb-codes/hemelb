//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_BORDERS_INTERACTION_H
#define HEMELB_REDBLOOD_BORDERS_INTERACTION_H

#include <cassert>
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    //! Names for each border
    enum class Borders : size_t
    {
      NONE = 0,
      TOP = 1,
      BOTTOM = 2,
      NORTH = 4,
      SOUTH = 8,
      WEST = 16,
      EAST = 32,
      LAST = 64
    };
    //! Converts border enum to actual lattice direction
    template<class T> util::Vector3D<T> direction(Borders border);
    //! Converts border enum to actual lattice direction
    template<class T> util::Vector3D<T> direction(size_t border);

    //! \brief Iterates only over nearest borders
    //! \details 
    class LoopThroughBorders
    {
      public:
        LoopThroughBorders(size_t const nearness)
          : current(static_cast<size_t>(Borders::NONE)), nearness(nearness)
        {
        }
        //! Current direction to look at
        LatticePosition operator*() const;
        //! Whether this iterator is valid
        operator bool() const;
        //! Go to next position
        void operator++(int);

      protected:
        //! Current direction;
        size_t current;
        //! Which directions to investigate
        size_t nearness;
    };

    template<class T> util::Vector3D<T> direction(Borders border)
    {
      switch (border)
      {
        case Borders::TOP:    return util::Vector3D<T>(1, 0, 0);
        case Borders::BOTTOM: return util::Vector3D<T>(-1, 0, 0);
        case Borders::NORTH:  return util::Vector3D<T>(0, 1, 0);
        case Borders::SOUTH:  return util::Vector3D<T>(0, -1, 0);
        case Borders::WEST:   return util::Vector3D<T>(0, 0, -1);
        case Borders::EAST:   return util::Vector3D<T>(0, 0, 1);
        default:     return util::Vector3D<T>(0, 0, 0);
      };
      return LatticeVector(0, 0, 0);
    }
    template<class T> util::Vector3D<T> direction(size_t border)
    {
      Borders const b(static_cast<Borders>(border));
      assert(b == Borders::TOP or Borders(b) == Borders::BOTTOM or b == Borders::NORTH or
          b == Borders::SOUTH or b == Borders::EAST or b == Borders::WEST);
      return direction<T>(b);
    }
  }
} // hemelb::redblood

#endif
