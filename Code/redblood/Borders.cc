//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/Borders.h"

namespace hemelb
{
  namespace redblood
  {
    BorderBoxIterator::BorderBoxIterator(size_t const nearness)
      : current(-1, -1, -1), isValid(true),
        isTop(nearness bitand static_cast<size_t>(Borders::TOP)),
        isBottom(nearness bitand static_cast<size_t>(Borders::BOTTOM)),
        isNorth(nearness bitand static_cast<size_t>(Borders::NORTH)),
        isSouth(nearness bitand static_cast<size_t>(Borders::SOUTH)),
        isWest(nearness bitand static_cast<size_t>(Borders::WEST)),
        isEast(nearness bitand static_cast<size_t>(Borders::EAST))
    {
      current.x = isBottom ? -1: 0;
      current.y = isSouth ? -1: 0;
      current.z = isWest ? -1: 0;
      isValid = isTop or isBottom or isSouth or isNorth or isWest or isEast;
      if(isValid and current.x == 0 and current.y == 0 and current.z == 0)
      {
        try
        {
          this->operator++();
        }
        catch(...)
        {
          isValid = false;
          throw;
        }
      }
    }

    void BorderBoxIterator::operator++()
    {
      if(not isValid)
      {
        return;
      }
      bool const hasCycled = true;
      auto increment = [this, hasCycled](value_type::value_type &x, value_type::value_type lower)
      {
        ++x;
        if(x == 2)
        {
          x = lower;
          return hasCycled;
        }
        return not hasCycled;
      };
      auto wannaSee = [this](value_type const &current) -> bool
      {
        return ((current.x == -1 and isBottom) or (current.x == 1 and isTop) or current.x == 0)
          and ((current.y == -1 and isSouth) or (current.y == 1 and isNorth) or current.y == 0)
          and ((current.z == -1 and isWest) or (current.z == 1 and isEast) or current.z == 0)
          and (current.x != 0 or current.y != 0 or current.z != 0);
      };

      value_type const lower(isBottom ? -1: 0, isSouth ? -1: 0, isWest ? -1: 0);
      do
      {
        bool const cycled0 = increment(current[0], lower[0]) == hasCycled;
        bool const cycled1 = cycled0 ? increment(current[1], lower[1]) == hasCycled: false;
        isValid = cycled1 ? increment(current[2], lower[2]) != hasCycled: true;
      }
      while(isValid and not wannaSee(current));
    }
  }
} // namespace hemelb::redblood
