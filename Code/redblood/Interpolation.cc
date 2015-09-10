//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/Interpolation.h"
#include <cmath>

namespace hemelb
{
  namespace redblood
  {
    void IndexIterator::operator++()
    {
#ifndef NDEBUG

      if (not IsValid())
      {
        throw Exception() << "Cannot increment invalid iterator\n";
      }

#endif

      if (current[2] < max[2])
      {
        ++current[2];
        return;
      }

      current[2] = min[2];

      if (current[1] < max[1])
      {
        ++current[1];
        return;
      }

      current[1] = min[1];
      ++current[0];
    }
  }
}
