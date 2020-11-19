// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
