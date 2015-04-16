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

    namespace
    {
      int minimumPosImpl(Dimensionless x, size_t range)
      {
        return static_cast<int>(std::floor(x - 0.5 * Dimensionless(range)) + 1);
      }
      int maximumPosImpl(Dimensionless x, size_t range)
      {
        return static_cast<int>(std::floor(x + 0.5 * Dimensionless(range)));
      }
    }

    LatticeVector InterpolationIterator::minimumPosition(LatticePosition const &node, size_t range)
    {
      return LatticeVector(minimumPosImpl(node.x, range),
                           minimumPosImpl(node.y, range),
                           minimumPosImpl(node.z, range));
    }
    LatticeVector InterpolationIterator::maximumPosition(LatticePosition const &node, size_t range)
    {
      return LatticeVector(maximumPosImpl(node.x, range),
                           maximumPosImpl(node.y, range),
                           maximumPosImpl(node.z, range));
    }

    InterpolationIterator interpolationIterator(LatticePosition const &where,
                                                stencil::types stencil)
    {
        // faff because intel
        const unsigned int a = static_cast<unsigned int>(stencil::types::FOUR_POINT);
        const unsigned int b = static_cast<unsigned int>(stencil::types::COSINE_APPROX);
        const unsigned int c = static_cast<unsigned int>(stencil::types::THREE_POINT);
        const unsigned int d = static_cast<unsigned int>(stencil::types::TWO_POINT);
        switch (static_cast<unsigned int>(stencil))
        {
            case a:
                return interpolationIterator<stencil::FourPoint>(where);
            case b:
                return interpolationIterator<stencil::CosineApprox>(where);
            case c:
                return interpolationIterator<stencil::ThreePoint>(where);
            case d:
                return interpolationIterator<stencil::TwoPoint>(where);
        }
        return interpolationIterator<stencil::TwoPoint>(where);
    }
  }
}
