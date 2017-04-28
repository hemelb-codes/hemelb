//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/RBCInserter.h"

namespace hemelb
{
  namespace redblood
  {
    CellContainer::value_type RBCInserterWithPerturbation::drop()
    {
      auto const result = RBCInserter::drop();

      // apply rotation
      auto const theta = dtheta * (2e0 * uniformDistribution(randomGenerator) - 1e0);
      auto const phi = dphi * (2e0 * uniformDistribution(randomGenerator) - 1e0);
      using std::cos;
      using std::sin;
      LatticePosition const z(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
      util::Matrix3D const rotation = initialRotation * rotationMatrix(LatticePosition(0, 0, 1), z)
          * initialRotation.transpose();
      *result *= rotation;

      // apply translation
      *result += dx * uniformDistribution(randomGenerator) + dy * uniformDistribution(randomGenerator);
      return result;
    }
  }
} // hemelb::redblood
