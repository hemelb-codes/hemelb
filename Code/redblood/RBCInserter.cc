// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
