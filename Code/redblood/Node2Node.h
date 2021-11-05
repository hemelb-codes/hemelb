// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_NODE2NODE_H
#define HEMELB_REDBLOOD_NODE2NODE_H

#include <cassert>
#include "Exception.h"
#include "units.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace redblood
  {
    //! \brief Force meant to prevent nodes from collapsing one onto the other
    //! \param[in] distance: Square of the distance between the two nodes
    //! \param[in] intensity: K_int, strength of the interaction
    //! \param[in] cutoffDistance: Maximum interaction distance, squared
    //! \param[in] exponent
    inline LatticeForce node2NodeForce(LatticeDistance const distance,
                                       LatticeModulus const intensity,
                                       LatticeDistance const cutoffDistance = 1.0, size_t exponent =
                                           2)
    {
      if (distance >= cutoffDistance)
      {
        return 0e0;
      }

      LatticeDistance const deltaX = 1;
      auto const d_to_pow = std::pow(deltaX / distance, exponent);
      auto const d_to_pow0 = std::pow(deltaX / cutoffDistance, exponent);
      if (std::isnan(d_to_pow) or std::isinf(d_to_pow))
      {
        std::string const message = "*** node-node interaction is NaN or infinite";
        log::Logger::Log<log::Critical, log::Singleton>(message);
        throw Exception() << message;
      }
      return -intensity * (d_to_pow - d_to_pow0);
    }

    // Repulsive force between two nodes
    inline LatticeForceVector node2NodeForce(LatticePosition distance,
                                             LatticeModulus const intensity,
                                             LatticeDistance const cutoffDistance = 1.0,
                                             size_t exponent = 2)
    {
      LatticeDistance const d = distance.GetMagnitude();
      assert(d > 1e-12); // can't determine direction of the force
      auto const magnitude = node2NodeForce(d, intensity, cutoffDistance, exponent);
      return distance * (magnitude / d);
    }

    // Repulsive force felt by A from interaction with B
    inline LatticeForceVector node2NodeForce(LatticePosition A, LatticePosition B,
                                             LatticeModulus const intensity,
                                             LatticeDistance const cutoffDistance = 1.0,
                                             size_t exponent = 2)
    {
      return node2NodeForce(B - A, intensity, cutoffDistance, exponent);
    }

    //! Holds node-node interaction parameters
    class Node2NodeForce
    {
      public:
        //! Strength of the interaction
        LatticeModulus intensity;
        //! Maximum distance of the interaction
        LatticeDistance cutoff;
        //! Power exponent
        size_t exponent;

        Node2NodeForce(LatticeModulus intensity = 0.0, LatticeDistance cutoff = 1.0,
                       size_t exponent = 2) :
            intensity(intensity), cutoff(cutoff), exponent(exponent)
        {
        }

        LatticeForce operator()(LatticeDistance const &distance) const
        {
          return node2NodeForce(distance, intensity, cutoff, exponent);
        }
        LatticeForceVector operator()(LatticePosition const &distance) const
        {
          return node2NodeForce(distance, intensity, cutoff, exponent);
        }
        LatticeForceVector operator()(LatticePosition const &A, LatticePosition const &B) const
        {
          return node2NodeForce(A, B, intensity, cutoff, exponent);
        }
    };
  }
} // hemelb::redblood

#endif
