//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_NODE2NODE_H
#define HEMELB_REDBLOOD_NODE2NODE_H

#include "units.h"

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
                                        LatticeDistance const cutoffDistance = 1.0,
                                        size_t exponent = 2)
    {
      if (distance >= cutoffDistance)
      {
        return 0e0;
      }

      LatticeDistance const deltaX = 1;
      return -intensity
          * (std::pow(deltaX / distance, exponent) - std::pow(deltaX / cutoffDistance, exponent));
    }

    // Repulsive force between two nodes
    inline LatticeForceVector node2NodeForce(LatticePosition distance,
                                             LatticeModulus const intensity,
                                             LatticeDistance const cutoffDistance = 1.0,
                                             size_t exponent = 2)
    {
      LatticeDistance const d = distance.GetMagnitude();
      return distance * (node2NodeForce(d, intensity, cutoffDistance, exponent) / d);
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
