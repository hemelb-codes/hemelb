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


namespace hemelb { namespace redblood {

  //! \brief Force meant to prevent nodes from collapsing one onto the other
  //! \param[in] _distance: Square of the distance between the two nodes
  //! \param[in] _intensity: K_int, strength of the interaction
  //! \param[in] _cutoffDistance: Maximum interaction distance, squared
  //! \param[in] _exponent
  inline PhysicalForce node2NodeForce(
      PhysicalDistance const _distance,
      PhysicalForce const _intensity,
      PhysicalDistance const _cutoffDistance = 1.0,
      size_t _exponent = 2
  ) {
    if(_distance >= _cutoffDistance) return 0e0;
    PhysicalDistance const deltaX = 1;
    return -_intensity * (
        std::pow(deltaX / _distance, _exponent)
        - std::pow(deltaX / _cutoffDistance, _exponent)
    );
  }

  // Repulsive force between two nodes
  inline LatticeForceVector node2NodeForce(
      LatticePosition _distance,
      PhysicalForce const _intensity,
      PhysicalDistance const _cutoffDistance = 1.0,
      size_t _exponent = 2
  ) {
    PhysicalDistance const d = _distance.GetMagnitude();
    return _distance * (
        node2NodeForce(d, _intensity, _cutoffDistance, _exponent) / d
    );
  }

  // Repulsive force felt by A from interaction with B
  inline LatticeForceVector node2NodeForce(
      LatticePosition _A, LatticePosition _B,
      PhysicalForce const _intensity,
      PhysicalDistance const _cutoffDistance = 1.0,
      size_t _exponent = 2
  ) {
    return node2NodeForce(_B - _A, _intensity, _cutoffDistance, _exponent);
  }


  //! Holds node-node interaction parameters
  struct Node2NodeForce {
    //! Strength of the interaction
    PhysicalForce intensity;
    //! Maximum distance of the interaction
    PhysicalDistance cutoff;
    //! Power exponent
    size_t exponent;

    Node2NodeForce(
        PhysicalForce _intensity = 0.0,
        PhysicalDistance _cutoff = 1.0,
        size_t _exponent = 2
    ) : intensity(_intensity), cutoff(_cutoff), exponent(_exponent) {}


    PhysicalForce operator()(PhysicalDistance const &_distance) const {
      return node2NodeForce(_distance, intensity, cutoff, exponent);
    }
    LatticeForceVector operator()(LatticePosition const &_distance) const {
      return node2NodeForce(_distance, intensity, cutoff, exponent);
    }
    LatticeForceVector operator()(
        LatticePosition const &_A, LatticePosition const &_B) const {
      return node2NodeForce(_A, _B, intensity, cutoff, exponent);
    }
  };
}} // hemelb::redblood

#endif
