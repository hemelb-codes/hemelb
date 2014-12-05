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
  //! \param[in] _distanceSquared: Square of the distance between the two nodes
  //! \param[in] _intensity: K_int, strength of the interaction
  //! \param[in] _cutoffDistanceSquared: Maximum interaction distance, squared
  //! \param[in] _exponent
  inline PhysicalForce node2NodeForce(
      PhysicalDistance const _distanceSquared,
      PhysicalForce const _intensity,
      PhysicalDistance const _cutoffDistanceSquared = 1.0,
      size_t _exponent = 1
  ) {
    if(_distanceSquared >= _cutoffDistanceSquared) return 0e0;
    PhysicalDistance const deltaXSquared = 1;
    return -_intensity * (
        std::pow(deltaXSquared / _distanceSquared, _exponent)
        - std::pow(deltaXSquared / _cutoffDistanceSquared, _exponent)
    );
  }

  // Repulsive force between two nodes
  inline LatticeForceVector node2NodeForce(
      LatticePosition _distance,
      PhysicalForce const _intensity,
      PhysicalDistance const _cutoffDistanceSquared = 1.0,
      size_t _exponent = 1
  ) {
    PhysicalDistance const d = _distance.GetMagnitudeSquared();
    return _distance * (
        node2NodeForce(d, _intensity, _cutoffDistanceSquared, _exponent)
        / std::sqrt(d)
    );
  }

  // Repulsive force felt by A from interaction with B
  inline LatticeForceVector node2NodeForce(
      LatticePosition _A, LatticePosition _B,
      PhysicalForce const _intensity,
      PhysicalDistance const _cutoffDistanceSquared = 1.0,
      size_t _exponent = 1
  ) {
    return node2NodeForce(
        _B - _A, _intensity, _cutoffDistanceSquared, _exponent);
  }

}} // hemelb::redblood

#endif
