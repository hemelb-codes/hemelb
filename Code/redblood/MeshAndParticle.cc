//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/MeshAndParticle.h"
#include "redblood/VelocityInterpolation.h"

namespace hemelb { namespace redblood {

  namespace {
    // Loops over nodes and does something
    template<class T_FUNCTOR>
      void nodeLoop(
        Particle const &_particle, geometry::LatticeData const &_latdat,
        T_FUNCTOR action,
        stencil::types _stencil = stencil::FOUR_POINT
      );

    struct UpdateVector {
      UpdateVector(
          stencil::types _stencil,
          Particle const &_particle,
          geometry::LatticeData const &_latdat,
          std::vector<LatticePosition> &_vector
      ) : stencil(_stencil), latDat(_latdat) {
        _vector.resize(_particle.GetNumberOfNodes());
        i_current = _vector.begin();
      }
      void operator()(MeshData::t_Vertices::const_reference _vertex) {
        *i_current++ = interpolateVelocity(latDat, _vertex, stencil);
      }
      stencil::types const stencil;
      geometry::LatticeData const &latDat;
      std::vector<LatticePosition> :: iterator i_current;
    };
  }

  bool compute_displacement(
      Particle const &_particle, geometry::LatticeData const &_latdat,
      std::vector<LatticePosition> &_displacement,
      stencil::types _stencil = stencil::FOUR_POINT
  ) {
  }

}} // hemelb::redblood
