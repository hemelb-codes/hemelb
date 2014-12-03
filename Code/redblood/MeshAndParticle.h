//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <vector>
#include "units.h"
#include "redblood/Particle.h"
#include "geometry/LatticeData.h"

namespace hemelb { namespace redblood {

  namespace {
    // Loops over nodes, computes velocity and does something
    // struct + member makes up for lack of partial function
    // specialization in c++ pre 11
    template<class T_KERNEL> struct VelocityNodeLoop {
      VelocityNodeLoop(
          stencil::types _stencil,
          Particle const &_particle,
          geometry::LatticeData const &_latDat
      ) : stencil(_stencil), particle(_particle), latticeData(_latDat) {}
      // Loop and does something
      template<class T_FUNCTOR> void loop(T_FUNCTOR apply) {
        typedef MeshData::t_Vertices::const_iterator const_iterator;
        const_iterator i_current = particle.GetVertices().begin();
        const_iterator const i_end = particle.GetVertices().end();
        for(; i_current != i_end; ++i_current) {
          PhysicalVelocity const velocity
            = interpolateVelocity<T_KERNEL>(latticeData, *i_current, stencil);
          std::cout << "   vel: " << velocity << " from " << *i_current << std::endl;
          apply(velocity);
        }
      }

      stencil::types const stencil;
      Particle const &particle;
      geometry::LatticeData const &latticeData;
    };

    //! Updates an assignable iterator of some kind
    template<class T_ITERATOR> struct TransformIterator {
      T_ITERATOR iterator;
      TransformIterator(T_ITERATOR _iterator) : iterator(_iterator) {}
      void operator()(typename T_ITERATOR::value_type const & _value) {
        *(iterator++) = _value;
      }
    };

    //! Updates an assignable iterator of some kind
    template<class T_ITERATOR>
      TransformIterator<T_ITERATOR> transform_iterator(T_ITERATOR iterator) {
        return TransformIterator<T_ITERATOR>(iterator);
      }
  }

  // Displacement/deformation of the particles inferred from lattice velocities
  template<class T_KERNEL> void compute_displacement(
      Particle const &_particle,
      geometry::LatticeData const &_latDat,
      stencil::types _stencil,
      std::vector<LatticePosition> &_displacements
  ) {
    _displacements.resize(_particle.GetNumberOfNodes());
    VelocityNodeLoop<T_KERNEL>(_stencil, _particle, _latDat)
      .loop(transform_iterator(_displacements.begin()));
  }
}} // hemelb::redblood
