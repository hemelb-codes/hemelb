//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
namespace details { namespace {

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

//! Iterates over vertices of a mesh and the nearby nodes of a grid
//! The functor argument is called with the current vertex index, the
//! global site index triplet, and the associated interpolation weight.
template<class T_FUNCTOR> void spreadForce2Grid(
    Particle const &_particle,
    T_FUNCTOR _functor,
    stencil::types _stencil
) {
  typedef MeshData::t_Vertices::const_iterator const_iterator;
  // Spread them onto lattice
  const_iterator i_vertex = _particle.GetVertices().begin();
  const_iterator const i_end = _particle.GetVertices().end();
  for(size_t i(0); i_vertex != i_end; ++i_vertex, ++i) {
    InterpolationIterator spreader
      = interpolationIterator(*i_vertex, _stencil);
    for(; spreader; ++spreader)
      _functor(i, *spreader, spreader.weight());
  }
}

class SpreadForces {
  public:
    SpreadForces(
        std::vector<LatticePosition> const &_forces,
        geometry::LatticeData &_latticeData
    ) : latticeData_(_latticeData), forces_(_forces) {}

    void operator()(
        size_t _vertex, LatticeVector const &_site,
        Dimensionless _weight) {
      proc_t procid;
      site_t siteid;
      if(latticeData_.GetContiguousSiteId(_site, procid, siteid))
        geometry::Site<geometry::LatticeData>(latticeData_.GetSite(_site))
          .AddToForce(forces_[_vertex] * _weight);
    }

  protected:
    geometry::LatticeData &latticeData_;
    std::vector<LatticeForceVector> const & forces_;
};

template<class LATTICE>
class SpreadForcesAndWallForces : public SpreadForces {
  public:
    SpreadForcesAndWallForces(
          Particle const &_particle,
          std::vector<LatticePosition> const &_forces,
          geometry::LatticeData &_latticeData
    ) : SpreadForces(_forces, _latticeData), particle_(_particle) {}
    void operator()(
          size_t _vertex, LatticeVector const &_site,
          Dimensionless _weight) {
      proc_t procid;
      site_t siteid;
      if(not latticeData_.GetContiguousSiteId(_site, procid, siteid))
        return;
      geometry::Site<geometry::LatticeData> site(latticeData_.GetSite(siteid));
      site.AddToForce(forces_[_vertex] * _weight);
      LatticePosition const vertex(particle_.GetVertex(_vertex));
      // std::cout << "v " << vertex << "\n";
      for(size_t i(1); i < LATTICE::NUMVECTORS; ++i) {
        PhysicalDistance const distance = site.GetWallDistance<LATTICE>(i);
        if(not site.HasWall(i)) continue;
        // Direction of streaming from wall to this site
        LatticePosition const direction = LatticePosition(
            LATTICE::CX[i], LATTICE::CY[i], LATTICE::CZ[i]
        );
        LatticePosition const wallnode = LatticePosition(_site)
          + direction.GetNormalised() * distance;
        site.AddToForce(particle_.nodeWall(vertex, wallnode) * _weight);
      }
    }
  protected:
    Particle const &particle_;
};

}} // namespace details::anonymous

