//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_PARTICLE_H
#define HEMELB_REDBLOOD_PARTICLE_H

#include "redblood/Mesh.h"
#include "units.h"

namespace hemelb { namespace redblood {


//! Deformable particle for which energy and forces can be computed
class Particle : public Mesh {
public:
  //! Holds all physical parameters
  struct Moduli {
    // Bending energy parameter
    PhysicalPressure bending;
    // Surface energy parameter
    PhysicalPressure surface;
    // Surface volume parameter
    PhysicalPressure volume;
    // Skalak dilation modulus
    PhysicalPressure dilation;
    // Skalak strain modulus
    PhysicalPressure strain;
    Moduli() : bending(0), surface(0), volume(0), dilation(0), strain(0) {};
  } moduli;

  //! \brief Initializes mesh from mesh data
  //! \details This version makes it possible to share the unmodified mesh
  //! across particles.
  //! \param [in] _mesh: Modifyiable mesh
  //! \param [in] _template: Original mesh
  Particle(Mesh const & _mesh, Mesh const &_template)
      : Mesh(_mesh), template_(_template.GetData()) {}

  //! \brief Initializes mesh from mesh data
  //! On top of taking ownership of the mesh, a unmodifiable copy of the mesh
  //! is also created.
  //! \param [in] _mesh: Modifyiable mesh
  //! \param [in] _template: Original mesh
  Particle   (boost::shared_ptr<MeshData> const & _mesh)
       : Mesh(_mesh), template_(new MeshData(*_mesh)) {}

  //! Unmodified mesh
  boost::shared_ptr<const MeshData> GetTemplateMesh() const
    { return template_; }
  //! Current mesh
  boost::shared_ptr<MeshData> GetCurrentMesh() const
    { return mesh_; }

  //! Facet bending energy
  PhysicalEnergy operator()() const;
  //! Facet bending energy
  PhysicalEnergy operator()(std::vector<LatticeForceVector> &_in) const;


protected:
 //! Unmodified original mesh
 boost::shared_ptr<MeshData const> template_;

private:
  // Computes facet bending energy over all facets
  PhysicalEnergy facetBending_() const;
  // Computes facet bending energy over all facets
  PhysicalEnergy facetBending_(std::vector<LatticeForceVector> &_forces) const;
};

}} // namespace hemelb::redblood
#endif
