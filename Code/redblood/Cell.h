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
#include "redblood/Node2Node.h"
#include "units.h"

namespace hemelb { namespace redblood {

class CellBase : public Mesh {

  public:
  //! \brief Initializes mesh from mesh data
  //! \details This version makes it possible to share the unmodified mesh
  //! across particles.
  //! \param [in] _mesh: Modifyiable mesh
  //! \param [in] _template: Original mesh
  CellBase(Mesh const & _mesh, Mesh const &_template)
      : Mesh(_mesh), template_(_template.GetData()) {}

  //! \brief Initializes mesh from mesh data
  //! \details This version makes it possible to share the unmodified mesh
  //! across particles.
  //! \param [in] _mesh: Modifyiable mesh
  CellBase(Mesh const & _mesh)
      : Mesh(_mesh), template_(new MeshData(*_mesh.GetData())) {}

  //! \brief Initializes mesh from mesh data
  //! On top of taking ownership of the mesh, a unmodifiable copy of the mesh
  //! is also created.
  //! \param [in] _mesh: Modifyiable mesh
  //! \param [in] _template: Original mesh
  CellBase(std::shared_ptr<MeshData> const & _mesh)
       : CellBase(Mesh(_mesh)) {}

  //! Because it is good practice
  virtual ~CellBase() {}

  //! Unmodified mesh
  std::shared_ptr<const MeshData> GetTemplateMesh() const
    { return template_; }
  //! Current mesh
  std::shared_ptr<MeshData> GetCurrentMesh() const
    { return mesh_; }

  //! Facet bending energy
  virtual PhysicalEnergy operator()() const = 0;
  //! Facet bending energy
  virtual PhysicalEnergy operator()(
      std::vector<LatticeForceVector> &_in) const = 0;

  protected:
   //! Unmodified original mesh
   std::shared_ptr<MeshData const> template_;
};

//! Deformable cell for which energy and forces can be computed
class Cell : public CellBase {
public:
  //! Holds all physical parameters
  struct Moduli {
    //! Bending energy parameter
    PhysicalPressure bending;
    //! Surface energy parameter
    PhysicalPressure surface;
    //! Surface volume parameter
    PhysicalPressure volume;
    //! Skalak dilation modulus
    PhysicalPressure dilation;
    //! Skalak strain modulus
    PhysicalPressure strain;
    Moduli() : bending(0), surface(0), volume(0), dilation(0), strain(0) {};
  } moduli;
  //! Node-wall interaction
  Node2NodeForce nodeWall;

  // inheriting constructors
  using CellBase::CellBase;

  //! Facet bending energy
  virtual PhysicalEnergy operator()() const override;
  //! Facet bending energy
  virtual PhysicalEnergy operator()(
      std::vector<LatticeForceVector> &_in) const override;

private:
  // Computes facet bending energy over all facets
  PhysicalEnergy facetBending_() const;
  // Computes facet bending energy over all facets
  PhysicalEnergy facetBending_(std::vector<LatticeForceVector> &_forces) const;
};

//! Typical cell container type
typedef std::vector<std::shared_ptr<CellBase>> CellContainer;

}} // namespace hemelb::redblood
#endif
