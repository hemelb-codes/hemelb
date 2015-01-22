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

class CellBase {

  public:
  //! \brief Initializes mesh from mesh data
  //! \param [in] _vertices: deformable vertices that define the cell. These
  //!    values are *not* modified by the scale.
  //! \param [in] _template: Original mesh. A shallow copy is made of this
  //!    object.
  //! \param [in] _scale: scales template by a given amount
  //!    The scale is added during internal operations. The template will still
  //!    refer to the same data in memory.
  CellBase(
      MeshData::t_Vertices && _vertices,
      Mesh const &_template,
      Dimensionless _scale = 1e0
  ) : vertices_(std::move(_vertices)), template_(_template), scale_(_scale) {
    assert(scale_ > 1e-12);
  }
  //! \brief Initializes mesh from mesh data
  //! \param [in] _vertices: deformable vertices that define the cell. These
  //!    values are *not* modified by the scale.
  //! \param [in] _template: Original mesh. A shallow copy is made of this
  //!    object.
  //! \param [in] _scale: scales template by a given amount
  //!    The scale is added during internal operations. The template will still
  //!    refer to the same data in memory.
  CellBase(
      MeshData::t_Vertices const & _vertices,
      Mesh const &_template,
      Dimensionless _scale = 1e0
  ) : vertices_(_vertices), template_(_template), scale_(_scale) {
    assert(scale_ > 1e-12);
  }

  //! \brief Initializes mesh from mesh data
  //! \param [in] _mesh: deformable vertices that define the cell are copied
  //!    from this mesh. These values are *not* modified by the scale.
  //! \param [in] _template: Original mesh. A shallow copy is made of this
  //!    object.
  //! \param [in] _scale: scales template by a given amount
  //!    The scale is added during internal operations. The template will still
  //!    refer to the same data in memory.
  CellBase(Mesh const & _mesh, Mesh const &_template, Dimensionless _scale=1e0)
      : CellBase(_mesh.GetVertices(), _template, _scale) {}

  //! \brief Initializes mesh from mesh data
  //! \param [in] _mesh: Modifyiable mesh and template. Deep copies are made of
  //!   both.
  CellBase(Mesh const & _mesh)
      : CellBase(_mesh.GetVertices(), _mesh.clone()) {}

  //! \brief Initializes mesh from mesh data
  //! \param [in] _mesh: Modifyiable mesh and template. Deep copies are made of
  //!   both
  CellBase(std::shared_ptr<MeshData> const & _mesh)
       : CellBase(_mesh->vertices, Mesh(*_mesh)) {}

  void operator=(Mesh const& _mesh) {
    template_ = _mesh;
    vertices_ = template_.GetVertices();
    scale_ = 1e0;
  }

  //! Because it is good practice
  virtual ~CellBase() {}

  //! Unmodified mesh
  Mesh const & GetTemplateMesh() const { return template_; }
  //! Facets for the mesh
  MeshData::t_Facets const & GetFacets() const {
    return template_.GetData()->facets;
  }
  //! Vertices of the cell
  MeshData::t_Vertices const & GetVertices() const { return vertices_; }
  //! Vertices of the cell
  MeshData::t_Vertices& GetVertices() { return vertices_; }
  //! Topology of the (template) mesh
  std::shared_ptr<MeshTopology const> GetTopology() const {
    return template_.GetTopology();
  }
  size_t GetNumberOfNodes() const { return vertices_.size(); }

  //! Facet bending energy
  virtual PhysicalEnergy operator()() const = 0;
  //! Facet bending energy
  virtual PhysicalEnergy operator()(
      std::vector<LatticeForceVector> &_in) const = 0;

  //! Scale mesh around barycenter
  void operator*=(Dimensionless const &_scale);
  //! Scale by matrix around barycenter
  void operator*=(util::Matrix3D const &_scale);
  //! Translate mesh
  void operator+=(LatticePosition const &_offset);
  //! Transform mesh
  void operator+=(std::vector<LatticePosition> const &_displacements);

  MeshData::t_Vertices::value_type GetBarycenter() const {
    return barycenter(vertices_);
  }

  //! Scale to apply to the template mesh
  void SetScale(Dimensionless _scale) {
    assert(_scale > 1e-12);
    scale_ = _scale;
  }
  //! Scale to apply to the template mesh
  Dimensionless GetScale() const { return scale_; }

  protected:
   //! Holds list of vertices for this cell
   MeshData::t_Vertices vertices_;
   //! Unmodified original mesh
   Mesh template_;
   //! Scale factor for the template;
   Dimensionless scale_;
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
typedef std::set<std::shared_ptr<CellBase>> CellContainer;

}} // namespace hemelb::redblood
#endif
