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

#include <set>
#include <utility>
#include "redblood/Mesh.h"
#include "redblood/Node2Node.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    class CellBase
    {
      public:
        //! \brief Initializes mesh from mesh data
        //! \param [in] verticesIn: deformable vertices that define the cell. These
        //!    values are *not* modified by the scale.
        //! \param [in] origMesh: Original mesh. A shallow copy is made of this
        //!    object.
        //! \param [in] scaleIn: scales template by a given amount
        //!    The scale is added during internal operations. The template will still
        //!    refer to the same data in memory.
        CellBase(MeshData::Vertices &&verticesIn, Mesh const &origMesh,
                 Dimensionless scaleIn = 1e0);
        //! \brief Initializes mesh from mesh data
        //! \param [in] verticesIn: deformable vertices that define the cell. These
        //!    values are *not* modified by the scale.
        //! \param [in] origMesh: Original mesh. A shallow copy is made of this
        //!    object.
        //! \param [in] scaleIn: scales template by a given amount
        //!    The scale is added during internal operations. The template will still
        //!    refer to the same data in memory.
        CellBase(MeshData::Vertices const &verticesIn, Mesh const &origMesh, Dimensionless scaleIn =
                     1e0);

        //! \brief Initializes mesh from mesh data
        //! \param [in] mesh: deformable vertices that define the cell are copied
        //!    from this mesh. These values are *not* modified by the scale.
        //! \param [in] origMesh: Original mesh. A shallow copy is made of this
        //!    object.
        //! \param [in] scale: scales template by a given amount
        //!    The scale is added during internal operations. The template will still
        //!    refer to the same data in memory.
        CellBase(Mesh const &mesh, Mesh const &origMesh, Dimensionless scaleIn = 1e0);

        //! \brief Initializes mesh from mesh data
        //! \param [in] mesh: Modifyiable mesh and template. Deep copies are made of
        //!   both.
        CellBase(Mesh const &mesh);
        //! \brief Initializes mesh from mesh data
        //! \param [in] mesh: Modifyiable mesh and template. Deep copies are made of
        //!   both
        CellBase(std::shared_ptr<MeshData> const &mesh);
        //! Copy constructor
        //! References same template mesh
        CellBase(CellBase const &cell);

        //! Tag to choose shallow clone constructor
        struct shallow_clone
        {
        };
        //! Shallow copy constructor
        //! References same data
        CellBase(CellBase const &cell, shallow_clone const&) :
            data(cell.data)
        {
        }
        //! Because it is good practice
        virtual ~CellBase()
        {
        }

        void operator=(Mesh const &mesh);

        //! Unmodified mesh
        Mesh const &GetTemplateMesh() const;
        //! Facets for the mesh
        MeshData::Facets const &GetFacets() const;
        //! Vertices of the cell
        MeshData::Vertices const &GetVertices() const;
        //! Vertices of the cell
        MeshData::Vertices &GetVertices();
        //! Topology of the (template) mesh
        std::shared_ptr<MeshTopology const> GetTopology() const;
        site_t GetNumberOfNodes() const;

        //! Facet bending energy
        virtual PhysicalEnergy operator()() const = 0;
        //! Facet bending energy - pretty printing for shared ptr
        PhysicalEnergy Energy() const
        {
          return operator()();
        }
        //! Facet bending energy
        virtual PhysicalEnergy operator()(std::vector<LatticeForceVector> &in) const = 0;
        //! Facet bending energy - pretty printing for shared ptr
        PhysicalEnergy Energy(std::vector<LatticeForceVector> &in) const
        {
          return operator()(in);
        }
        //! Interaction between wall and a node
        virtual LatticeForceVector WallInteractionForce(LatticePosition const &vertex,
                                                        LatticePosition const &wall) const
        {
          return LatticeForceVector(0, 0, 0);
        }
        virtual bool HasWallForces() const
        {
          return false;
        }

        //! Scale mesh around barycenter
        void operator*=(Dimensionless const &);
        //! Linear transform of each vertex, centered around barycenter
        void operator*=(util::Matrix3D const &);
        //! Translate mesh
        void operator+=(LatticePosition const &offset);
        //! Translate mesh
        void operator-=(LatticePosition const &offset)
        {
          return operator+=(-offset);
        }
        //! Transform mesh
        void operator+=(std::vector<LatticePosition> const &displacements);

        MeshData::Vertices::value_type GetBarycenter() const;

        //! Scale to apply to the template mesh
        void SetScale(Dimensionless scaleIn);
        //! Scale to apply to the template mesh
        Dimensionless GetScale() const;

        // cloneImpl is virtual and returns a pointer to abstract class
        // clone will be overriden. It will call cloneImpl and cast it to derived type.
        // The loop jumping is necessary since we are returning managed pointers.
        //! Clones: shallow copy reference mesh, deep-copy everything else
        std::unique_ptr<CellBase> clone() const
        {
          return cloneImpl();
        }

      protected:
        //! Clones: shallow copy reference mesh, deep-copy everything else
        std::unique_ptr<CellBase> virtual cloneImpl() const = 0;
        //! allows separation of data and behaviors
        class CellData;
        //! Holds data
        std::shared_ptr<CellData> data;
    };

    //! Deformable cell for which energy and forces can be computed
    class Cell : public CellBase
    {
      public:
        //! Holds all physical parameters
        class Moduli
        {
          public:
            //! Bending energy parameter
            PhysicalPressure bending;
            //! Surface energy parameter
            PhysicalPressure surface;
            //! Volume energy parameter
            PhysicalPressure volume;
            //! Skalak dilation modulus
            PhysicalPressure dilation;
            //! Skalak strain modulus
            PhysicalPressure strain;

            Moduli(PhysicalPressure b = 0, PhysicalPressure s = 0, PhysicalPressure v = 0,
                   PhysicalPressure d = 0, PhysicalPressure st = 0) :
                bending(b), surface(s), volume(v), dilation(d), strain(s)
            {
            }
            Moduli(std::initializer_list<PhysicalPressure> const &l)
            {
              bending = l.size() > 0 ?
                *l.begin() :
                0;
              surface = l.size() > 1 ?
                * (l.begin() + 1) :
                0;
              volume = l.size() > 2 ?
                * (l.begin() + 2) :
                0;
              dilation = l.size() > 3 ?
                * (l.begin() + 3) :
                0;
              strain = l.size() > 4 ?
                * (l.begin() + 4) :
                0;
            }
        } moduli;
        //! Node-wall interaction
        Node2NodeForce nodeWall;

#       ifndef CPP11_HAS_CONSTRUCTOR_INHERITANCE
        Cell(MeshData::Vertices &&verticesIn, Mesh const &origMesh, Dimensionless scaleIn = 1e0) :
            CellBase(std::move(verticesIn), origMesh, scaleIn)
        {
        }
        Cell(MeshData::Vertices const &verticesIn, Mesh const &origMesh,
             Dimensionless scaleIn = 1e0) :
            CellBase(verticesIn, origMesh, scaleIn)
        {
        }
        Cell(Mesh const &mesh, Mesh const &origMesh, Dimensionless scaleIn = 1e0) :
            CellBase(mesh, origMesh, scaleIn)
        {
        }
        Cell(Mesh const &mesh) :
            CellBase(mesh)
        {
        }
        Cell(std::shared_ptr<MeshData> const &mesh) :
            CellBase(mesh)
        {
        }
        Cell(Cell const &cell, CellBase::shallow_clone const&) :
            CellBase(cell, CellBase::shallow_clone())
        {
        }
#       else
        // inheriting constructors
        using CellBase::CellBase;
#       endif
        //! Copy constructor
        //! Copy refers to the same template mesh
        Cell(Cell const &cell) :
            CellBase(cell), moduli(cell.moduli), nodeWall(cell.nodeWall)
        {
        }

        //! Facet bending energy
        virtual PhysicalEnergy operator()() const override;
        //! Facet bending energy
        virtual PhysicalEnergy operator()(std::vector<LatticeForceVector> &in) const override;
        //! Node-Wall interaction
        virtual LatticeForceVector WallInteractionForce(LatticePosition const &vertex,
                                                        LatticePosition const &wall) const override
        {
          return nodeWall(vertex, wall);
        }
        virtual bool HasWallForces() const override
        {
          return true;
        }

        //! Clones: shallow copy reference mesh, deep-copy everything else
        std::unique_ptr<Cell> clone() const
        {
          return std::unique_ptr<Cell>(static_cast<Cell*>(cloneImpl().release()));
        }

      private:
        //! Clones: shallow copy reference mesh, deep-copy everything else
        std::unique_ptr<CellBase> cloneImpl() const override;

        // Computes facet bending energy over all facets
        PhysicalEnergy facetBending() const;
        // Computes facet bending energy over all facets
        PhysicalEnergy facetBending(std::vector<LatticeForceVector> &forces) const;
    };

    //! Typical cell container type
    typedef std::set<std::shared_ptr<CellBase>> CellContainer;
  }
} // namespace hemelb::redblood
#endif
