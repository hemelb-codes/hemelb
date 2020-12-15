// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_CELL_H
#define HEMELB_REDBLOOD_CELL_H

#include <boost/uuid/uuid.hpp>

#include <set>
#include <utility>

#include "redblood/CellBase.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    class Mesh;

    //! Deformable cell for which energy and forces can be computed
    class Cell : public CellBase
    {
      public:
        //! Holds all physical parameters
        class Moduli
        {
          public:
            //! Bending energy parameter
            LatticeModulus bending;
            //! Surface energy parameter
            LatticeModulus surface;
            //! Volume energy parameter
            LatticeModulus volume;
            //! Skalak dilation modulus
            LatticeModulus dilation;
            //! Skalak strain modulus
            LatticeModulus strain;

            Moduli(LatticeModulus b = 0, LatticeModulus s = 0, LatticeModulus v = 0,
                   LatticeModulus d = 0, LatticeModulus st = 0) :
                bending(b), surface(s), volume(v), dilation(d), strain(st)
            {
            }
            Moduli(std::initializer_list<LatticeModulus> const &l)
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

#       ifndef CPP11_HAS_CONSTRUCTOR_INHERITANCE
        Cell(MeshData::Vertices &&verticesIn, Mesh const &origMesh, LatticeDistance scale = 1e0,
             std::string const & templateName = "default") :
            CellBase(std::move(verticesIn), origMesh, scale, templateName)
        {
        }
        Cell(MeshData::Vertices const &verticesIn, Mesh const &origMesh,
             LatticeDistance scale = 1e0, std::string const & templateName = "default") :
            CellBase(verticesIn, origMesh, scale, templateName)
        {
        }
        Cell(Mesh const &mesh, Mesh const &origMesh, LatticeDistance scale = 1e0,
             std::string const & templateName = "default") :
            CellBase(mesh, origMesh, scale, templateName)
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
            CellBase(cell), moduli(cell.moduli)
        {
        }

        //! Facet bending energy
        virtual LatticeEnergy operator()() const override;
        //! Facet bending energy
        virtual LatticeEnergy operator()(std::vector<LatticeForceVector> &in) const override;

        //! Clones: shallow copy reference mesh, deep-copy everything else
        std::unique_ptr<Cell> clone() const
        {
          return std::unique_ptr<Cell>(static_cast<Cell*>(cloneImpl().release()));
        }

        // Computes facet bending energy over all facets
        LatticeEnergy facetBending() const;
        // Computes facet bending energy over all facets
        LatticeEnergy facetBending(std::vector<LatticeForceVector> &forces) const;

      private:
        //! Clones: shallow copy reference mesh, deep-copy everything else
        std::unique_ptr<CellBase> cloneImpl() const override;
    };
    static_assert(
        (not std::is_default_constructible<Cell>::value)
        and (not std::is_nothrow_default_constructible<Cell>::value)
        and std::is_move_constructible<Cell>::value
        and (not std::is_nothrow_move_constructible<Cell>::value)
        and std::is_copy_constructible<Cell>::value
        and std::is_copy_assignable<Cell>::value
        and std::is_nothrow_copy_assignable<Cell>::value
        and (not std::is_standard_layout<Cell>::value)
        and (not std::is_pod<Cell>::value)

        and std::is_default_constructible<Cell::Moduli>::value
        and (not std::is_nothrow_default_constructible<Cell::Moduli>::value)
        and std::is_move_constructible<Cell::Moduli>::value
        and std::is_nothrow_move_constructible<Cell::Moduli>::value
        and std::is_copy_constructible<Cell::Moduli>::value
        and std::is_copy_assignable<Cell::Moduli>::value
        and std::is_nothrow_copy_assignable<Cell::Moduli>::value
        and std::is_standard_layout<Cell::Moduli>::value
        and (not std::is_pod<Cell::Moduli>::value),
        "Explicit type characteristics"
    );

  }
} // namespace hemelb::redblood
#endif
