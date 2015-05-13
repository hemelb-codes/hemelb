//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_BUFFER_COLUMNS_H
#define HEMELB_REDBLOOD_BUFFER_COLUMNS_H

#include <memory>

#include "redblood/Cell.h"
#include "redblood/Interpolation.h"
#include "redblood/FlowExtension.h"
#include "Exception.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    namespace buffer
    {
      //! Iterates over positions of cells arranged as columns in a cylinder
      class ColumnPositionIterator : private IndexIterator
      {
        public:
          //! Renders columns of cells
          //! \param[in] cylinder: Cylinder in which to pack the columns
          //! \param[in] vertices: vertices of the cell that will be dropped
          //!   Used to figure out size of the cell.
          //! \param[in] cellAxis: In conjunctions with colAxis, defines a rotation of the template
          //!      mesh. In practice, it should be vector going from one dimple to the other.
          //! \param[in] colAxis: Axis of the column of cells.
          //! \param[in] separation: separation between the cells. Defaults to 10% of smallest
          //! distance.
          ColumnPositionIterator(
              std::shared_ptr<Cylinder const> cylinder, MeshData::Vertices const& vertices,
              LatticePosition cellAxis, LatticePosition colAxis,
              LatticeDistance wallSeparation);
          //! Destroys buffer
          virtual ~ColumnPositionIterator()
          {
          }

          //! Goes to next cell position
          void operator++();
          //! Current position
          LatticePosition operator*() const
          {
            return major * static_cast<LatticeDistance>(current.z)
                + minor * static_cast<LatticeDistance>(current.y)
                + depth * static_cast<LatticeDistance>(current.x)
                + cylinder->origin;
          }

        protected:
          //! Cylinder within which columns are placed
          std::shared_ptr<Cylinder const> cylinder;
          //! Rotation matrix from template cell to cell as stacked in columns
          util::Matrix3D cellRotation;
          //! Grid spacings
          LatticePosition major, minor, depth;
          //! Spacing between wall and cells
          LatticePosition spacing;

          bool IsInside() const;
      };

      class ColumnCellDrop
      {
        public:
          ColumnCellDrop(
              std::shared_ptr<Cylinder const> cylinder, CellContainer::value_type cell,
              LatticePosition cellAxis, LatticePosition colAxis,
              LatticeDistance wallSeparation);

          //! Returns next cell
          CellContainer::value_type operator()();

        protected:
          //! Iterates over cell positions
          ColumnPositionIterator iterator;
          //! Reference to template cell.
          std::unique_ptr<CellBase> templateCell;
      };

    }
  }
} // namespace hemelb::redblood
#endif

