//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_BUFFER_BUFFER_H
#define HEMELB_REDBLOOD_BUFFER_BUFFER_H

#include "redblood/Cell.h"
#include "redblood/FlowExtension.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    namespace buffer
    {
      //! Container of virtual cells
      class Buffer
      {
        public:
          //! Constructs buffer from geometry and cells
          //! \param[in] cyl: Geometric description of the buffer
          //!     The origin of the cylinder is the drop point in the coordinates of LB. The normal
          //!     should point from the region nearest to drop point to furthest from drop point. It
          //!     is co-linear with the axis of the cylinder.
          //! \param[in] cells: Current cells in the buffer
          Buffer(std::shared_ptr<Cylinder> cyl, CellContainer const& cells = CellContainer())
            : geometry(cyl), virtuals(cells), offset(cyl->origin)
          {
          }
          //! Constructs buffer from geometry only
          Buffer(Cylinder const & cyl, CellContainer const& cells = CellContainer())
            : geometry(new Cylinder(cyl)), virtuals(cells), offset(cyl.origin)
          {
          }
          //! Destroys buffer
          virtual ~Buffer()
          {
          }

          //! Inserts a cell. Does not check whether it is inside the geometry of the buffer, nor
          //! whether it overlaps with other cells.
          auto insert(CellContainer::value_type cell) -> decltype(CellContainer().insert(cell))
          {
            return virtuals.insert(cell);
          }

          //! Drops the next nearest cell
          CellContainer::value_type drop();

          LatticePosition const & GetOffset() const
          {
            return offset;
          }
          void SetOffset(LatticePosition const &off)
          {
            offset = off;
          }
#         ifdef HEMELB_DOING_UNITTESTS
          CellContainer::value_type GetJustDropped() const
          {
            return justDropped;
          }
#         endif


        protected:
          //! Geometry of the buffer
          std::shared_ptr<Cylinder> geometry;
          //! Container of virtual cells
          CellContainer virtuals;
          //! Offset between LB and buffer coordinate systems.
          LatticePosition offset;
          //! Last dropped
          CellContainer::value_type justDropped;
      };
    }
  }
} // namespace hemelb::redblood
#endif
