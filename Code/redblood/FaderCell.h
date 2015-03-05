//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_METACELL_H
#define HEMELB_REDBLOOD_METACELL_H

#include <list>
#include <functional>

#include "redblood/FlowExtension.h"
#include "redblood/Cell.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    //! A cell that computes forces differently depending on position
    //! Makes it easier to implement fade-in fade-out
    class FaderCell : public CellBase
    {
      public:
        //! Construct fader cell over other CellBase type object
        //! \param[in] wrappee object to which we are adding fade-in fade-out behaviour
        //! \param[in] flowExtensions iolet descriptions. A copy is made of this object.
        FaderCell(std::shared_ptr<CellBase> cell,
          std::vector<FlowExtension> const &flowExtensions = {})
          : CellBase(*cell, shallow_clone()),
            iolets(std::make_shared<std::vector<FlowExtension>>(flowExtensions)),
            wrappee(cell)
        {
        }
        FaderCell(std::shared_ptr<CellBase> cell,
          std::vector<FlowExtension> const &&flowExtensions)
          : CellBase(*cell, shallow_clone()),
            iolets(std::make_shared<std::vector<FlowExtension>>(std::move(flowExtensions))),
            wrappee(cell)
        {
        }
        //! Construct fader cell over other CellBase type object
        //! \param[in] wrappee object to which we are adding fade-in fade-out behaviour
        //! \param[in] flowExtensions iolet descriptions. A reference(shared pointer) is kept to
        //! this object.
        FaderCell(std::shared_ptr<CellBase> cell,
          std::shared_ptr<std::vector<FlowExtension>> flowExtensions)
          : CellBase(*cell, shallow_clone()), iolets(flowExtensions), wrappee(cell)
        {
        }
        // inheriting constructors
        using CellBase::CellBase;
        //! Facet bending energy
        virtual PhysicalEnergy operator()() const override;
        //! Facet bending energy
        virtual PhysicalEnergy operator()(std::vector<LatticeForceVector> &in) const override;
        //! Node-Wall interaction
        virtual LatticeForceVector WallInteractionForce(
            LatticePosition const &vertex, LatticePosition const &wall) const override
        {
          return wrappee->WallInteractionForce(vertex, wall);
        }

        std::shared_ptr<std::vector<FlowExtension> const> GetIOlets() const
        {
          return iolets;
        }
        std::shared_ptr<std::vector<FlowExtension>> GetIOlets()
        {
          return iolets;
        }

      private:
        //! Pointer to shared data
        std::shared_ptr<std::vector<FlowExtension>> iolets;
        //! Pointer to wrappee cell object
        std::shared_ptr<CellBase> wrappee;
    };
  }
} // namespace hemelb::redblood
#endif
