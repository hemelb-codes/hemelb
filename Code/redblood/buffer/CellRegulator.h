//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_BUFFER_CELLREGULATOR_H
#define HEMELB_REDBLOOD_BUFFER_CELLREGULATOR_H

#include <vector>
#include "redblood/Cell.h"
#include "redblood/FlowExtension.h"
#include "redblood/buffer/Buffer.h"

namespace hemelb
{
  namespace redblood
  {
    namespace buffer
    {
      //! \brief Regulates the number of cells in simulation
      //! \details This mixin is responsible for making requests to buffers for new cells.
      //! This is implemented as a *Mix-in*. Currently, its sole requirement is that its derived
      //! class contains an addCell function taking a shared pointer to a cell.
      template<class DERIVED> class CellRegulator
      {
        public:
          //! \brief Signature of functions queried about cell insertion
          //! \details For constant flow, this function will return 1 once every x call, and 0
          //! otherwise.
          typedef std::function<site_t()> FlowFunction;
          //! Signature of a function returning a new cell at each call
          typedef Buffer::CellDistributionFunction CellDistributionFunction;
          //! \brief Regulates cells coming from X inlets
          CellRegulator();
          //! \brief virtual because style points
          virtual ~CellRegulator();

          //! Adds a new inlet
          void addInlet(
              std::shared_ptr<Cylinder>, CellDistributionFunction const &, FlowFunction const &);

          //! Makes request to the buffers
          void requestNewCells();
          //! drop all cells that can be dropped
          void dropCells();

        private:
          //! Container of inlets where cells are dropped
          std::vector<Buffer> buffers;
          //! Container of inlets where cells are dropped
          std::vector<FlowFunction> flows;
      };

      template<class DERIVED> CellRegulator<DERIVED>::CellRegulator()
      {
      }

      template<class DERIVED>
        void CellRegulator<DERIVED>::addInlet(
            std::shared_ptr<Cylinder> cyl,
            CellDistributionFunction const &density,
            FlowFunction const &flow)
        {
          buffers.emplace_back(cyl);
          buffers.back().SetCellDistributionFunction(density);
          flows.emplace_back(flow);
        }

      template<class DERIVED>
        void CellRegulator<DERIVED>::requestNewCells()
        {
          for(size_t i(0); i < buffers.size(); ++i)
          {
            if(site_t const n = flows[i]())
            {
              buffers[i].requestNewCells(n);
            }
          }
        }

      template<class DERIVED>
        void CellRegulator<DERIVED>::dropCells()
        {
          auto *const derived = static_cast<DERIVED*>(this);
          auto callback = [derived](CellContainer::value_type cell)
          {
            derived->AddCell(cell);
          };
          for(auto &buffer: buffers)
          {
            buffer(callback);
          }
        }
    } // buffer
  } // redblood
} // hemelb
#endif
