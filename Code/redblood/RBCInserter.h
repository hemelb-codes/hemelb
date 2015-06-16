#ifndef HEMELB_REDBLOOD_RBC_INSERTER_H
#define HEMELB_REDBLOOD_RBC_INSERTER_H

#include <iostream>
#include <memory>
#include <functional>
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/iolets/InOutLet.h"
#include "Mesh.h"
#include "Cell.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {

    /**
     * The RBC Inserter inserts cells into the simulation while some condition
     * evaluates to true.  The cell is inserted into an inlet or the start of a
     * flow extension (if the inlet has a flow extension defined).
     * The shape of the cell is read from a text file.
     */
    class RBCInserter
    {
      public:

        /**
         * Creates an RBC Inserter.
         *
         * @param condition a cell will only be inserted on a LB step if this
         * condition evaluates to true
         * @param cell to insert. Inserted as is. E.g. same position etc.
         * @param scale the scale of the cell to insert
         */
        RBCInserter(std::function<bool()> condition, std::unique_ptr<CellBase const> cell)
          : condition(condition), cell(std::move(cell))
        {
        }

        /**
         * Cell insertion callback called on each step of the simulation.  Cells
         * are inserted into the inlets while condition() evaluates to true.
         * Multiple cells are potentially inserted on each iteration.
         *
         * @param insertFn the function to insert a new cell into the simulation
         *
         * @see hemelb::redblood::CellArmy::SetCellInsertion
         * @see hemelb::redblood::CellArmy::CallCellInsertion
         */
        void operator()(CellInserter insertFn) const
        {
          if(condition())
          {
            insertFn(CellContainer::value_type(cell->clone().release()));
          }
        }

      private:
        //! When to insert cells
        std::function<bool()> condition;
        //! The shape of the cells to insert
        std::unique_ptr<CellBase const> cell;
    };


    //! Reads and inserts cell from XML
    std::shared_ptr<CellBase> read_cell(io::xml::Element const&, util::UnitConverter const&);
    // RBCInserter inserterFactory(io::xml:Element const *cellDef,  
  }
}

#endif
