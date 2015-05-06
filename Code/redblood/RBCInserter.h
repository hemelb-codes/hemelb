#ifndef HEMELB_REDBLOOD_RBC_INSERTER_H
#define HEMELB_REDBLOOD_RBC_INSERTER_H

#include <iostream>
#include <memory>
#include <functional>
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
         * @param mesh_path the path to the mesh text file to read the cell
         * shape from
         * @param scale the scale of the cell to insert
         */
        RBCInserter(std::function<bool()> condition, const std::string & mesh_path,
                    std::vector<lb::iolets::InOutLet *> inlets =
                        std::vector<lb::iolets::InOutLet *>(),
                    Cell::Moduli moduli = Cell::Moduli(), Dimensionless scale = 1.0);

        /**
         * Creates an RBC Inserter.
         *
         * @param condition a cell will only be inserted on a LB step if this
         * condition evaluates to true
         * @param mesh_stream an input stream to read the cell shape from
         * @param scale the scale of the cell to insert
         */
        RBCInserter(std::function<bool()> condition, std::istream & mesh_stream,
                    std::vector<lb::iolets::InOutLet *> inlets =
                        std::vector<lb::iolets::InOutLet *>(),
                    Cell::Moduli moduli = Cell::Moduli(), Dimensionless scale = 1.0);

        /**
         * Creates an RBC Inserter.
         *
         * @param condition a cell will only be inserted on a LB step if this
         * condition evaluates to true
         * @param shape the shape of the cells to create
         * @param scale the scale of the cell to insert
         */
        RBCInserter(std::function<bool()> condition, const MeshData & shape,
                    std::vector<lb::iolets::InOutLet *> inlets =
                        std::vector<lb::iolets::InOutLet *>(),
                    Cell::Moduli moduli = Cell::Moduli(), Dimensionless scale = 1.0);

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
        void operator()(std::function<void(CellContainer::value_type)> insertFn) const;

        void SetShape(const MeshData & shape);
        std::shared_ptr<const MeshData> GetShape() const;

        void AddInLet(lb::iolets::InOutLet *);
        void RemoveInLet(lb::iolets::InOutLet *);

        void SetScale(Dimensionless scale);
        Dimensionless GetScale() const;

        void SetModuli(Cell::Moduli & moduli);
        const Cell::Moduli & GetModuli() const;

        void SetCondition(std::function<bool()> condition);
        std::function<bool()> GetCondition() const;

      private:

        //! When to insert cells
        std::function<bool()> condition;

        //! The shape of the cells to insert
        std::shared_ptr<const MeshData> shape;

        //! The iolets to insert cells into
        std::vector<lb::iolets::InOutLet *> inlets;

        //! The cell moduli
        Cell::Moduli moduli;

        //! The initial scale of the new cells
        Dimensionless scale;
    };

  }
}

#endif
