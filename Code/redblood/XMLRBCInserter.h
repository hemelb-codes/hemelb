#ifndef HEMELB_REDBLOOD_XML_RBC_INSERTER_H
#define HEMELB_REDBLOOD_XML_RBC_INSERTER_H

#include <iostream>
#include <memory>
#include "Mesh.h"
#include "Cell.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {

    /**
     * The XML RBC Inserter listens for red blood cells leaving outlets and
     * inserts an equal number of cells into an inlet.  The cell position and
     * scale are read from XML while the cell shape is read from a text file or
     * input stream.
     */
    class XMLRBCInserter
    {
      public:

        /**
         * Creates an XML RBC Inserter.
         *
         * @param xml_path the path to the XML file to read the cell scale and
         * position from
         * @param mesh_path the path to the mesh text file to read the cell
         * shape from
         */
        XMLRBCInserter(const std::string & xml_path, const std::string & mesh_path);

        /**
         * Creates an XML RBC Inserter.
         *
         * @param xml_path the path to the XML file to read the cell scale and
         * position from
         * @param mesh_stream an input stream to read the cell shape from
         */
        XMLRBCInserter(const std::string & xml_path, std::istream & mesh_stream);

        /**
         * Creates an XML RBC Inserter.
         *
         * @param xml_path the path to the XML file to read the cell scale and
         * position from
         * @param shape the shape of the cells to create
         */
        XMLRBCInserter(const std::string & xml_path, const MeshData & shape);

        /**
         * Creates an XML RBC Inserter.
         *
         * @param position the position to insert new cells into the simulation
         * @param scale the initial scale of the new cells
         * @param shape the shape of the cells to create
         */
        XMLRBCInserter(const MeshData::Vertices::value_type & position,
                       Dimensionless scale, const MeshData & shape);

        /**
         * Cell insertion callback called on each step of the simulation.  If
         * a cell has left the flow domain through an outlet then this method
         * will add a new cell via an inlet.
         *
         * @param insertFn the function to insert a new cell into the simulation
         *
         * @see hemelb::redblood::CellArmy::SetCellInsertion
         * @see hemelb::redblood::CellArmy::CallCellInsertion
         */
        void InsertCell(std::function<void(CellContainer::value_type)> insertFn);

        void SetShape(const MeshData & shape);
        std::shared_ptr<MeshData> GetShape();
        std::shared_ptr<const MeshData> GetShape() const;

        void SetPosition(const MeshData::Vertices::value_type & position);
        MeshData::Vertices::value_type & GetPosition();
        const MeshData::Vertices::value_type & GetPosition() const;

        void SetScale(Dimensionless scale);
        Dimensionless GetScale() const;

      private:
        /**
         * Reads the position and scale from an XML file.
         *
         * @param [in]  xml_path the xml file to read
         * @param [out] position the position to insert cells into the simulation
         * @param [out] scale the initial scale of the new cells
         */
        void read_position_and_scale_from_xml(const std::string & xml_path,
                                              MeshData::Vertices::value_type & position,
                                              Dimensionless & scale);

        //! The shape of the cells to insert
        std::shared_ptr<MeshData> shape;

        //! The position to insert cells into the simulation
        MeshData::Vertices::value_type position;

        //! The initial scale of the new cells
        Dimensionless scale;
    };

  }
}

#endif
