#ifndef HEMELB_REDBLOOD_IO_H
#define HEMELB_REDBLOOD_IO_H

#include <iostream>
#include <memory>
#include <functional>
#include "io/xml/XmlAbstractionLayer.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace redblood
  {
    class CellBase;
    class RBCInserter;
    class FlowExtension;
    //! Reads cell from XML
    std::unique_ptr<CellBase> read_cell(io::xml::Element const&, util::UnitConverter const&);
    //! Reads flow extension from XML
    FlowExtension read_flow_extension(io::xml::Element const&, util::UnitConverter const&);
    //! Reads an RBCInserter from XML
    std::unique_ptr<RBCInserter> read_rbcinserter(
        io::xml::Element const&, util::UnitConverter const&);
  }
}

#endif
