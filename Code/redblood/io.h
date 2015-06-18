#ifndef HEMELB_REDBLOOD_IO_H
#define HEMELB_REDBLOOD_IO_H

#include <iostream>
#include <memory>
#include <functional>
#include "redblood/Cell.h"
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
    std::unique_ptr<CellBase> readCell(io::xml::Element const&, util::UnitConverter const&);
    //! Reads flow extension from XML
    FlowExtension readFlowExtension(io::xml::Element const&, util::UnitConverter const&);
    //! Reads all flow extensions from XML
    std::unique_ptr<std::vector<FlowExtension>> readFlowExtensions(
        io::xml::Element const& inletsNode, util::UnitConverter const& converter);
    //! Reads an RBCInserter from XML
    std::unique_ptr<RBCInserter> readSingleRBCInserter(
        io::xml::Element const&, util::UnitConverter const&);
    // std::function<void(CellInserter)> readRBCInserter(
    //     io::xml::Element const&, util::UnitConverter const&);
  }
}

#endif
