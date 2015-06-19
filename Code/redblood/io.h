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
    class RBCInserter;
    class FlowExtension;
    //! Reads cell from XML
    std::unique_ptr<CellBase> readCell(io::xml::Element const&, util::UnitConverter const&);
    //! Reads flow extension from XML
    FlowExtension readFlowExtension(io::xml::Element const&, util::UnitConverter const&);
    //! Reads all flow extensions from XML
    std::unique_ptr<std::vector<FlowExtension>> readFlowExtensions(
        io::xml::Element const& inletsNode, util::UnitConverter const& converter);
    //! Reads template meshes from XML
    std::unique_ptr<TemplateCellContainer> readTemplateCells(io::xml::Element const&,
                                                             util::UnitConverter const&);
    //! Reads all RBC inserters from XML
    std::function<void(CellInserter const&)> readRBCInserters(io::xml::Element const&,
                                                              util::UnitConverter const&,
                                                              TemplateCellContainer const&);
  }
}

#endif
