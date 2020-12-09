// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_RBCCONFIG_H
#define HEMELB_REDBLOOD_RBCCONFIG_H

#include <functional>
#include <memory>
#include <vector>

#include "units.h"
#include "redblood/types_fwd.h"

namespace hemelb {
  namespace configuration {
    class SimConfig;
  }

  namespace io {
    namespace xml {
      class Element;
    }
  }
  namespace util {
    class UnitConverter;
  }

  namespace redblood {
    class FlowExtension;
    class Node2NodeForce;

    //! Reads flow extension from XML
    FlowExtension readFlowExtension(io::xml::Element const&, util::UnitConverter const&);

    // Hold all the configuration information from the XML concerning
    // the redblood module.
    class RBCConfig {
    public:
      using inserter_func_t = std::function<void(CellInserter const&)>;
    private:
      inserter_func_t rbcinserter;
      LatticeDistance boxSize;
      std::shared_ptr<TemplateCellContainer> rbcMeshes;
      std::shared_ptr<std::vector<FlowExtension>> rbcOutlets;
      std::shared_ptr<Node2NodeForce> cell2Cell;
      std::shared_ptr<Node2NodeForce> cell2Wall;
      LatticeTimeStep rbcOutputPeriod;

    public:
      // Called from configuration::SimConfig to handle our part of the XML.
      void DoIOForRedBloodCells(const io::xml::Element& topNode,
				const configuration::SimConfig& fullconfig,
				util::UnitConverter const&);

      inline std::shared_ptr<TemplateCellContainer> GetRBCMeshes() const
      {
	return rbcMeshes;
      }
      inline std::shared_ptr<std::vector<FlowExtension>> GetRBCOutlets() const
      {
	return rbcOutlets;
      }
      const Node2NodeForce& GetCell2Cell() const;
      const Node2NodeForce& GetCell2Wall() const;

      inline LatticeTimeStep const & GetRBCOutputPeriod() const
      {
	return rbcOutputPeriod;
      }
      /**
       * Returns the object used to insert red blood cells into the simulation.
       * @return
       */
      inline inserter_func_t GetInserter() const
      {
	return rbcinserter;
      }

      /**
       * Gets the box size for the RBC CellController.
       * @return
       */
      inline LatticeDistance GetBoxSize() const
      {
	return boxSize;
      }

    };
  }
}
#endif
