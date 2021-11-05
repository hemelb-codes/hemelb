// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "redblood/RBCConfig.h"

#include "configuration/SimConfig.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "redblood/FlowExtension.h"
#include "redblood/Node2Node.h"
#include "util/UnitConverter.h"

#ifdef HEMELB_BUILD_RBC
#include "redblood/xmlIO.h"
#endif

namespace hemelb {
  namespace redblood {
    using configuration::GetDimensionalValue;
    using configuration::GetDimensionalValueInLatticeUnits;

    FlowExtension readFlowExtension(io::xml::Element const& node,
                                    util::UnitConverter const& converter)
    {
      if (node == node.Missing())
      {
        throw Exception() << "Expected non-empty XML node";
      }

      auto const flowXML = node.GetChildOrThrow("flowextension");
      FlowExtension result;
      GetDimensionalValueInLatticeUnits(flowXML.GetChildOrThrow("length"), "m", converter, result.length);
      GetDimensionalValueInLatticeUnits(flowXML.GetChildOrThrow("radius"), "m", converter, result.radius);
      const auto fadeEl = flowXML.GetChildOrNull("fadelength");
      if (fadeEl) {
	GetDimensionalValueInLatticeUnits(fadeEl, "m", converter, result.fadeLength);
      }	else {
	result.fadeLength = result.length;
      }

      // Infer normal and position from inlet
      // However, normals point in *opposite* direction, and, as a result, origin are at opposite
      // end of the cylinder
      GetDimensionalValue(node.GetChildOrThrow("normal"), "dimensionless", result.normal);
      result.normal *= -1;
      result.normal.Normalise();

      {
	// Position has to be offset by origin, so can't use the simple converter here
	PhysicalPosition tmp;
	GetDimensionalValue(node.GetChildOrThrow("position"), "m", tmp);
	result.origin = converter.ConvertPositionToLatticeUnits(tmp);
      }
      result.origin -= result.normal * result.length;

      return result;
    }

    const Node2NodeForce& RBCConfig::GetCell2Cell() const
    {
      return *cell2Cell;
    }
    const Node2NodeForce& RBCConfig::GetCell2Wall() const
    {
      return *cell2Wall;
    }

    void RBCConfig::DoIOForRedBloodCells(const io::xml::Element& topNode,
					 const configuration::SimConfig& fullconfig,
					 util::UnitConverter const& units)
    {
#ifdef HEMELB_BUILD_RBC
      auto&& rbcNode = topNode.GetChildOrThrow("redbloodcells");
      const io::xml::Element controllerNode = rbcNode.GetChildOrThrow("controller");
      GetDimensionalValue(controllerNode.GetChildOrThrow("boxsize"), "lattice", boxSize);

      auto&& inletsNode = topNode.GetChildOrNull("inlets");
      rbcMeshes.reset(readTemplateCells(topNode, fullconfig, units).release());
      rbcinserter = readRBCInserters(inletsNode, units, *rbcMeshes);
      rbcOutlets = readRBCOutlets(topNode, units);
      cell2Cell = readNode2NodeForce(
          rbcNode.GetChildOrNull("cell2Cell"), units);
      cell2Wall = readNode2NodeForce(
          rbcNode.GetChildOrNull("cell2Wall"), units);
      if(boxSize < cell2Wall->cutoff)
      {
        throw Exception() << "Box-size < cell-wall interaction size: "
          "cell-wall interactions cannot be all accounted for.";
      }
      if(boxSize < cell2Cell->cutoff)
      {
        throw Exception() << "Box-size < cell-cell interaction size: "
          "cell-cell interactions cannot be all accounted for.";
      }
      const io::xml::Element outputNode = rbcNode.GetChildOrThrow("output");
      GetDimensionalValue(outputNode.GetChildOrThrow("period"), "lattice", rbcOutputPeriod);
#else
      throw Exception() << "Trying to use redblood cells when HEMELB_BUILD_RBC=OFF";
#endif
    }

  }
}
