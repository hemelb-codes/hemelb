#include "redblood/Cell.h"
#include "redblood/FaderCell.h"
#include "redblood/RBCInserter.h"
#include "redblood/FlowExtension.h"
#include "redblood/io.h"

namespace hemelb
{
  namespace redblood
  {
    //! Throws if input does not have units
    template<typename T>
    T GetDimensionalValue(const io::xml::Element& elem, const std::string& units)
    {
      T value;
      const std::string& got = elem.GetAttributeOrThrow("units");
      if (got != units)
      {
        throw Exception() << "Invalid units for element " << elem.GetPath() << ". Expected '"
            << units << "', got '" << got << "'";
      }

      elem.GetAttributeOrThrow("value", value);
      return value;
    }
    //! Defaults to some value if parent, or its child elemname are not present
    template<typename T>
    T GetDimensionalValue(
        const io::xml::Element& parent, const std::string &elemname, const std::string& units,
        T default_)
    {
      if(parent == parent.Missing())
      {
        return default_;
      }
      auto const element = parent.GetChildOrNull(elemname);
      if(element == element.Missing())
      {
        return default_;
      }
      return GetDimensionalValue<T>(element, units);
    }
    //! Gets value and convert to LB units. Default value should be in physical units.
    template<typename T>
    T GetDimensionalValue(
        const io::xml::Element& parent, const std::string &elemname, const std::string& units,
        util::UnitConverter const &converter, T default_)
    {
      T const value = GetDimensionalValue<T>(parent, elemname, units, default_);
      return units == "LB" ? value: converter.ConvertToLatticeUnits(units, value);
    }
    //! Gets value and convert to LB units
    template<typename T>
    T GetDimensionalValue(
        const io::xml::Element& parent, const std::string &elemname, const std::string& units,
        util::UnitConverter const &converter)
    {
      T const value = GetDimensionalValue<T>(parent.GetChildOrThrow(elemname), units);
      return units == "LB" ? value: converter.ConvertToLatticeUnits(units, value);
    }
    //! Gets position and convert to LB units
    LatticePosition GetPosition(
        const io::xml::Element& parent, const std::string &elemname,
        util::UnitConverter const &converter)
    {
      auto value = GetDimensionalValue<LatticePosition>(parent.GetChildOrThrow(elemname), "m");
      return converter.ConvertPositionToLatticeUnits(value);
    }


    Cell::Moduli readModuli(io::xml::Element const& node, util::UnitConverter const &converter)
    {
      redblood::Cell::Moduli moduli;
      auto const moduliNode = node.GetChildOrNull("moduli");
      moduli.bending = GetDimensionalValue(moduliNode, "bending", "Nm", converter, 2e-19);
      moduli.surface = GetDimensionalValue(moduliNode, "surface", "LB", converter, 1e0);
      moduli.volume = GetDimensionalValue(moduliNode, "volume", "LB", converter, 1e0);
      moduli.dilation = GetDimensionalValue(moduliNode, "dilation", "N/m", converter, 5e-1);
      moduli.strain = GetDimensionalValue(moduliNode, "strain", "N/m", converter, 5e-6);
      return moduli;
    }

    Node2NodeForce readNode2NodeForce(
        io::xml::Element const& parent, util::UnitConverter const & converter)
    {
      Node2NodeForce result(
          1e0/converter.ConvertToLatticeUnits("N", 1e0), converter.GetVoxelSize(), 2);
      if(parent == parent.Missing())
      {
        return result;
      }
      auto const node = parent.GetChildOrNull("interaction");
      result.intensity = GetDimensionalValue(node, "intensity", "N", converter, result.intensity);
      result.cutoff = GetDimensionalValue(node, "cutoffdistance", "m", converter, result.cutoff);
      auto const exponentNode = node != node.Missing() ?
        node.GetChildOrNull("exponent"): node.Missing();
      if(exponentNode != node.Missing())
      {
        exponentNode.GetAttributeOrThrow("value", result.exponent);
      }
      return result;
    }

    std::unique_ptr<CellBase> readCell(
        io::xml::Element const& node, util::UnitConverter const& converter)
    {
      // auto const node = topNode.GetChildOrThrow("redbloodcells");
      if(node  == node.Missing())
      {
        throw Exception() << "Expected non-empty XML node";
      }
      const io::xml::Element & cellNode = node.GetChildOrThrow("cell");
      std::string const mesh_path
        = cellNode.GetChildOrThrow("shape").GetAttributeOrThrow("mesh_path");
      auto const mesh_data = read_mesh(mesh_path);
      auto const scale = GetDimensionalValue<LatticeDistance>(cellNode, "scale", "m", converter);
      std::unique_ptr<Cell> cell(new Cell(mesh_data->vertices, Mesh(mesh_data), scale));
      *cell *= scale;
      cell->moduli = readModuli(cellNode, converter);
      cell->nodeWall = readNode2NodeForce(cellNode, converter);

      std::unique_ptr<CellBase> cellbase(static_cast<CellBase*>(cell.release()));
      return cellbase;
      // // Check if there are iolets with flow extensions
      // auto const flowExtensions = readFlowExtensions(topNode.GetChildOrNull("inlets"), converter);
      //
      // return flowExtensions ? FaderCell(std::move(cellbase), flowExtensions).clone(): cellbase;
    }

    std::unique_ptr<std::vector<FlowExtension>> readFlowExtensions(
        io::xml::Element const& inletsNode, util::UnitConverter const& converter)
    {
      std::vector<FlowExtension> result;
      auto inletNode = inletsNode.GetChildOrNull("inlet");
      for(; inletNode != inletNode.Missing(); inletNode = inletNode.NextSiblingOrNull("inlet"))
      {
      }
      if(result.size() == 0)
      {
        return nullptr;
      }
      return std::unique_ptr<decltype(result)>(new decltype(result)(std::move(result)));
    }

    FlowExtension readFlowExtension(
        io::xml::Element const& node, util::UnitConverter const& converter)
    {
      if(node  == node.Missing())
      {
        throw Exception() << "Expected non-empty XML node";
      }

      auto const flowXML = node.GetChildOrThrow("flowextension");
      FlowExtension result;

      result.length = GetDimensionalValue<PhysicalDistance>(flowXML, "length", "m", converter);
      result.radius = GetDimensionalValue<PhysicalDistance>(flowXML, "radius", "m", converter);
      result.fadeLength = GetDimensionalValue<PhysicalDistance>(
          flowXML, "fadelength", "m", converter,
          converter.ConvertDistanceToPhysicalUnits(result.length));

      // Infer normal and position from inlet
      // However, normals point in *opposite* direction, and, as a result, origin are at opposite
      // end of the cylinder
      result.normal = -GetDimensionalValue<LatticePosition>(
          node, "normal", "dimensionless", converter);
      result.normal.Normalise();
      result.origin = GetPosition(node, "position", converter) - result.normal * result.length;

      return result;
    }


    //! Finds first inlet with Cell insertion
    io::xml::Element findFirstInletWithCellInsertion(io::xml::Element const &inlet)
    {
      io::xml::Element sibling(inlet);
      for(; sibling.Missing() != sibling; sibling = sibling.NextSiblingOrNull("inlet"))
      {
        auto hasRequisiteXMLTags =
            sibling.GetChildOrNull("flowextension") != sibling.Missing()
            and sibling.GetChildOrNull("insertcell") != sibling.Missing();
        if(hasRequisiteXMLTags)
        {
          return sibling;
        }
      }
      return sibling.Missing();
    }

    // std::function<void(CellInserter)> readRBCInserter(
    //     io::xml::Element const& node, util::UnitConverter const& converter)
    // {
    //   // Check if we have any cell insertion action going on at all
    //   auto const inlets = node.GetChildOrThrow("inlets");
    //   auto inlet = findFirstInletWithCellInsertion(inlets.GetChildOrThrow("inlet"));
    //   if(inlet.Missing() == inlet)
    //   {
    //     return nullptr;
    //   }
    // }

    std::unique_ptr<RBCInserter> readSingleRBCInserter(
        io::xml::Element const& node, util::UnitConverter const& converter)
    {
      // First look for  input with a flow extension and iteration
      auto const inlets = node.GetChildOrThrow("inlets");
      auto inlet = findFirstInletWithCellInsertion(inlets.GetChildOrThrow("inlet"));
      if(inlet.Missing() == inlet)
      {
        return nullptr;
      }
      // Check there is only one valid inlet for inserting cells
      if(findFirstInletWithCellInsertion(inlet.NextSiblingOrNull("inlet")) != inlet.Missing())
      {
        throw Exception() << "Can't handle multiple inlets dropping cells";
      }
      // Now gets data for cell insertion
      auto const insNode = inlet.GetChildOrThrow("insertcell");
      auto const deltaTime = GetDimensionalValue<PhysicalTime>(insNode, "every", "s", converter);
      auto const offset = GetDimensionalValue<PhysicalTime>(insNode, "offset", "s", converter, 0);
      auto const flowExtension = readFlowExtension(inlet, converter);
      auto cellbase = readCell(node.GetChildOrThrow("redbloodcells"), converter);
      std::unique_ptr<Cell> cell(static_cast<Cell*>(cellbase.release()));

      // Figure out size of cell alongst cylinder axis
      auto const barycenter = cell->GetBarycenter();
      auto minExtent = [barycenter, &flowExtension](LatticePosition const pos)
      {
        return std::max((pos - barycenter).Dot(flowExtension.normal), 0e0);
      };
      auto const minZ = *std::min_element(
          cell->GetVertices().begin(), cell->GetVertices().end(),
          [&minExtent](LatticePosition const &a, LatticePosition const& b)
          {
            return minExtent(a) < minExtent(b);
          }
      );
      // Place cell as close as possible to 0 of fade length
      *cell += flowExtension.origin
        + flowExtension.normal * (flowExtension.fadeLength - minExtent(minZ))
        - barycenter;

      // fail if any node outside flow extension
      for(auto const &vertex: cell->GetVertices())
      {
        if(not contains(flowExtension, vertex))
        {
          throw Exception() << "BAD INPUT: Cell not contained within flow extension";
        }
      }

      // Drops first cell when time reaches offset, and then every deltaTime thereafter.
      auto condition = [deltaTime, offset]()
      {
        static PhysicalTime time
          = deltaTime - 1e0 + std::numeric_limits<PhysicalTime>::epsilon() - offset;
        time += 1e0;
        if(time >= deltaTime)
        {
          time -= deltaTime;
          return true;
        }
        return false;
      };
      return std::unique_ptr<RBCInserter>(new RBCInserter(condition, std::move(cell)));
    }
  }
}
