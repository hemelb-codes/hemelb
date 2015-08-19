#include <utility>
#include "redblood/Cell.h"
#include "redblood/FaderCell.h"
#include "redblood/RBCInserter.h"
#include "redblood/FlowExtension.h"
#include "redblood/io.h"

namespace hemelb
{
  namespace redblood
  {
    namespace
    {
      //! Throws if input does not have units
      template<typename T>
      T GetNonDimensionalValue(const io::xml::Element& elem, const std::string& units)
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
      T GetNonDimensionalValue(const io::xml::Element& parent, const std::string &elemname,
                            const std::string& units, T default_)
      {
        if (parent == parent.Missing())
        {
          return default_;
        }
        auto const element = parent.GetChildOrNull(elemname);
        if (element == element.Missing())
        {
          return default_;
        }
        return GetNonDimensionalValue<T>(element, units);
      }
      //! Gets value and convert to LB units. Default value should be in physical units.
      template<typename T>
      T GetNonDimensionalValue(const io::xml::Element& parent, const std::string &elemname,
                            const std::string& units, util::UnitConverter const &converter,
                            T default_)
      {
        T const value = GetNonDimensionalValue<T>(parent, elemname, units, default_);
        return units == "LB" ?
          value :
          converter.ConvertToLatticeUnits(units, value);
      }
      //! Gets value and convert to LB units
      template<typename T>
      T GetNonDimensionalValue(const io::xml::Element& parent, const std::string &elemname,
                            const std::string& units, util::UnitConverter const &converter)
      {
        T const value = GetNonDimensionalValue<T>(parent.GetChildOrThrow(elemname), units);
        return units == "LB" ?
          value :
          converter.ConvertToLatticeUnits(units, value);
      }
      //! Gets position and convert to LB units
      LatticePosition GetPosition(const io::xml::Element& parent, const std::string &elemname,
                                  util::UnitConverter const &converter)
      {
        auto value = GetNonDimensionalValue<LatticePosition>(parent.GetChildOrThrow(elemname), "m");
        return converter.ConvertPositionToLatticeUnits(value);
      }

      void readFlowExtensions(io::xml::Element const& ioletsNode,
                              util::UnitConverter const& converter,
                              std::vector<FlowExtension> &results,
                              bool mustHaveFlowExtension = false)
      {
        if (ioletsNode == ioletsNode.Missing())
        {
          return;
        }
        auto const name = ioletsNode.GetName().substr(0, ioletsNode.GetName().size() - 1);
        auto ioletNode = ioletsNode.GetChildOrNull(name);
        for (; ioletNode != ioletNode.Missing(); ioletNode = ioletNode.NextSiblingOrNull(name))
        {
          if (ioletNode.GetChildOrNull("flowextension") != ioletNode.Missing())
          {
            results.emplace_back(readFlowExtension(ioletNode, converter));
          }
          else if(mustHaveFlowExtension)
          {
            throw Exception() << "Could not find flow extension in iolet";
          }
        }
      }

      std::function<void(CellInserter const&)> readSingleRBCInserter(
          io::xml::Element const& node, util::UnitConverter const& converter,
          TemplateCellContainer const &templateCells)
      {
        // Now gets data for cell insertion
        auto const insNode = node.GetChildOrThrow("insertcell");
        auto const offset = GetNonDimensionalValue<LatticeTime>(insNode, "offset", "s", converter, 0);
        auto const templateName = insNode.GetAttributeOrThrow("template");
        if (templateCells.count(templateName) == 0)
        {
          throw Exception() << "Template cell name does not match a known template cell";
        }
        auto cell = templateCells.find(templateName)->second->clone();
        auto const flowExtension = readFlowExtension(node, converter);

        // Rotate cell to align z axis with given position, and then z axis with flow
        // If phi == 0, then cell symmetry axis is aligned with the flow
        using std::cos; using std::sin;
        auto const theta = GetNonDimensionalValue<Angle>(insNode, "theta", "rad", converter, 0e0);
        auto const phi = GetNonDimensionalValue<Angle>(insNode, "theta", "rad", converter, 0e0);
        LatticePosition const z(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
        auto const rotateToFlow = rotationMatrix(LatticePosition(0, 0, 1), flowExtension.normal);
        auto const rotation = rotateToFlow * rotationMatrix(LatticePosition(0, 0, 1), z);
        *cell *= rotation;

        // Figure out size of cell alongst cylinder axis
        auto const barycenter = cell->GetBarycenter();
        auto maxExtent = [barycenter, &flowExtension](LatticePosition const pos)
        {
          return std::max((pos - barycenter).Dot(flowExtension.normal), 0e0);
        };
        auto const maxZ =
            *std::max_element(cell->GetVertices().begin(),
                              cell->GetVertices().end(),
                              [&maxExtent](LatticePosition const &a, LatticePosition const& b)
                              {
                                return maxExtent(a) < maxExtent(b);
                              });
        // Place cell as close as possible to 0 of fade length
        *cell += flowExtension.origin
            + flowExtension.normal * (flowExtension.fadeLength - maxExtent(maxZ)) - barycenter;

        // fail if any node outside flow extension
        for (auto const &vertex : cell->GetVertices())
        {
          if (not contains(flowExtension, vertex))
          {
            HEMELB_CAPTURE(flowExtension.normal);
            HEMELB_CAPTURE(flowExtension.origin);
            HEMELB_CAPTURE(flowExtension.radius);
            HEMELB_CAPTURE(flowExtension.length);
            HEMELB_CAPTURE(vertex);
            throw Exception() << "BAD INPUT: Cell not contained within flow extension";
          }
        }

        // Drops first cell when time reaches offset, and then every deltaTime thereafter.
        // Note: c++14 will allow more complex captures. Until then, we will need to create
        // semi-local lambda variables on the stack as shared pointers. Where semi-local means the
        // variables should live as long as the lambda. But longuer than a single call.
        auto const timeStep = GetNonDimensionalValue<LatticeTime>(insNode, "every", "s", converter);
        auto const dt = GetNonDimensionalValue<LatticeTime>(insNode, "delta_t", "s", converter, 0e0);
        auto time = std::make_shared<LatticeTime>
        (
          timeStep - 1e0 + std::numeric_limits<LatticeTime>::epsilon() - offset
        );
        auto condition = [time, timeStep, dt, offset]()
        {
          *time += 1e0;
          if(*time >= timeStep)
          {
            *time -= timeStep + dt * (double(rand() % 10000) / 10000.e0);
            return true;
          }
          return false;
        };
        auto const dtheta
          = GetNonDimensionalValue<Angle>(insNode, "delta_theta", "rad", converter, 0e0);
        auto const dphi
          = GetNonDimensionalValue<Angle>(insNode, "delta_phi", "rad", converter, 0e0);
        auto const dx
          = GetNonDimensionalValue<LatticeDistance>(insNode, "delta_x", "m", converter, 0e0);
        auto const dy
          = GetNonDimensionalValue<LatticeDistance>(insNode, "delta_y", "m", converter, 0e0);
        return RBCInserterWithPerturbation
        (
            condition, std::move(cell),
            rotation,
            dtheta, dphi,
            rotateToFlow * LatticePosition(dx, 0, 0),
            rotateToFlow * LatticePosition(0, dy, 0)
        );
      }
    }

    Cell::Moduli readModuli(io::xml::Element const& node, util::UnitConverter const &converter)
    {
      redblood::Cell::Moduli moduli;
      auto const moduliNode = node.GetChildOrNull("moduli");
      moduli.bending = GetNonDimensionalValue(moduliNode, "bending", "Nm", converter, 2e-19);
      moduli.surface = GetNonDimensionalValue(moduliNode, "surface", "LB", converter, 1e0);
      moduli.volume = GetNonDimensionalValue(moduliNode, "volume", "LB", converter, 1e0);
      moduli.dilation = GetNonDimensionalValue(moduliNode, "dilation", "N/m", converter, 5e-1);
      moduli.strain = GetNonDimensionalValue(moduliNode, "strain", "N/m", converter, 5e-6);
      return moduli;
    }

    Node2NodeForce readNode2NodeForce(io::xml::Element const& parent,
                                      util::UnitConverter const & converter)
    {
      Node2NodeForce result(1e0 / converter.ConvertToLatticeUnits("N", 1e0), 1, 2);
      if (parent == parent.Missing())
      {
        return result;
      }
      auto const node = parent.GetChildOrNull("interaction");
      result.intensity = GetNonDimensionalValue(node, "intensity", "Nm", converter, result.intensity);
      result.cutoff = GetNonDimensionalValue(node, "cutoffdistance", "LB", converter, result.cutoff);
      auto const exponentNode = node != node.Missing() ?
        node.GetChildOrNull("exponent") :
        node.Missing();
      if (exponentNode != node.Missing())
      {
        exponentNode.GetAttributeOrThrow("value", result.exponent);
      }
      return result;
    }

    //! Reads multiple templates from XML and stores in container
    std::unique_ptr<TemplateCellContainer> readTemplateCells(io::xml::Element const& topNode,
                                                             util::UnitConverter const& converter)
    {
      std::unique_ptr<TemplateCellContainer> result(new TemplateCellContainer);
      // read flow extensions, if they exist
      std::shared_ptr<std::vector<FlowExtension>> flowExtensions(readFlowExtensions(topNode,
                                                                                    converter).release());
      // Then read template cells
      auto const cellsNode = topNode.GetChildOrThrow("redbloodcells").GetChildOrThrow("cells");
      auto cellNode = cellsNode.GetChildOrThrow("cell");
      for (; cellNode != cellNode.Missing(); cellNode = cellNode.NextSiblingOrNull("cell"))
      {
        auto const name = cellNode.GetAttributeOrNull("name");
        auto const key = name != nullptr ?
          *name :
          "default";
        if (result->count(key) != 0)
        {
          throw Exception() << "Multiple template mesh with same name";
        }
        auto cell = readCell(cellNode, converter);
        if (flowExtensions)
        {
          auto fader = FaderCell(std::move(cell), flowExtensions).clone();
          cell = std::move(fader);
        }
        result->emplace(key, std::shared_ptr<CellBase>(cell.release()));
      }
      return result->size() > 0 ?
        std::move(result) :
        nullptr;
    }

    std::unique_ptr<CellBase> readCell(io::xml::Element const& node,
                                       util::UnitConverter const& converter)
    {
      // auto const node = topNode.GetChildOrThrow("redbloodcells");
      if (node == node.Missing())
      {
        throw Exception() << "Expected non-empty XML node";
      }
      const auto cellNode = node.GetName() == "cell" ?
        node :
        node.GetChildOrThrow("cell");
      auto const name = cellNode.GetAttributeOrNull("name") == nullptr ?
        "default" :
        cellNode.GetAttributeOrThrow("name");
      std::string const mesh_path =
          cellNode.GetChildOrThrow("shape").GetAttributeOrThrow("mesh_path");
      auto const mesh_data = readMesh(mesh_path);
      auto const scale = GetNonDimensionalValue<LatticeDistance>(cellNode, "scale", "m", converter);
      std::unique_ptr<Cell> cell(new Cell(mesh_data->vertices, Mesh(mesh_data), scale, name));
      *cell *= scale;
      cell->moduli = readModuli(cellNode, converter);
      cell->nodeWall = readNode2NodeForce(cellNode, converter);

      std::unique_ptr<CellBase> cellbase(static_cast<CellBase*>(cell.release()));
      return cellbase;
    }

    std::unique_ptr<std::vector<FlowExtension>> readFlowExtensions(
        io::xml::Element const& topNode, util::UnitConverter const& converter)
    {
      if (topNode == topNode.Missing())
      {
        return nullptr;
      }
      std::vector<FlowExtension> result;
      auto inletsNode = topNode.GetChildOrNull("inlets");
      if (inletsNode != inletsNode.Missing())
      {
        readFlowExtensions(inletsNode, converter, result);
      }
      auto outletsNode = topNode.GetChildOrNull("outlets");
      if (outletsNode != outletsNode.Missing())
      {
        readFlowExtensions(outletsNode, converter, result);
      }
      if (result.size() == 0)
      {
        return nullptr;
      }
      return std::unique_ptr<decltype(result)>(new decltype(result)(std::move(result)));
    }

    FlowExtension readFlowExtension(io::xml::Element const& node,
                                    util::UnitConverter const& converter)
    {
      if (node == node.Missing())
      {
        throw Exception() << "Expected non-empty XML node";
      }

      auto const flowXML = node.GetChildOrThrow("flowextension");
      FlowExtension result;

      result.length = GetNonDimensionalValue<LatticeDistance>(flowXML, "length", "m", converter);
      result.radius = GetNonDimensionalValue<LatticeDistance>(flowXML, "radius", "m", converter);
      result.fadeLength =
          GetNonDimensionalValue<LatticeDistance>(flowXML,
                                                "fadelength",
                                                "m",
                                                converter,
                                                converter.ConvertDistanceToPhysicalUnits(result.length));

      // Infer normal and position from inlet
      // However, normals point in *opposite* direction, and, as a result, origin are at opposite
      // end of the cylinder
      result.normal = -GetNonDimensionalValue<LatticePosition>(node,
                                                            "normal",
                                                            "dimensionless",
                                                            converter);
      result.normal.Normalise();
      result.origin = GetPosition(node, "position", converter) - result.normal * result.length;

      return result;
    }

    //! Finds first inlet with Cell insertion
    io::xml::Element findFirstInletWithCellInsertion(io::xml::Element const &inlet)
    {
      io::xml::Element sibling(inlet);
      for (; sibling.Missing() != sibling; sibling = sibling.NextSiblingOrNull("inlet"))
      {
        auto hasRequisiteXMLTags = sibling.GetChildOrNull("flowextension") != sibling.Missing()
            and sibling.GetChildOrNull("insertcell") != sibling.Missing();
        if (hasRequisiteXMLTags)
        {
          return sibling;
        }
      }
      return sibling.Missing();
    }

    std::function<void(CellInserter const&)> readRBCInserters(
        io::xml::Element const& node, util::UnitConverter const& converter,
        TemplateCellContainer const& templateCells)
    {
      // Check if we have any cell insertion action going on at all
      auto const inlets = node.GetChildOrThrow("inlets");
      auto inlet = findFirstInletWithCellInsertion(inlets.GetChildOrThrow("inlet"));
      std::vector<std::function<void(CellInserter const&)>> results;
      while (inlet != inlet.Missing())
      {
        results.emplace_back(std::move(readSingleRBCInserter(inlet, converter, templateCells)));
        inlet = findFirstInletWithCellInsertion(inlet.NextSiblingOrNull("inlet"));
      }

      // do all insertion functions in one go
      if (results.size() > 1)
      {
        // go to shared pointer to avoid copies
        std::shared_ptr<decltype(results)> functions(new decltype(results)(std::move(results)));
        // return a lambda that loops over functions
        return [functions](CellInserter const& inserter)
        {
          for(auto const & func: *functions)
          {
            func(inserter);
          }
        };
      }
      // return null if no insertion, and just the function if only one
      return results.size() == 0 ?
        nullptr :
        results.front();
    }

    std::shared_ptr<std::vector<FlowExtension>> readRBCOutlets(io::xml::Element const& topNode,
                                                               util::UnitConverter const& converter)
    {
      // First read outlets from XML
      auto const result = std::make_shared<std::vector<FlowExtension>>();
      auto outletsNode = topNode.GetChildOrThrow("outlets");
      readFlowExtensions(outletsNode, converter, *result, true);
      // Then transforms them to cell outlets: should start somewhere near the end of fadelength
      for(auto &flowExt: *result)
      {
        // Minimum length should probably be sufficiently larger than cells.
        auto length = std::min
        (
            flowExt.length,
            std::max(flowExt.length - 0.9*flowExt.fadeLength, 2.)
        );
        flowExt.origin += flowExt.normal * (flowExt.length - length);
        flowExt.length = length;
        flowExt.fadeLength = length;
      }
      return result;
    }
  }
}
