// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <utility>
#include <chrono>
#include "configuration/SimConfig.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "redblood/stencil.h"
#include "redblood/xmlIO.h"
#include "redblood/Cell.h"
#include "redblood/FaderCell.h"
#include "redblood/MeshIO.h"
#include "redblood/Node2Node.h"
#include "redblood/RBCConfig.h"
#include "redblood/RBCInserter.h"
#include "redblood/FlowExtension.h"
#include "Traits.h"

namespace hemelb
{
  namespace redblood
  {
    namespace
    {

      //! Throws if input does not have units
      template<typename T>
      std::pair<T, std::string> GetNonDimensionalValue(const io::xml::Element& elem,
                                                       const std::string& units)
      {
        T value;
        const std::string& got = elem.GetAttributeOrThrow("units");
        if (got != units && got != "lattice")
        {
          throw Exception() << "Invalid units for element " << elem.GetPath() << ". Expected '"
              << units << "', got '" << got << "'";
        }

        elem.GetAttributeOrThrow("value", value);
        return std::pair<T, std::string>(value, got);
      }

      //! Defaults to some value if parent, or its child elemname are not present
      template<typename T>
      std::pair<T, std::string> GetNonDimensionalValue(const io::xml::Element& parent,
                                                       const std::string &elemname,
                                                       const std::string& units, T default_)
      {
        if (parent == parent.Missing())
        {
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::Singleton>("Using internal default value for RBC parameter: %s\n",
                                                                               elemname.c_str());
          return std::pair<T, std::string>(default_, units);
        }
        auto const element = parent.GetChildOrNull(elemname);
        if (element == element.Missing())
        {
          hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::Singleton>("Using internal default value for RBC parameter: %s\n",
                                                                               elemname.c_str());
          return std::pair<T, std::string>(default_, units);
        }
        return GetNonDimensionalValue<T>(element, units);
      }

      //! Gets value and convert to LB units. Default value should be in physical units.
      template<typename T>
      T GetNonDimensionalValue(const io::xml::Element& parent, const std::string &elemname,
                               const std::string& units, util::UnitConverter const &converter,
                               T default_)
      {
        auto const value_unit = GetNonDimensionalValue<T>(parent, elemname, units, default_);
        return value_unit.second == "lattice" ?
          value_unit.first :
          converter.ConvertToLatticeUnits(value_unit.second, value_unit.first);
      }

      //! Gets value and convert to LB units
      template<typename T>
      T GetNonDimensionalValue(const io::xml::Element& parent, const std::string &elemname,
                               const std::string& units, util::UnitConverter const &converter)
      {
        auto const value_unit = GetNonDimensionalValue<T>(parent.GetChildOrThrow(elemname), units);
        return value_unit.second == "lattice" ?
          value_unit.first :
          converter.ConvertToLatticeUnits(value_unit.second, value_unit.first);
      }

      //! Gets position and convert to LB units
      LatticePosition GetPosition(const io::xml::Element& parent, const std::string &elemname,
                                  util::UnitConverter const &converter)
      {
        std::pair<LatticePosition, std::string> const value_unit = GetNonDimensionalValue<
            LatticePosition>(parent.GetChildOrThrow(elemname), "m");
        return value_unit.second == "lattice" ?
          value_unit.first :
          converter.ConvertPositionToLatticeUnits(value_unit.first);
      }

      void readFlowExtensions(io::xml::Element const& ioletsNode,
                              util::UnitConverter const& converter,
                              std::vector<FlowExtension> &results, bool mustHaveFlowExtension =
                                  false)
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
          else if (mustHaveFlowExtension)
          {
            throw Exception() << "Could not find flow extension in iolet";
          }
        }
      }

      // Read physical outlets that are somehow configured as numerical inlets in the XML parameter file
      void readFlowExtensionsWithoutInsertElement(io::xml::Element const& ioletsNode,
                              util::UnitConverter const& converter,
                              std::vector<FlowExtension> &results, bool mustHaveFlowExtension =
                                  false)
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
            if (ioletNode.GetChildOrNull("insertcell") != ioletNode.Missing())
            {
              continue;
            }
            else
            {
              results.emplace_back(readFlowExtension(ioletNode, converter));
            }
          }
          else if (mustHaveFlowExtension)
          {
            throw Exception() << "Could not find flow extension in iolet";
          }
        }
      }

      //! Rotates a cell to be aligned with the flow and translates it to the start of the flow extension fade length
      void rotateTranslateCellToFlow(std::unique_ptr<CellBase> & cell, const Angle theta,
                                     const Angle phi, const FlowExtension & flowExtension,
                                     LatticePosition const & translation,
                                     util::Matrix3D & rotateToFlow, util::Matrix3D & rotation)
      {
        // Rotate cell to align z axis with given position, and then z axis with flow
        // If phi == 0, then cell symmetry axis is aligned with the flow
        using std::cos;
        using std::sin;
        LatticePosition const z(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
        rotateToFlow = rotationMatrix(LatticePosition(0, 0, 1), flowExtension.normal);
        rotation = rotateToFlow * rotationMatrix(LatticePosition(0, 0, 1), z);
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
            + flowExtension.normal * (flowExtension.fadeLength - maxExtent(maxZ)) - barycenter
            + rotateToFlow * translation;

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
      }

      std::function<void(CellInserter const&)> readSingleRBCInserter(
          io::xml::Element const& node, util::UnitConverter const& converter,
          TemplateCellContainer const &templateCells)
      {
        // We need to seed each of the RBCInserterWithPerturbation objects consistently across MPI processes
        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        auto comm_world = hemelb::net::MpiCommunicator::World();
        comm_world.Broadcast(seed, 0);
        std::stringstream message;
        message << "RBC insertion random seed: " << std::hex << std::showbase << seed;
        hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>(message.str());

        // Now gets data for cell insertion

        // There are potentially multiple cell inserters each with their own
        // insertion criteria
        CompositeRBCInserter composite;
        for (auto insNode = node.GetChildOrThrow("insertcell"); insNode != insNode.Missing();
            insNode = insNode.NextSiblingOrNull("insertcell"))
        {

          auto const offset = GetNonDimensionalValue<LatticeTime>(insNode,
                                                                  "offset",
                                                                  "s",
                                                                  converter,
                                                                  0);
          auto const templateName = insNode.GetAttributeOrThrow("template");
          if (templateCells.count(templateName) == 0)
          {
            throw Exception() << "Template cell name does not match a known template cell";
          }
          auto cell = templateCells.find(templateName)->second->clone();
          auto const flowExtension = readFlowExtension(node, converter);

          // Rotate cell to align z axis with given position, and then z axis with flow
          // If phi == 0, then cell symmetry axis is aligned with the flow
          auto const theta = GetNonDimensionalValue<Angle>(insNode, "theta", "rad", converter, 0e0);
          auto const phi = GetNonDimensionalValue<Angle>(insNode, "phi", "rad", converter, 0e0);
          auto const x = GetNonDimensionalValue<LatticeDistance>(insNode, "x", "m", converter, 0e0);
          auto const y = GetNonDimensionalValue<LatticeDistance>(insNode, "y", "m", converter, 0e0);
          auto const z = GetNonDimensionalValue<LatticeDistance>(insNode, "z", "m", converter, 0e0);

          util::Matrix3D rotateToFlow, rotation;
          rotateTranslateCellToFlow(cell,
                                    theta,
                                    phi,
                                    flowExtension,
                                    LatticePosition(x, y, z),
                                    rotateToFlow,
                                    rotation);

          // Drops first cell when time reaches offset, and then every deltaTime thereafter.
          // Note: c++14 will allow more complex captures. Until then, we will need to create
          // semi-local lambda variables on the stack as shared pointers. Where semi-local means the
          // variables should live as long as the lambda. But longuer than a single call.
          auto const timeStep = GetNonDimensionalValue<LatticeTime>(insNode,
                                                                    "every",
                                                                    "s",
                                                                    converter);
          auto const dt = GetNonDimensionalValue<LatticeTime>(insNode,
                                                              "delta_t",
                                                              "s",
                                                              converter,
                                                              0e0);
          auto time = std::make_shared<LatticeTime>(timeStep - 1e0
              + std::numeric_limits<LatticeTime>::epsilon() - offset);

          std::default_random_engine randomGenerator(seed);
          std::uniform_real_distribution<double> uniformDistribution(-1.0, 1.0);

          auto condition = [time, timeStep, dt, uniformDistribution, randomGenerator]() mutable
          {
            *time += 1e0;
            if(*time >= timeStep)
            {
              *time -= timeStep + dt * uniformDistribution(randomGenerator);
              return true;
            }
            return false;
          };
          auto const dtheta = GetNonDimensionalValue<Angle>(insNode,
                                                            "delta_theta",
                                                            "rad",
                                                            converter,
                                                            0e0);
          auto const dphi = GetNonDimensionalValue<Angle>(insNode,
                                                          "delta_phi",
                                                          "rad",
                                                          converter,
                                                          0e0);
          auto const dx = GetNonDimensionalValue<LatticeDistance>(insNode,
                                                                  "delta_x",
                                                                  "m",
                                                                  converter,
                                                                  0e0);
          auto const dy = GetNonDimensionalValue<LatticeDistance>(insNode,
                                                                  "delta_y",
                                                                  "m",
                                                                  converter,
                                                                  0e0);

          composite.AddInserter(std::static_pointer_cast<RBCInserter>(std::make_shared<
              RBCInserterWithPerturbation>(condition,
                                           std::move(cell),
                                           rotation,
                                           dtheta,
                                           dphi,
                                           rotateToFlow * LatticePosition(dx, 0, 0),
                                           rotateToFlow * LatticePosition(0, dy, 0),
                                           seed)));
          seed++;
        }

        return composite;
      }
    }

    Cell::Moduli readModuli(io::xml::Element const& node, util::UnitConverter const &converter)
    {
      redblood::Cell::Moduli moduli;
      auto const moduliNode = node.GetChildOrNull("moduli");
      moduli.bending = GetNonDimensionalValue(moduliNode, "bending", "Nm", converter, 2e-19);
      moduli.surface = GetNonDimensionalValue(moduliNode, "surface", "lattice", converter, 1e0);
      moduli.volume = GetNonDimensionalValue(moduliNode, "volume", "lattice", converter, 1e0);
      moduli.dilation = GetNonDimensionalValue(moduliNode, "dilation", "lattice", converter, 0.75);
      if (1e0 < moduli.dilation or moduli.dilation < 0.5)
      {
        log::Logger::Log<log::Critical, log::Singleton>("Dilation moduli is outside the recommended range 1e0 >= m >= 0.5");
      }
      moduli.strain = GetNonDimensionalValue(moduliNode, "strain", "N/m", converter, 5e-6);
      return moduli;
    }

    std::shared_ptr<Node2NodeForce> readNode2NodeForce(io::xml::Element const& node,
						       util::UnitConverter const & converter)
    {
      auto result = std::make_shared<Node2NodeForce>(1e0 / converter.ConvertToLatticeUnits("Nm", 1e0), 1, 2);
      if (node == node.Missing())
      {
        return result;
      }
      result->intensity = GetNonDimensionalValue(node,
                                                "intensity",
                                                "Nm",
                                                converter,
                                                result->intensity);
      result->cutoff = GetNonDimensionalValue(node,
                                             "cutoffdistance",
                                             "lattice",
                                             converter,
                                             result->cutoff);
      if (2e0 * result->cutoff > Dimensionless(Traits<>::Stencil::GetRange()))
      {
        log::Logger::Log<log::Warning, log::Singleton>("Input inconsistency: cell-cell and cell-wall interactions larger then stencil size\n"
                                                       "See issue #586.");
        throw Exception() << "Cell-cell interaction longuer that stencil size permits";
      }
      auto const exponentNode = node != node.Missing() ?
        node.GetChildOrNull("exponent") :
        node.Missing();
      if (exponentNode != node.Missing())
      {
        exponentNode.GetAttributeOrThrow("value", result->exponent);
      }
      return result;
    }

    //! Reads multiple templates from XML and stores in container
    std::unique_ptr<TemplateCellContainer> readTemplateCells(io::xml::Element const& topNode,
							     const configuration::SimConfig& fullconfig,
                                                             util::UnitConverter const& converter)
    {
      std::unique_ptr<TemplateCellContainer> result(new TemplateCellContainer);
      // read flow extensions, if they exist
      std::shared_ptr<std::vector<FlowExtension>> flowExtensions(readFlowExtensions(topNode,
                                                                                    converter).release());
      // Then read template cells
      auto const redbloodcellsNode = topNode.GetChildOrNull("redbloodcells");
      if (not redbloodcellsNode)
      {
        return result;
      }
      auto const cellsNode = redbloodcellsNode.GetChildOrNull("cells");
      if (not cellsNode)
      {
        return result;
      }
      auto cellNode = cellsNode.GetChildOrNull("cell");
      for (; cellNode; cellNode = cellNode.NextSiblingOrNull("cell"))
      {
        auto const name = cellNode.GetAttributeOrNull("name");
        auto const key = name != nullptr ?
          *name :
          "default";
        if (result->count(key) != 0)
        {
          throw Exception() << "Multiple template mesh with same name";
        }
        auto cell = readCell(cellNode, fullconfig, converter);
        if (!validateCellEdgeLengths(*cell))
        {
          log::Logger::Log<log::Critical, log::Singleton>("Average edge length in cell mesh not in [0.7, 1.3]");
        }
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

    bool validateCellEdgeLengths(const CellBase& cell)
    {
      auto edgeLength = cell.GetAverageEdgeLength();

      // Acceptable average edge length to voxel size ratio is [0.7, 1.3].
      // Note that cell vertices location is given in lattice units, therefore
      // no need to normalise again.
      return (edgeLength >= 0.7 && edgeLength <= 1.3);
    }

    std::unique_ptr<MeshIO> make_meshio(io::xml::Element const& shape, const std::string& format_attr_name) {
      auto const format = shape.GetAttributeOrThrow(format_attr_name);
      if (format == "VTK") {
	return std::make_unique<VTKMeshIO>();
      } else if (format == "Krueger") {
	log::Logger::Log<log::Warning, log::Singleton>("Krueger format meshes are deprecated, move to VTK when you can.");
	return std::make_unique<KruegerMeshIO>();
      } else {
	throw Exception() << "Invalid " << format_attr_name << " '" << format << "' on element " << shape.GetPath();
      }
    }

    std::unique_ptr<CellBase> readCell(io::xml::Element const& node,
				       const configuration::SimConfig& fullconfig,
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
      auto const shape = cellNode.GetChildOrThrow("shape");
      std::string const mesh_path = fullconfig.RelPathToFullPath(
          shape.GetAttributeOrThrow("mesh_path")
	);
      auto io = make_meshio(shape, "mesh_format");
      auto const mesh_data = io->readFile(mesh_path, true);
      auto const scale = GetNonDimensionalValue<LatticeDistance>(cellNode, "scale", "m", converter);
      std::string const * reference_mesh_path =
          shape.GetAttributeOrNull("reference_mesh_path");
      std::unique_ptr<Cell> cell;
      if (reference_mesh_path)
      {
	auto io = make_meshio(shape, "reference_mesh_format");
        auto const reference_mesh_data = io->readFile(*reference_mesh_path, true);
        assert(mesh_data->facets == reference_mesh_data->facets);
        assert(volume(mesh_data->vertices, reference_mesh_data->facets) > 0.0);
        cell = std::unique_ptr<Cell>(new Cell(mesh_data->vertices, reference_mesh_data, scale, name));
      }
      else
      {
        cell = std::unique_ptr<Cell>(new Cell(mesh_data->vertices, Mesh(mesh_data), scale, name));
      }
      *cell *= scale;
      cell->moduli = readModuli(cellNode, converter);

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
        io::xml::Element const& inletsNode, util::UnitConverter const& converter,
        TemplateCellContainer const& templateCells)
    {
      // Check if we have any cell insertion action going on at all
      auto inlet = findFirstInletWithCellInsertion(inletsNode.GetChildOrThrow("inlet"));
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
      readFlowExtensions(outletsNode, converter, *result);

      auto inletsNode = topNode.GetChildOrThrow("inlets");
      readFlowExtensionsWithoutInsertElement(inletsNode, converter, *result);

      // Then transforms them to cell outlets: should start somewhere near the end of fadelength
      for (auto &flowExt : *result)
      {
        auto length = flowExt.length - flowExt.fadeLength;
        flowExt.origin += flowExt.normal * (flowExt.length - length);
        flowExt.length = length;
        flowExt.fadeLength = length;
      }
      return result;
    }
  }
}
