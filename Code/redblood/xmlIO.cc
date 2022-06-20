// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <utility>
#include <chrono>
#include "configuration/SimConfig.h"
#include "io/xml.h"
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
	return parent.and_then([&](io::xml::Element const& el) {
	    return el.GetChildOrNull(elemname);
	  }).transform([&](io::xml::Element const& el) {
	      return GetNonDimensionalValue<T>(el, units);
	    }).value_or(std::make_pair(default_, units));
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

      std::vector<FlowExtension> readFlowExtensionsOfOneType(io::xml::Element const& ioletsNode,
                              util::UnitConverter const& converter,
                              bool mustHaveFlowExtension = false)
      {
	std::vector<FlowExtension> results;
        if (ioletsNode) {
	  auto const name = ioletsNode.GetName().substr(0, ioletsNode.GetName().size() - 1);
	  for (
	     auto ioletNode = ioletsNode.GetChildOrNull(name);
	     ioletNode != ioletNode.Missing();
	     ioletNode = ioletNode.NextSiblingOrNull(name)
	     ) {
	    if (ioletNode.GetChildOrNull("flowextension"))
            {
	      results.push_back(readFlowExtension(ioletNode, converter));
            }
            else if (mustHaveFlowExtension)
            {
              throw Exception() << "Could not find flow extension in iolet";
            }
          }
	}
	return results;
      }

      // Read physical outlets that are somehow configured as numerical inlets in the XML parameter file
      std::vector<FlowExtension> readFlowExtensionsWithoutInsertElement(io::xml::Element const& ioletsNode,
                              util::UnitConverter const& converter,
                              bool mustHaveFlowExtension = false)
      {
	std::vector<FlowExtension> results;
        if (ioletsNode)
        {
	  auto const name = ioletsNode.GetName().substr(0, ioletsNode.GetName().size() - 1);
	  for (
	    auto ioletNode = ioletsNode.GetChildOrNull(name);
	    ioletNode != ioletNode.Missing();
	    ioletNode = ioletNode.NextSiblingOrNull(name)
	       ) {
	    if (ioletNode.GetChildOrNull("flowextension"))
	    {
              if (ioletNode.GetChildOrNull("insertcell"))
	      {
		continue;
	      }
	      else
	      {
		results.push_back(readFlowExtension(ioletNode, converter));
	      }
	    }
	    else if (mustHaveFlowExtension)
	    {
	      throw Exception() << "Could not find flow extension in iolet";
	    }
	  }
	}
	return results;
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
      auto def_intensity = 1e0 / converter.ConvertToLatticeUnits("Nm", 1e0);
      auto intensity = GetNonDimensionalValue(node, "intensity", "Nm", converter, def_intensity);

      auto cutoffdist = GetNonDimensionalValue(node, "cutoffdistance", "lattice", converter, 1.0);

      // exponent doesnt have units apparently
      auto exponent = node.and_then(
        [](io::xml::Element const& el) {return el.GetChildOrNull("exponent");}
      ).transform(
	[](io::xml::Element const& el) { return el.GetAttributeOrThrow<size_t>("value"); }
      ).value_or(2);

      if (2e0 * cutoffdist > Dimensionless(Traits<>::Stencil::GetRange()))
      {
        log::Logger::Log<log::Warning, log::Singleton>("Input inconsistency: cell-cell and cell-wall interactions larger then stencil size\n"
                                                       "See issue #586.");
        throw Exception() << "Cell-cell interaction longuer that stencil size permits";
      }
      return std::make_shared<Node2NodeForce>(intensity, cutoffdist, exponent);
    }

    //! Reads multiple templates from XML and stores in container
    std::unique_ptr<TemplateCellContainer> readTemplateCells(io::xml::Element const& topNode,
							     const configuration::SimConfig& fullconfig,
                                                             util::UnitConverter const& converter)
    {
      auto result = std::make_unique<TemplateCellContainer>();

      // Read flow extensions, if they exist.
      // Needs to be a shared pointer so FaderCell can share ownership.
      auto flowExtensions = [&]() {
	using VFE = std::vector<FlowExtension>;
	std::shared_ptr<VFE> ans;
	VFE tmp = readFlowExtensions(topNode, converter);
	if (tmp.size()) {
	  ans = std::make_shared<VFE>(std::move(tmp));
	}
	return ans;
      }();

      // Then read template cells
      auto cellNode = topNode.GetChildOrNull("redbloodcells").and_then([](auto&& _) { return _.GetChildOrNull("cells"); }).and_then([](auto&&_) { return _.GetChildOrNull("cell"); });
      for (; cellNode; cellNode = cellNode.NextSiblingOrNull("cell"))
      {
        auto const key = cellNode.GetAttributeMaybe("name").value_or("default");
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
      return result;
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
      auto const name = cellNode.GetAttributeMaybe("name").value_or("default");
      auto const shape = cellNode.GetChildOrThrow("shape");
      std::string const mesh_path = fullconfig.RelPathToFullPath(
          shape.GetAttributeOrThrow("mesh_path")
	);
      auto io = make_meshio(shape, "mesh_format");
      auto const mesh_data = io->readFile(mesh_path, true);
      auto const scale = GetNonDimensionalValue<LatticeDistance>(cellNode, "scale", "m", converter);
      std::optional<std::string> reference_mesh_path = shape.GetAttributeMaybe("reference_mesh_path");
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

    std::vector<FlowExtension> readFlowExtensions(
        io::xml::Element const& topNode, util::UnitConverter const& converter)
    {
      using VFE = std::vector<FlowExtension>;
      return topNode.transform([&](io::xml::Element const& topNode) {
	  auto do_iolets = [&](io::xml::Element const& node) {
	    return readFlowExtensionsOfOneType(node, converter);
	  };
	  auto ins = topNode.GetChildOrNull("inlets").transform(do_iolets).value_or(VFE{});
	  auto outs = topNode.GetChildOrNull("outlets").transform(do_iolets).value_or(VFE{});
	  ins.insert(ins.end(), outs.begin(), outs.end());
	  return ins;
	}).value_or<VFE>({});
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
        results.emplace_back(readSingleRBCInserter(inlet, converter, templateCells));
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
      auto outs = readFlowExtensionsOfOneType(outletsNode, converter);

      auto inletsNode = topNode.GetChildOrThrow("inlets");
      auto ins = readFlowExtensionsWithoutInsertElement(inletsNode, converter);

      // Then transforms them to cell outlets: should start somewhere near the end of fadelength
      auto trans = [](FlowExtension flowExt) {
        auto length = flowExt.length - flowExt.fadeLength;
        flowExt.origin += flowExt.normal * (flowExt.length - length);
        flowExt.length = length;
        flowExt.fadeLength = length;
	return flowExt;
      };
      std::transform(outs.begin(), outs.end(), std::back_inserter(*result), trans);
      std::transform(ins.begin(), ins.end(), std::back_inserter(*result), trans);
      return result;
    }
  }
}
