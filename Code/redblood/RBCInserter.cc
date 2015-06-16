#include "RBCInserter.h"

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


    Cell::Moduli read_moduli(io::xml::Element const& node, util::UnitConverter const &converter)
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

    Node2NodeForce read_node2node_force(
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

    std::shared_ptr<CellBase> read_cell(
        io::xml::Element const& node, util::UnitConverter const& converter)
    {
      if(node  == node.Missing())
      {
        throw Exception() << "Expected non-empty XML node";
      }
      const io::xml::Element & cellNode = node.GetChildOrThrow("cell");
      std::string const mesh_path
        = cellNode.GetChildOrThrow("shape").GetAttributeOrThrow("mesh_path");
      auto const mesh_data = read_mesh(mesh_path);
      auto const scale = GetDimensionalValue<LatticeDistance>(cellNode, "scale", "m", converter);
      auto result = std::make_shared<Cell>(mesh_data->vertices, Mesh(mesh_data), scale);
      *result *= scale;
      result->moduli = read_moduli(cellNode, converter);
      result->nodeWall = read_node2node_force(cellNode, converter);
      return result;
    }
  }
}
