// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_SIMCONFIG_H
#define HEMELB_CONFIGURATION_SIMCONFIG_H

#include <optional>
#include <variant>
#include <vector>

#include "configuration/MonitoringConfig.h"
#include "util/Vector3D.h"
#include "lb/LbmParameters.h"
#include "lb/iolets/InOutLets.h"
#include "extraction/GeometrySelectors.h"
#include "extraction/PropertyOutputFile.h"
#include "io/xml.h"
#include "quantity.h"

namespace hemelb::configuration
{
    // Note on specifying physical quantities.
    //
    // In the XML, these are encoded in attributes of an element (the name
    // of the element doesn't matter for these purposes). The two required
    // attributes are:
    // 1. "units" (e.g. "m" for a physical length or "lattice" for something specified in scaled units)
    // 2. "value" which is a string encoding of the representation (typically a double)
    //
    // Now some parts of the code (e.g. RBC) can specify units as multiple options (e.g. physical OR lattice)
    // for these, use the quantity union types. The consumer (i.e. SimBuilder etc) will have
    // to take care around this.
    inline void check_unit_spec(const io::xml::Element& elem, std::string_view actual, std::string_view const& expected) {
        if (actual != expected)
            throw Exception() << "Invalid units for element " << elem.GetPath()
                              << "."" Expected '" << expected
                              << "', got '" << actual << "'";
    }

    //! Check the units of the quantity and decode the value into @param value
    template<typename T>
    void GetDimensionalValue(const io::xml::Element& elem, std::string_view units, T& value)
    {
        auto got = elem.GetAttributeOrThrow("units");
        check_unit_spec(elem, got, units);
        elem.GetAttributeOrThrow("value", value);
    }
    //! Check the units of the quantity and return the value
    template<typename T>
    T GetDimensionalValue(const io::xml::Element& elem, std::string_view units)
    {
        auto got = elem.GetAttributeOrThrow("units");
        check_unit_spec(elem, got, units);

        return elem.GetAttributeOrThrow<T>("value");
    }

    //! Given an element (@param elem), check for a child with the given @param name.
    //! If it exists, return the unit-checked (against @param unit) value.
    //! If it doesn't exist, return @param default_value.
    template <typename T>
    T GetDimensionalValueWithDefault(const io::xml::Element& elem,
                                     char const* name, char const* unit, T default_value) {
        using Element = io::xml::Element;
        return elem
                .and_then([&](const Element& _) { return _.GetChildOrNull(name); })
                .transform([&](const Element& _) { return GetDimensionalValue<T>(_, unit); })
                .value_or(default_value);
    }

    template <QuantityUnion QUnion, IsVariantAlternative<QUnion> DefaultQ>
    QUnion GetDimensionalValueWithDefault(const io::xml::Element& elem, char const* name, DefaultQ const& default_q) {
        using RepT = typename quantity_union_traits<QUnion>::representation_type;
        using io::xml::Element;
        return elem
            .and_then([&](Element const& _) { return _.GetChildOrNull(name); })
            .transform([](Element const& _) {
                auto xml_units = _.GetAttributeOrThrow("units");
                auto xml_val = _.GetAttributeOrThrow<RepT>("value");
                return quantity_union_factory<QUnion>()(xml_val, xml_units);
            })
            .value_or(default_q);
    }

    // Base for initial conditions configuration
    struct ICConfigBase {
      ICConfigBase(std::optional<LatticeTimeStep> t);
      std::optional<LatticeTimeStep> t0;
    };

    // Uniform equilibrium IC
    struct EquilibriumIC : ICConfigBase {
      EquilibriumIC(std::optional<LatticeTimeStep> t, PhysicalPressure p);
      EquilibriumIC(std::optional<LatticeTimeStep> t, PhysicalPressure p, const PhysicalVelocity& v);
      PhysicalPressure p_mmHg;
      PhysicalVelocity v_ms;
    };

    // Read from checkpoint IC
    struct CheckpointIC : ICConfigBase {
      CheckpointIC(std::optional<LatticeTimeStep> t, std::filesystem::path cp, std::optional<std::filesystem::path> maybeOff);
      std::filesystem::path cpFile;
      std::optional<std::filesystem::path> maybeOffFile;
    };

    // Variant including null state
    using ICConfig = std::variant<std::monostate, EquilibriumIC, CheckpointIC>;

    struct TimeInfo {
        std::uint64_t total_steps = 0;
        std::uint64_t warmup_steps = 0;
        PhysicalTime step_s = 0.0;
    };

    struct SpaceInfo {
        PhysicalDistance step_m = 0;
        PhysicalPosition geometry_origin_m = {0.0, 0.0, 0.0};
    };

    struct FluidInfo {
        PhysicalDensity density_kgm3 = DEFAULT_FLUID_DENSITY_Kg_per_m3;
        PhysicalDynamicViscosity viscosity_Pas = DEFAULT_FLUID_VISCOSITY_Pas;
        PhysicalPressure reference_pressure_mmHg = 0.0;
    };

    struct CheckpointInfo {
        LatticeTimeStep period;
    };

    struct GlobalSimInfo {
        lb::StressTypes stress_type;
        TimeInfo time;
        SpaceInfo space;
        FluidInfo fluid;
        std::optional<CheckpointInfo> checkpoint;
    };

    struct FlowExtensionConfig {
        PhysicalDistance length_m;
        PhysicalDistance radius_m;
        PhysicalDistance fadelength_m;
        util::Vector3D<double> normal;
        PhysicalPosition origin_m;
    };

    struct CellInserterConfig {
        std::int64_t seed;
        std::string template_name;
        quantity_union<double, "s", "lattice"> offset;
        Angle theta_rad;
        Angle phi_rad;
        PhysicalDisplacement translation_m;
        PhysicalTime drop_period_s;
        PhysicalTime dt_s;
        Angle dtheta_rad;
        Angle dphi_rad;
        PhysicalDistance dx_m;
        PhysicalDistance dy_m;
    };

    struct IoletConfigBase {
        PhysicalPosition position;
        util::Vector3D<double> normal;
        std::optional<std::uint64_t> warmup_steps;
        std::optional<FlowExtensionConfig> flow_extension;
        std::vector<CellInserterConfig> cell_inserters;
    };

    // Consider removing
    struct PressureIoletConfig : IoletConfigBase {
    };

    struct CosinePressureIoletConfig : PressureIoletConfig {
        // Note pressure difference
        PhysicalPressure amp_mmHg;
        // Note absolute pressure
        PhysicalPressure mean_mmHg;
        Angle phase_rad;
        PhysicalTime period_s;
    };

    struct FilePressureIoletConfig : PressureIoletConfig {
        std::filesystem::path file_path;
    };

    struct MultiscalePressureIoletConfig : PressureIoletConfig {
        PhysicalPressure pressure_reference_mmHg;
        PhysicalVelocity velocity_reference_ms;
        std::string label;
    };

    struct VelocityIoletConfig : IoletConfigBase {
    };

    struct ParabolicVelocityIoletConfig : VelocityIoletConfig {
        PhysicalDistance radius_m;
        PhysicalSpeed max_speed_ms;
    };

    struct WomersleyVelocityIoletConfig : VelocityIoletConfig {
        PhysicalDistance radius_m;
        PhysicalPressureGradient pgrad_amp_mmHgm;
        PhysicalTime period_s;
        double womersley;
    };

    struct FileVelocityIoletConfig : VelocityIoletConfig {
        std::filesystem::path file_path;
        PhysicalDistance radius_m;
    };

    using IoletConfig = std::variant<std::monostate,
            CosinePressureIoletConfig, FilePressureIoletConfig, MultiscalePressureIoletConfig,
            ParabolicVelocityIoletConfig, WomersleyVelocityIoletConfig, FileVelocityIoletConfig
    >;

    struct VTKMeshFormat {};
    struct KruegerMeshFormat {};
    using MeshFormat = std::variant<std::monostate, VTKMeshFormat, KruegerMeshFormat>;

    struct CellModuli {
        //! Bending energy parameter
        PhysicalEnergy bending_Nm;
        //! Surface energy parameter
        LatticeModulus surface_lat;
        //! Volume energy parameter
        LatticeModulus volume_lat;
        //! Skalak dilation modulus
        LatticeModulus dilation_lat;
        //! Skalak strain modulus
        PhysicalModulus strain_Npm;
    };

    struct TemplateCellConfig {
        std::string name;
        std::filesystem::path mesh_path;
        MeshFormat format;
        PhysicalDistance scale_m;
        std::optional<std::filesystem::path> reference_mesh_path;
        MeshFormat reference_mesh_format;
        CellModuli moduli;
    };

    struct NodeForceConfig {
        quantity_union<double, "Nm", "lattice"> intensity;
        LatticeDistance cutoffdist;
        std::size_t exponent;
    };

    struct RBCConfig {
        LatticeDistance boxSize;
        std::map<std::string, TemplateCellConfig> meshes;
        NodeForceConfig cell2cell;
        NodeForceConfig cell2wall;
        LatticeTimeStep output_period;
    };

    class SimConfig
    {
        friend class SimBuilder;
        friend class SimConfigReader;

    public:
        using path = std::filesystem::path;
        static SimConfig New(const path& p);

        inline GlobalSimInfo const& GetSimInfo() const {
            return sim_info;
        }

        const std::vector<IoletConfig> & GetInlets() const
        {
          return inlets;
        }
        const std::vector<IoletConfig> & GetOutlets() const
        {
          return outlets;
        }
        lb::StressTypes GetStressType() const
        {
          return sim_info.stress_type;
        }
        const path& GetDataFilePath() const
        {
          return dataFilePath;
        }
        LatticeTimeStep GetTotalTimeSteps() const
        {
          return sim_info.time.total_steps;
        }
        LatticeTimeStep GetWarmUpSteps() const
        {
          return sim_info.time.warmup_steps;
        }
        PhysicalTime GetTimeStepLength() const
        {
          return sim_info.time.step_s;
        }
        PhysicalDistance GetVoxelSize() const
        {
          return sim_info.space.step_m;
        }
        PhysicalPosition GetGeometryOrigin() const
        {
          return sim_info.space.geometry_origin_m;
        }
        unsigned int PropertyOutputCount() const
        {
          return propertyOutputs.size();
        }
        extraction::PropertyOutputFile& GetPropertyOutput(unsigned int index)
        {
          return propertyOutputs[index];
        }
        std::vector<extraction::PropertyOutputFile> const& GetPropertyOutputs() const
        {
          return propertyOutputs;
        }

        inline path const& GetColloidConfigPath() const
        {
          return colloid_xml_path.value();
        }

        inline bool HasColloidSection() const {
            return colloid_xml_path.has_value();
        }

        // Get the initial condtion config
        inline const ICConfig& GetInitialCondition() const {
	  return initial_condition;
	}

        /**
         * Return the configuration of various checks/test
         * @return monitoring configuration
         */
        const MonitoringConfig& GetMonitoringConfiguration() const;

        /**
         * True if the XML file has a section specifying red blood cells.
         * @return
         */
        inline bool HasRBCSection() const
        {
          return rbcConf.has_value();
        }

        inline const RBCConfig& GetRBCConfig() const {
            return *rbcConf;
        }

    protected:
        std::shared_ptr<io::xml::Document> source_xml;
        path dataFilePath;

        std::vector<extraction::PropertyOutputFile> propertyOutputs;

        // Path of the colloid XML if colloids section present
        std::optional<path> colloid_xml_path;

        MonitoringConfig monitoringConfig; ///< Configuration of various checks/tests

        std::optional<RBCConfig> rbcConf;

        GlobalSimInfo sim_info;
        ICConfig initial_condition;
        std::vector<IoletConfig> inlets;
        std::vector<IoletConfig> outlets;
    };

    class SimConfigReader {
        using path = std::filesystem::path;
        using Element = io::xml::Element;
        path xmlFilePath;

    public:
        SimConfigReader(path);

        virtual ~SimConfigReader() = default;

        [[nodiscard]] virtual SimConfig Read() const;
        // Turn an input XML-relative path into a full path
        [[nodiscard]] path RelPathToFullPath(std::string_view path) const;

        /**
         * Check that the iolet is OK for the CMake configuration.
         * @param ioletEl
         * @param requiredBC
         */
        virtual void CheckIoletMatchesCMake(const Element& ioletEl,
                                            std::string_view requiredBC) const;
        [[nodiscard]] virtual SimConfig DoIO(const Element xmlNode) const;
        [[nodiscard]] virtual GlobalSimInfo DoIOForSimulation(const Element simEl) const;
        [[nodiscard]] virtual path DoIOForGeometry(const Element geometryEl) const;

        [[nodiscard]] virtual std::vector<IoletConfig> DoIOForInOutlets(GlobalSimInfo const& sim_info, const Element xmlNode) const;

        void DoIOForBaseInOutlet(GlobalSimInfo const& sim_info, const Element& ioletEl, IoletConfigBase& ioletConf) const;

        [[nodiscard]] IoletConfig DoIOForPressureInOutlet(const Element& ioletEl) const;
        [[nodiscard]] IoletConfig DoIOForCosinePressureInOutlet(const Element& ioletEl) const;
        [[nodiscard]] IoletConfig DoIOForFilePressureInOutlet(const Element& ioletEl) const;
        [[nodiscard]] IoletConfig DoIOForMultiscalePressureInOutlet(
                const Element& ioletEl) const;

        [[nodiscard]] IoletConfig DoIOForVelocityInOutlet(const Element& ioletEl) const;
        [[nodiscard]] IoletConfig DoIOForParabolicVelocityInOutlet(
                const Element& ioletEl) const;
        /**
         * Reads a Womersley velocity iolet definition from the XML config file and returns
         * an InOutLetWomersleyVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetWomersleyVelocity object
         */
        [[nodiscard]] IoletConfig DoIOForWomersleyVelocityInOutlet(const Element& ioletEl) const;

        /**
         * Reads a file velocity iolet definition from the XML config file and returns
         * an InOutLetFileVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetFileVelocity object
         */
        [[nodiscard]] IoletConfig DoIOForFileVelocityInOutlet(const Element& ioletEl) const;

        [[nodiscard]] virtual std::vector<extraction::PropertyOutputFile> DoIOForProperties(GlobalSimInfo const& sim_info, const Element& xmlNode) const;
        [[nodiscard]] virtual extraction::OutputField DoIOForPropertyField(GlobalSimInfo const& sim_info, const Element& xmlNode) const;
        [[nodiscard]] virtual extraction::PropertyOutputFile DoIOForPropertyOutputFile(GlobalSimInfo const& sim_info, const Element& propertyoutputEl) const;
        [[nodiscard]] extraction::StraightLineGeometrySelector* DoIOForLineGeometry(
                const Element& xmlNode) const;
        [[nodiscard]] extraction::PlaneGeometrySelector* DoIOForPlaneGeometry(const Element&) const;
        [[nodiscard]] extraction::SurfacePointSelector* DoIOForSurfacePoint(const Element&) const;

        [[nodiscard]] virtual ICConfig DoIOForInitialConditions(Element parent) const;
        //virtual void DoIOForCheckpointFile(const Element& checkpointEl) const;

        /**
         * Reads monitoring configuration from XML file
         *
         * @param monEl in memory representation of <monitoring> xml element
         */
        [[nodiscard]] virtual MonitoringConfig DoIOForMonitoring(const Element& monEl) const;

        /**
         * Reads configuration of steady state flow convergence check from XML file
         *
         * @param convEl in memory representation of the <steady_flow_convergence> XML element
         */
        void DoIOForSteadyFlowConvergence(const Element& convEl, MonitoringConfig& monitoringConfig) const;

        /**
         * Reads the configuration of one of the possible several converge criteria provided
         *
         * @param criterionEl in memory representation of the <criterion> XML element
         */
        void DoIOForConvergenceCriterion(const Element& criterionEl, MonitoringConfig& monitoringConfig) const;

        [[nodiscard]] virtual TemplateCellConfig readCell(const Element& cellNode) const;
        [[nodiscard]] virtual std::map<std::string, TemplateCellConfig> readTemplateCells(Element const& cellsEl) const;
        [[nodiscard]] virtual RBCConfig DoIOForRedBloodCells(SimConfig const&, const Element& rbcEl) const;

    };

}

#endif /* HEMELB_CONFIGURATION_SIMCONFIG_H */
