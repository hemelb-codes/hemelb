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
    // This will be given as a tuple of string like values. The consumer (i.e. SimBuilder etc) will have
    // to take care around this.

    using UnitUnion = std::vector<std::string_view>;

    inline void check_unit_spec(const io::xml::Element& elem, std::string_view actual, std::string_view const& expected) {
        if (actual != expected)
            throw Exception() << "Invalid units for element " << elem.GetPath()
                              << "."" Expected '" << expected
                              << "', got '" << actual << "'";
    }
    inline void check_unit_spec(const io::xml::Element& elem, std::string_view actual, UnitUnion const& expected) {
        if (expected.empty())
            throw Exception() << "Invalid empty unit spec for element " << elem.GetPath();

        if (std::find(expected.begin(), expected.end(), actual) == expected.end()) {
            auto e = Exception() << "Invalid units for element " << elem.GetPath()
                                 << ". Expected one of ";
            char sep = '(';
            for (auto& exp: expected) {
                e << sep << '"' << exp << '"';
                sep = ',';
            }
            e << "), got \"" << actual << '"';
        }
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

    //! Check the units of the quantity and return the value and the actual units
    template<typename T>
    auto GetDimensionalValue(const io::xml::Element& elem, UnitUnion units)
    {
        auto got = elem.GetAttributeOrThrow("units");
        check_unit_spec(elem, got, units);

        return std::pair<T, std::string>{elem.GetAttributeOrThrow<T>("value"), got};
    }

    //! Given an element (@param elem), check for a child with the given @param name.
    //! If it exists, return the unit-checked (against @param unit) value.
    //! If it doesn't exist, return @param default_value.
    template <typename T>
    T GetDimensionalValueWithDefault(const io::xml::Element& elem,
                                     std::string_view name, std::string_view unit, T default_value) {
        return elem
                .and_then([&](const io::xml::Element& _) { return _.GetChildOrNull(name); })
                .transform([&](const io::xml::Element& _) { return GetDimensionalValue<T>(_, unit); })
                .value_or(default_value);
    }
    //! Given an element (@param elem), check for a child with the given @param name.
    //! If it exists, return the unit-checked (against @param unit) value.
    //! If it doesn't exist, return @param default_value.
    template <typename T>
    std::pair<T, std::string> GetDimensionalValueWithDefault(
            const io::xml::Element& elem, std::string_view name,
            UnitUnion unit, std::pair<T, std::string> default_value_units
    ) {
        return elem
                .and_then([&](const io::xml::Element& _) { return _.GetChildOrNull(name); })
                .transform([&](const io::xml::Element& _) { return GetDimensionalValue<T>(_, unit); })
                .value_or(default_value_units);
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
      CheckpointIC(std::optional<LatticeTimeStep> t, std::string cp, std::optional<std::string> maybeOff);
      std::string cpFile;
      std::optional<std::string> maybeOffFile;
    };

    // Variant including null state
    using ICConfig = std::variant<std::monostate, EquilibriumIC, CheckpointIC>;

    struct TimeInfo {
        std::uint64_t total_steps;
        std::uint64_t warmup_steps;
        PhysicalTime step_s;
    };

    struct SpaceInfo {
        PhysicalDistance step_m;
        PhysicalPosition geometry_origin_m;
    };

    struct FluidInfo {
        PhysicalDensity density_kgm3;
        PhysicalDynamicViscosity viscosity_Pas;
        PhysicalPressure reference_pressure_mmHg;
    };

    struct GlobalSimInfo {
        lb::StressTypes stress_type;
        TimeInfo time;
        SpaceInfo space;
        FluidInfo fluid;
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
        PhysicalTime offset_s;
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
        double intensity;
        std::string intensity_units;
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
      public:
        using path = std::filesystem::path;

        static std::unique_ptr<SimConfig> New(const path& p);

      protected:
    	explicit SimConfig(const path& p);
        void Init();

      public:
        virtual ~SimConfig() = default;

        // Turn an input XML-relative path into a full path
        [[nodiscard]] path RelPathToFullPath(std::string_view path) const;

        void Save(std::string path); // TODO this method should be able to be CONST
        // but because it uses DoIo, which uses one function signature for both reading and writing, it cannot be.

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
        path const& GetColloidConfigPath() const
        {
          return xmlFilePath;
        }
        /**
         * True if the XML file has a section specifying colloids.
         * @return
         */
        bool HasColloidSection() const;

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

        /**
         * Check that the iolet is OK for the CMake configuration.
         * @param ioletEl
         * @param requiredBC
         */
        virtual void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                            const std::string& requiredBC) const;

      public:
        void DoIO(const io::xml::Element xmlNode);
        void DoIOForSimulation(const io::xml::Element simEl);
        void DoIOForGeometry(const io::xml::Element geometryEl);

        std::vector<IoletConfig> DoIOForInOutlets(const io::xml::Element xmlNode) const;

        void DoIOForBaseInOutlet(const io::xml::Element& ioletEl, IoletConfigBase& ioletConf) const;

        IoletConfig DoIOForPressureInOutlet(const io::xml::Element& ioletEl) const;
        IoletConfig DoIOForCosinePressureInOutlet(const io::xml::Element& ioletEl) const;
        IoletConfig DoIOForFilePressureInOutlet(const io::xml::Element& ioletEl) const;
        IoletConfig DoIOForMultiscalePressureInOutlet(
            const io::xml::Element& ioletEl) const;

        IoletConfig DoIOForVelocityInOutlet(const io::xml::Element& ioletEl) const;
        IoletConfig DoIOForParabolicVelocityInOutlet(
            const io::xml::Element& ioletEl) const;
        /**
         * Reads a Womersley velocity iolet definition from the XML config file and returns
         * an InOutLetWomersleyVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetWomersleyVelocity object
         */
        IoletConfig DoIOForWomersleyVelocityInOutlet(const io::xml::Element& ioletEl) const;

        /**
         * Reads a file velocity iolet definition from the XML config file and returns
         * an InOutLetFileVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetFileVelocity object
         */
        IoletConfig DoIOForFileVelocityInOutlet(const io::xml::Element& ioletEl) const;

        void DoIOForProperties(const io::xml::Element& xmlNode);
        void DoIOForProperty(io::xml::Element xmlNode, bool isLoading);
        extraction::OutputField DoIOForPropertyField(const io::xml::Element& xmlNode);
        extraction::PropertyOutputFile DoIOForPropertyOutputFile(
            const io::xml::Element& propertyoutputEl);
        extraction::StraightLineGeometrySelector* DoIOForLineGeometry(
            const io::xml::Element& xmlNode);
        extraction::PlaneGeometrySelector* DoIOForPlaneGeometry(const io::xml::Element&);
        extraction::SurfacePointSelector* DoIOForSurfacePoint(const io::xml::Element&);

        void DoIOForInitialConditions(io::xml::Element parent);
	void DoIOForCheckpointFile(const io::xml::Element& checkpointEl);

        /**
         * Reads monitoring configuration from XML file
         *
         * @param monEl in memory representation of <monitoring> xml element
         */
        void DoIOForMonitoring(const io::xml::Element& monEl);

        /**
         * Reads configuration of steady state flow convergence check from XML file
         *
         * @param convEl in memory representation of the <steady_flow_convergence> XML element
         */
        void DoIOForSteadyFlowConvergence(const io::xml::Element& convEl);

        /**
         * Reads the configuration of one of the possible several converge criteria provided
         *
         * @param criterionEl in memory representation of the <criterion> XML element
         */
        void DoIOForConvergenceCriterion(const io::xml::Element& criterionEl);

        TemplateCellConfig readCell(const io::xml::Element& cellNode) const;
        std::map<std::string, TemplateCellConfig> readTemplateCells(io::xml::Element const& cellsEl) const;
        RBCConfig DoIOForRedBloodCells(const io::xml::Element& rbcEl) const;

    private:
        path xmlFilePath;
        path dataFilePath;

        std::vector<extraction::PropertyOutputFile> propertyOutputs;
        /**
         * True if the file has a colloids section.
         */
        bool hasColloidSection = false;

        MonitoringConfig monitoringConfig; ///< Configuration of various checks/tests

        std::optional<RBCConfig> rbcConf;

      protected:
        GlobalSimInfo sim_info;
        ICConfig initial_condition;
        // These have to contain pointers because there are multiple derived types that might be
        // instantiated.
        std::vector<IoletConfig> inlets;
        std::vector<IoletConfig> outlets;
      private:
    };
}

#endif /* HEMELB_CONFIGURATION_SIMCONFIG_H */
