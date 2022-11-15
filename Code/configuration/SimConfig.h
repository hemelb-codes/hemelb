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

namespace hemelb::redblood {
    class RBCConfig;
}

namespace hemelb::configuration
{
    template<typename T>
    void GetDimensionalValue(const io::xml::Element& elem, const std::string& units, T& value)
    {
      const std::string& got = elem.GetAttributeOrThrow("units");
      if (got != units)
      {
        throw Exception() << "Invalid units for element " << elem.GetPath() << ". Expected '"
            << units << "', got '" << got << "'";
      }

      elem.GetAttributeOrThrow("value", value);
    }
    template<typename T>
    T GetDimensionalValue(const io::xml::Element& elem, const std::string& units)
    {
      const std::string& got = elem.GetAttributeOrThrow("units");
      if (got != units)
      {
        throw Exception() << "Invalid units for element " << elem.GetPath() << ". Expected '"
            << units << "', got '" << got << "'";
      }

      return elem.GetAttributeOrThrow<T>("value");
    }

    template<typename T>
    void GetDimensionalValueInLatticeUnits(const io::xml::Element& elem, const std::string& units, const util::UnitConverter& converter, T& value)
    {
      if (units == "lattice")
	return GetDimensionalValue(elem, units, value);

      T phys;
      GetDimensionalValue(elem, units, phys);
      value = converter.ConvertToLatticeUnits(units, phys);
    }

    template<typename T>
    T GetDimensionalValueInLatticeUnits(const io::xml::Element& elem, const std::string& units, const util::UnitConverter& converter)
    {
      auto phys = GetDimensionalValue<T>(elem, units);

      if (units == "lattice")
	return phys;
      else
	return converter.ConvertToLatticeUnits(units, phys);
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
      CheckpointIC(std::optional<LatticeTimeStep> t, const std::string& cp, std::optional<std::string> const & maybeOff);
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

    struct IoletConfigBase {
        PhysicalPosition position;
        util::Vector3D<double> normal;
        std::optional<std::uint64_t> warmup_steps;
        std::optional<FlowExtensionConfig> flow_extension;
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
        virtual ~SimConfig();

        // Turn an input XML-relative path into a full path
        [[nodiscard]] path RelPathToFullPath(const std::string& path) const;

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
          return rbcConf != nullptr;
        }

	inline const redblood::RBCConfig* GetRBCConfig() const {
	  return rbcConf;
	}
      protected:

        /**
         * Check that the iolet is OK for the CMake configuration.
         * @param ioletEl
         * @param requiredBC
         */
        virtual void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                            const std::string& requiredBC) const;

      private:
        void DoIO(const io::xml::Element xmlNode);
        void DoIOForSimulation(const io::xml::Element simEl);
        void DoIOForGeometry(const io::xml::Element geometryEl);

        std::vector<IoletConfig> DoIOForInOutlets(const io::xml::Element xmlNode) const;

        void DoIOForBaseInOutlet(const io::xml::Element& ioletEl, IoletConfigBase& value) const;

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

        path xmlFilePath;
        path dataFilePath;

        std::vector<extraction::PropertyOutputFile> propertyOutputs;
        /**
         * True if the file has a colloids section.
         */
        bool hasColloidSection = false;

        MonitoringConfig monitoringConfig; ///< Configuration of various checks/tests

        // We want to keep the RBC types isolated in their own library
        // so put all the configuration in a type that is opaque to
        // the rest of the code.
        //
        // nullptr => no RBC
        redblood::RBCConfig* rbcConf = nullptr;

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
