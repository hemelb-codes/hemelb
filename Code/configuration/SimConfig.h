// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_SIMCONFIG_H
#define HEMELB_CONFIGURATION_SIMCONFIG_H

#include <vector>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "configuration/MonitoringConfig.h"
#include "util/Vector3D.h"
#include "lb/LbmParameters.h"
#include "lb/iolets/InOutLets.h"
#include "extraction/GeometrySelectors.h"
#include "extraction/PropertyOutputFile.h"
#include "io/xml/XmlAbstractionLayer.h"

namespace hemelb
{
  namespace redblood {
    class RBCConfig;
  }

  namespace configuration
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
    void GetDimensionalValueInLatticeUnits(const io::xml::Element& elem, const std::string& units, const util::UnitConverter& converter, T& value)
    {
      if (units == "lattice")
	return GetDimensionalValue(elem, units, value);

      T phys;
      GetDimensionalValue(elem, units, phys);
      value = converter.ConvertToLatticeUnits(units, phys);
    }

    // Base for initial conditions configuration
    struct ICConfigBase {
      ICConfigBase(const util::UnitConverter* units, boost::optional<LatticeTimeStep> t);

      const util::UnitConverter* unitConverter;
      boost::optional<LatticeTimeStep> t0;
    };

    // Uniform equilibrium IC
    struct EquilibriumIC : ICConfigBase {
      EquilibriumIC(const util::UnitConverter* units, boost::optional<LatticeTimeStep> t, PhysicalPressure p);
      EquilibriumIC(const util::UnitConverter* units, boost::optional<LatticeTimeStep> t, PhysicalPressure p, const PhysicalVelocity& v);
      PhysicalPressure p_mmHg;
      PhysicalVelocity v_ms;
    };

    // Read from checkpoint IC
    struct CheckpointIC : ICConfigBase {
      CheckpointIC(const util::UnitConverter* units, boost::optional<LatticeTimeStep> t, const std::string& cp);
      std::string cpFile;
    };

    // Variant including null state
    using ICConfig = boost::variant<std::nullptr_t, EquilibriumIC, CheckpointIC>;

    class SimConfig
    {
      public:
	static SimConfig* New(const std::string& path);

      protected:
	SimConfig(const std::string& path);
        void Init();

      public:
        virtual ~SimConfig();

        // Turn an input XML-relative path into a full path
        std::string RelPathToFullPath(const std::string& path) const;

        void Save(std::string path); // TODO this method should be able to be CONST
        // but because it uses DoIo, which uses one function signature for both reading and writing, it cannot be.

        const util::Vector3D<float> & GetVisualisationCentre() const
        {
          return visualisationCentre;
        }
        const std::vector<lb::iolets::InOutLet*> & GetInlets() const
        {
          return inlets;
        }
        const std::vector<lb::iolets::InOutLet*> & GetOutlets() const
        {
          return outlets;
        }
        lb::StressTypes GetStressType() const
        {
          return stressType;
        }
        float GetVisualisationLongitude() const
        {
          return visualisationLongitude;
        }
        float GetVisualisationLatitude() const
        {
          return visualisationLatitude;
        }
        float GetVisualisationZoom() const
        {
          return visualisationZoom;
        }
        float GetVisualisationBrightness() const
        {
          return visualisationBrightness;
        }
        float GetMaximumVelocity() const
        {
          return maxVelocity;
        }
        float GetMaximumStress() const
        {
          return maxStress;
        }
        const std::string& GetDataFilePath() const
        {
          return dataFilePath;
        }
        LatticeTimeStep GetTotalTimeSteps() const
        {
          return totalTimeSteps;
        }
        LatticeTimeStep GetWarmUpSteps() const
        {
          return warmUpSteps;
        }
        PhysicalTime GetTimeStepLength() const
        {
          return timeStepSeconds;
        }
        PhysicalDistance GetVoxelSize() const
        {
          return voxelSizeMetres;
        }
        PhysicalPosition GetGeometryOrigin() const
        {
          return geometryOriginMetres;
        }
        unsigned int PropertyOutputCount() const
        {
          return propertyOutputs.size();
        }
        extraction::PropertyOutputFile * GetPropertyOutput(unsigned int index) const
        {
          return propertyOutputs[index];
        }
        std::vector<extraction::PropertyOutputFile*> const GetPropertyOutputs() const
        {
          return propertyOutputs;
        }
        const std::string GetColloidConfigPath() const
        {
          return colloidConfigPath;
        }
        /**
         * True if the XML file has a section specifying colloids.
         * @return
         */
        bool HasColloidSection() const;

        // Get the initial condtion config
        inline const ICConfig& GetInitialCondition() const {
	  return icConfig;
	}

        const util::UnitConverter& GetUnitConverter() const;

        /**
         * Return the configuration of various checks/test
         * @return monitoring configuration
         */
        const MonitoringConfig* GetMonitoringConfiguration() const;

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
         * Create the unit converter - virtual so that mocks can override it.
         */
        virtual void CreateUnitConverter();

        /**
         * Check that the iolet is OK for the CMake configuration.
         * @param ioletEl
         * @param requiredBC
         */
        virtual void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                            const std::string& requiredBC);

        template<typename T>
        void GetDimensionalValueInLatticeUnits(const io::xml::Element& elem,
                                               const std::string& units, T& value)
        {
          GetDimensionalValue(elem, units, value);
          value = unitConverter->ConvertToLatticeUnits(units, value);
        }

        template<typename T>
        T GetDimensionalValueInLatticeUnits(const io::xml::Element& elem, const std::string& units)
        {
          T ans;
          GetDimensionalValueInLatticeUnits(elem, units, ans);
          return ans;
        }

      private:
        void DoIO(const io::xml::Element xmlNode);
        void DoIOForSimulation(const io::xml::Element simEl);
        void DoIOForGeometry(const io::xml::Element geometryEl);

        std::vector<lb::iolets::InOutLet*> DoIOForInOutlets(const io::xml::Element xmlNode);
        void DoIOForFlowExtension(lb::iolets::InOutLet *, const io::xml::Element &);

        void DoIOForBaseInOutlet(const io::xml::Element& ioletEl, lb::iolets::InOutLet* value);

        lb::iolets::InOutLet* DoIOForPressureInOutlet(const io::xml::Element& ioletEl);
        lb::iolets::InOutLetCosine* DoIOForCosinePressureInOutlet(const io::xml::Element& ioletEl);
        lb::iolets::InOutLetFile* DoIOForFilePressureInOutlet(const io::xml::Element& ioletEl);
        lb::iolets::InOutLetMultiscale* DoIOForMultiscalePressureInOutlet(
            const io::xml::Element& ioletEl);

        lb::iolets::InOutLet* DoIOForVelocityInOutlet(const io::xml::Element& ioletEl);
        lb::iolets::InOutLetParabolicVelocity* DoIOForParabolicVelocityInOutlet(
            const io::xml::Element& ioletEl);
        /**
         * Reads a Womersley velocity iolet definition from the XML config file and returns
         * an InOutLetWomersleyVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetWomersleyVelocity object
         */
        lb::iolets::InOutLetWomersleyVelocity* DoIOForWomersleyVelocityInOutlet(
            const io::xml::Element& ioletEl);

        /**
         * Reads a file velocity iolet definition from the XML config file and returns
         * an InOutLetFileVelocity object
         *
         * @param ioletEl in memory representation of <inlet> or <outlet> xml element
         * @return InOutLetFileVelocity object
         */
        lb::iolets::InOutLetFileVelocity* DoIOForFileVelocityInOutlet(
            const io::xml::Element& ioletEl);

        void DoIOForProperties(const io::xml::Element& xmlNode);
        void DoIOForProperty(io::xml::Element xmlNode, bool isLoading);
        extraction::OutputField DoIOForPropertyField(const io::xml::Element& xmlNode);
        extraction::PropertyOutputFile* DoIOForPropertyOutputFile(
            const io::xml::Element& propertyoutputEl);
        extraction::StraightLineGeometrySelector* DoIOForLineGeometry(
            const io::xml::Element& xmlNode);
        extraction::PlaneGeometrySelector* DoIOForPlaneGeometry(const io::xml::Element&);
        extraction::SurfacePointSelector* DoIOForSurfacePoint(const io::xml::Element&);

        void DoIOForInitialConditions(io::xml::Element parent);
	void DoIOForCheckpointFile(const io::xml::Element& checkpointEl);
        void DoIOForVisualisation(const io::xml::Element& visEl);

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

        std::string xmlFilePath;
        io::xml::Document* rawXmlDoc;
        std::string dataFilePath;

        util::Vector3D<float> visualisationCentre;
        float visualisationLongitude;
        float visualisationLatitude;
        float visualisationZoom;
        float visualisationBrightness;
        float maxVelocity;
        float maxStress;
        lb::StressTypes stressType;
        std::vector<extraction::PropertyOutputFile*> propertyOutputs;
        std::string colloidConfigPath;
        /**
         * True if the file has a colloids section.
         */
        bool hasColloidSection;

        MonitoringConfig monitoringConfig; ///< Configuration of various checks/tests

        // We want to keep the RBC types isolated in their own library
        // so put all the configuration in a type that is opaque to
        // the rest of the code.
        //
        // nullptr => no RBC
        redblood::RBCConfig* rbcConf = nullptr;

      protected:
        // These have to contain pointers because there are multiple derived types that might be
        // instantiated.
        std::vector<lb::iolets::InOutLet*> inlets;
        std::vector<lb::iolets::InOutLet*> outlets;
        PhysicalTime timeStepSeconds;
        unsigned long totalTimeSteps;
        unsigned long warmUpSteps;
        PhysicalDistance voxelSizeMetres;
        PhysicalPosition geometryOriginMetres;
	PhysicalDensity fluidDensityKgm3;
	PhysicalDynamicViscosity fluidViscosityPas;
	PhysicalPressure reference_pressure_mmHg;
        util::UnitConverter* unitConverter;
        ICConfig icConfig;
      private:
    };
  }
}

#endif /* HEMELB_CONFIGURATION_SIMCONFIG_H */
