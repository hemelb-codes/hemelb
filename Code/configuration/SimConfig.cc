// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <string>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>

#include "configuration/SimConfig.h"
#include "reporting/BuildInfo.h"
#include "log/Logger.h"
#include "lb/InitialCondition.h"
#include "redblood/FlowExtension.h"
#include "redblood/RBCConfig.h"

namespace hemelb
{
  namespace configuration
  {
    // Base IC
    ICConfigBase::ICConfigBase(const util::UnitConverter* units, std::optional<LatticeTimeStep> t) : unitConverter(units), t0(t) {
    }

    // Uniform equilibrium IC
    EquilibriumIC::EquilibriumIC(const util::UnitConverter* units, std::optional<LatticeTimeStep> t, PhysicalPressure p) : ICConfigBase(units, t), p_mmHg(p), v_ms(0.0) {
    }

    EquilibriumIC::EquilibriumIC(const util::UnitConverter* units, std::optional<LatticeTimeStep> t, PhysicalPressure p, const PhysicalVelocity& v) : ICConfigBase(units, t), p_mmHg(p), v_ms(v) {
    }

    // checkpoint IC
    CheckpointIC::CheckpointIC(const util::UnitConverter* units, std::optional<LatticeTimeStep> t, const std::string& cp, std::optional<std::string> const& maybeOff) : ICConfigBase(units, t), cpFile(cp), maybeOffFile(maybeOff) {
    }


    std::unique_ptr<SimConfig> SimConfig::New(const path& path)
    {
      auto ans = std::unique_ptr<SimConfig>{new SimConfig{path}};
      ans->Init();
      return ans;
    }


    SimConfig::SimConfig(const path& path) :
        xmlFilePath(path)
    {
    }

    void SimConfig::Init()
    {
      if (!std::filesystem::exists(xmlFilePath))
      {
        throw Exception() << "Config file '" << xmlFilePath << "' does not exist";
      }
      rawXmlDoc = std::make_unique<io::xml::Document>(xmlFilePath);
      DoIO(rawXmlDoc->GetRoot());
    }

    SimConfig::~SimConfig()
    {
      delete unitConverter;
      unitConverter = nullptr;

      delete rbcConf;
      rbcConf = nullptr;
    }

    // Turn an input XML-relative path into a full path
    std::filesystem::path SimConfig::RelPathToFullPath(const std::string& path) const {
        auto xml_dir = xmlFilePath.parent_path();
        return std::filesystem::absolute(xml_dir / path);
    }

    void SimConfig::DoIO(io::xml::Element topNode)
    {
      // Top element must be:
      // <hemelbsettings version="5" />
      if (topNode.GetName() != "hemelbsettings")
        throw Exception() << "Invalid root element: " << topNode.GetPath();

      auto version = topNode.GetAttributeOrThrow<unsigned>("version");
      if (version != 5U)
        throw Exception() << "Unrecognised XML version. Expected 5, got " << version;

      DoIOForSimulation(topNode.GetChildOrThrow("simulation"));
      CreateUnitConverter();

      DoIOForGeometry(topNode.GetChildOrThrow("geometry"));

      if (topNode.GetChildOrNull("colloids") != io::xml::Element::Missing())
      {
        hasColloidSection = true;
      }

      DoIOForInitialConditions(topNode.GetChildOrThrow("initialconditions"));

      inlets = DoIOForInOutlets(topNode.GetChildOrThrow("inlets"));
      outlets = DoIOForInOutlets(topNode.GetChildOrThrow("outlets"));

      // Optional element <properties>
      if (auto propertiesEl = topNode.GetChildOrNull("properties"))
        DoIOForProperties(propertiesEl);

      // Optional element <monitoring>
      if (auto monitoringEl = topNode.GetChildOrNull("monitoring"))
        DoIOForMonitoring(monitoringEl);

      // The RBC section must be parsed *after* the inlets and outlets have been
      // defined
      if (auto rbcEl = topNode.GetChildOrNull("redbloodcells")) {
#ifdef HEMELB_BUILD_RBC
	rbcConf = new redblood::RBCConfig;
	rbcConf->DoIOForRedBloodCells(topNode, *this, *unitConverter);
#else
	throw Exception() << "Input XML has redbloodcells section but HEMELB_BUILD_RBC=OFF";
#endif
      }
    }

    void SimConfig::DoIOForSimulation(const io::xml::Element simEl)
    {
      // Required element
      // <stresstype value="enum lb::StressTypes" />
      auto tmp = simEl.GetChildOrThrow("stresstype").GetAttributeOrThrow<unsigned>("value");
      switch (tmp)
      {
        case lb::IgnoreStress:
          stressType = lb::IgnoreStress;
          break;
        case lb::ShearStress:
          stressType = lb::IgnoreStress;
          break;
        case lb::VonMises:
          stressType = lb::IgnoreStress;
          break;
        default:
          throw Exception() << "Invalid stresstype: " << tmp;
      }

      // Required element
      // <steps value="unsigned" units="lattice />
      const io::xml::Element stepsEl = simEl.GetChildOrThrow("steps");
      GetDimensionalValue(stepsEl, "lattice", totalTimeSteps);

      // Required element
      // <step_length value="float" units="s" />
      const io::xml::Element tsEl = simEl.GetChildOrThrow("step_length");
      GetDimensionalValue(tsEl, "s", timeStepSeconds);

      // Optional element
      // <extra_warmup_steps value="unsigned" units="lattice" />
      if (auto wuEl = simEl.GetChildOrNull("extra_warmup_steps"))
      {
        GetDimensionalValue(wuEl, "lattice", warmUpSteps);
        totalTimeSteps += warmUpSteps;
      }

      // Required element
      // <voxel_size value="float" units="m" />
      const io::xml::Element vsEl = simEl.GetChildOrThrow("voxel_size");
      GetDimensionalValue(vsEl, "m", voxelSizeMetres);

      // Required element
      // <origin value="(x,y,z)" units="m" />
      const io::xml::Element originEl = simEl.GetChildOrThrow("origin");
      GetDimensionalValue(originEl, "m", geometryOriginMetres);

      // Optional element
      // <fluid_density value="float" units="kg/m3" />
      fluidDensityKgm3 = simEl.GetChildOrNull("fluid_density").transform(
        [](io::xml::Element const& el) {
	  return GetDimensionalValue<PhysicalDensity>(el, "kg/m3");
	}).value_or(DEFAULT_FLUID_DENSITY_Kg_per_m3);

      // Optional element
      // <fluid_viscosity value="float" units="Pa.s" />
      fluidViscosityPas = simEl.GetChildOrNull("fluid_viscosity").transform(
        [](io::xml::Element const& el) {
	  return GetDimensionalValue<PhysicalDynamicViscosity>(el, "Pa.s");
	}).value_or(DEFAULT_FLUID_VISCOSITY_Pas);
     
      // Optional element (default = 0)
      // <reference_pressure value="float" units="mmHg" />
      reference_pressure_mmHg = simEl.GetChildOrNull("reference_pressure").transform(
        [](io::xml::Element const& el) {
	  return GetDimensionalValue<PhysicalPressure>(el, "mmHg");
	}).value_or(0);
    }

    void SimConfig::DoIOForGeometry(const io::xml::Element geometryEl)
    {
      // Required element
      // <geometry>
      //  <datafile path="relative path to GMY" />
      // </geometry>
      dataFilePath = RelPathToFullPath(geometryEl.GetChildOrThrow("datafile").GetAttributeOrThrow("path"));
    }

    void SimConfig::CreateUnitConverter()
    {
      unitConverter = new util::UnitConverter(timeStepSeconds,
                                              voxelSizeMetres,
                                              geometryOriginMetres,
					      fluidDensityKgm3,
					      reference_pressure_mmHg);
    }

    /**
     * Helper function to ensure that the iolet being created matches the compiled
     * iolet BC.
     * @param ioletEl
     * @param requiredBC
     */
    void SimConfig::CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                           const std::string& requiredBC)
    {
      // Check that HEMELB_*LET_BOUNDARY is consistent with this
      const std::string& ioletTypeName = ioletEl.GetName();
      std::string hemeIoletBC;

      if (ioletTypeName == "inlet")
        hemeIoletBC = reporting::inlet_boundary_condition;
      else if (ioletTypeName == "outlet")
        hemeIoletBC = reporting::outlet_boundary_condition;
      else
        throw Exception() << "Unexpected element name '" << ioletTypeName
            << "'. Expected 'inlet' or 'outlet'";

      if (requiredBC != hemeIoletBC)
      {
        throw Exception() << "XML configuration for " << ioletTypeName << " (line "
            << ioletEl.GetLine()
            << ") not consistent with compile-time choice of boundary condition '" << hemeIoletBC
            << "'";
      }
    }

    auto SimConfig::DoIOForInOutlets(const io::xml::Element ioletsEl) -> std::vector<IoletPtr>
    {
      const std::string& nodeName = ioletsEl.GetName();

      const std::string childNodeName = nodeName.substr(0, nodeName.size() - 1);
      std::vector<IoletPtr> ioletList;
      for (io::xml::Element currentIoletNode = ioletsEl.GetChildOrNull(childNodeName);
          currentIoletNode != io::xml::Element::Missing();
          currentIoletNode = currentIoletNode.NextSiblingOrNull(childNodeName))
      {
        // Determine which InOutlet to create
        io::xml::Element conditionEl = currentIoletNode.GetChildOrThrow("condition");
        const std::string& conditionType = conditionEl.GetAttributeOrThrow("type");

        IoletPtr newIolet;

        if (conditionType == "pressure")
        {
          newIolet = DoIOForPressureInOutlet(currentIoletNode);
        }
        else if (conditionType == "velocity")
        {
          newIolet = DoIOForVelocityInOutlet(currentIoletNode);
        }
        else
        {
          throw Exception() << "Invalid boundary condition type '" << conditionType << "' in "
              << conditionEl.GetPath();
        }
        DoIOForFlowExtension(newIolet.get(), currentIoletNode);
        ioletList.push_back(std::move(newIolet));
      }
      return ioletList;
    }

    void SimConfig::DoIOForFlowExtension(lb::iolets::InOutLet * iolet,
                                         const io::xml::Element & ioletNode)
    {
      if (ioletNode.GetChildOrNull("flowextension") == io::xml::Element::Missing())
        return;
      auto const flowExtension = redblood::readFlowExtension(ioletNode, GetUnitConverter());
      iolet->SetFlowExtension(std::make_shared<hemelb::redblood::FlowExtension>(flowExtension));
    }

    auto SimConfig::DoIOForPressureInOutlet(const io::xml::Element& ioletEl) -> IoletPtr
    {
      CheckIoletMatchesCMake(ioletEl, "NASHZEROTHORDERPRESSUREIOLET");
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      if (conditionSubtype == "cosine")
      {
        return DoIOForCosinePressureInOutlet(ioletEl);
      }
      else if (conditionSubtype == "file")
      {
        return DoIOForFilePressureInOutlet(ioletEl);
      }
      else if (conditionSubtype == "multiscale")
      {
        return DoIOForMultiscalePressureInOutlet(ioletEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << ioletEl.GetPath();
      }
    }

    auto SimConfig::DoIOForVelocityInOutlet(const io::xml::Element& ioletEl) -> IoletPtr
    {
      CheckIoletMatchesCMake(ioletEl, "LADDIOLET");
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      if (conditionSubtype == "parabolic")
      {
        return DoIOForParabolicVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "womersley")
      {
        return DoIOForWomersleyVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "file")
      {
        return DoIOForFileVelocityInOutlet(ioletEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << ioletEl.GetPath();
      }
    }

    void SimConfig::DoIOForProperties(const io::xml::Element& propertiesEl)
    {
      for (io::xml::ChildIterator poPtr = propertiesEl.IterChildren("propertyoutput");
          !poPtr.AtEnd(); ++poPtr)
      {
        propertyOutputs.push_back(DoIOForPropertyOutputFile(*poPtr));
      }

      if (auto cpEl = propertiesEl.GetChildOrNull("checkpoint")) {
	// Create a checkpoint property extractor.
	//
	// This is just a normal one, but fixed to be whole geometry,
	// only distributions, at double precision.
	//
	// Get stuff from XML
	auto file = extraction::PropertyOutputFile{};
	file.filename = cpEl.GetAttributeOrThrow("file");
	cpEl.GetAttributeOrThrow("period", file.frequency);
	// Configure the file
	file.geometry.reset(new extraction::WholeGeometrySelector());
	file.ts_mode = extraction::single_timestep_files{};

	extraction::OutputField field;
	field.name = "distributions";
	field.noffsets = 0;
	field.offset = {};
	field.typecode = distribn_t{0};
	field.src = extraction::source::Distributions{};
	file.fields.push_back(field);
	// Add to outputs
	propertyOutputs.push_back(std::move(file));
      }
    }

    extraction::PropertyOutputFile SimConfig::DoIOForPropertyOutputFile(
        const io::xml::Element& propertyoutputEl)
    {
      auto file = extraction::PropertyOutputFile{};

      auto&& ts_mode = propertyoutputEl.GetAttributeMaybe("timestep_mode").value_or("multi");
      if (ts_mode == "multi") {
	file.ts_mode = extraction::multi_timestep_file{};
      } else if (ts_mode == "single") {
	file.ts_mode = extraction::single_timestep_files{};
      } else {
	throw Exception()
	  << "Invalid value of timestep_mode attribute '" << ts_mode
	  << "' at: " << propertyoutputEl.GetPath();
      }

      file.filename = propertyoutputEl.GetAttributeOrThrow("file");
      if (std::holds_alternative<extraction::single_timestep_files>(file.ts_mode)) {
	auto const& p = file.filename.native();
	// Path must contain exactly one printf conversion specifier
	// for an integer.
	auto const errmsg = "For single timestep output files, "
	  "the path must contain exactly one '%d' and no other '%' characters";
	auto n_pc = std::count(p.begin(), p.end(), '%');
	auto i_pcd = p.find("%d", 0, 2);
	if (n_pc != 1 || i_pcd == std::string::npos) {
	  throw Exception() << errmsg;
	}
      }

      propertyoutputEl.GetAttributeOrThrow("period", file.frequency);

      io::xml::Element geometryEl = propertyoutputEl.GetChildOrThrow("geometry");
      const std::string& type = geometryEl.GetAttributeOrThrow("type");

      if (type == "plane")
      {
        file.geometry.reset(DoIOForPlaneGeometry(geometryEl));
      }
      else if (type == "line")
      {
        file.geometry.reset(DoIOForLineGeometry(geometryEl));
      }
      else if (type == "whole")
      {
        file.geometry.reset(new extraction::WholeGeometrySelector());
      }
      else if (type == "surface")
      {
        file.geometry.reset(new extraction::GeometrySurfaceSelector());
      }
      else if (type == "surfacepoint")
      {
        file.geometry.reset(DoIOForSurfacePoint(geometryEl));
      }
      else
      {
        throw Exception() << "Unrecognised property output geometry selector '" << type
            << "' in element " << geometryEl.GetPath();
      }

      for (io::xml::ChildIterator fieldPtr = propertyoutputEl.IterChildren("field");
          !fieldPtr.AtEnd(); ++fieldPtr)
        file.fields.push_back(DoIOForPropertyField(*fieldPtr));

      return file;
    }

    extraction::StraightLineGeometrySelector* SimConfig::DoIOForLineGeometry(
        const io::xml::Element& geometryEl)
    {
      io::xml::Element point1El = geometryEl.GetChildOrThrow("point");
      io::xml::Element point2El = point1El.NextSiblingOrThrow("point");

      PhysicalPosition point1;
      PhysicalPosition point2;

      GetDimensionalValue(point1El, "m", point1);
      GetDimensionalValue(point2El, "m", point2);

      return new extraction::StraightLineGeometrySelector(point1.as<float>(), point2.as<float>());
    }

    extraction::PlaneGeometrySelector* SimConfig::DoIOForPlaneGeometry(
        const io::xml::Element& geometryEl)
    {
      io::xml::Element pointEl = geometryEl.GetChildOrThrow("point");
      io::xml::Element normalEl = geometryEl.GetChildOrThrow("normal");

      PhysicalPosition point;
      util::Vector3D<float> normal;

      GetDimensionalValue(pointEl, "m", point);
      GetDimensionalValue(normalEl, "dimensionless", normal);

      auto radius = geometryEl.GetChildOrNull("radius").transform(
	[](io::xml::Element const& el) {
	  return GetDimensionalValue<PhysicalDistance>(el, "m");
	});
      if (radius)
        return new extraction::PlaneGeometrySelector(point.as<float>(), normal, *radius);
      else	
        return new extraction::PlaneGeometrySelector(point.as<float>(), normal);
    }

    extraction::SurfacePointSelector* SimConfig::DoIOForSurfacePoint(
        const io::xml::Element& geometryEl)
    {
      io::xml::Element pointEl = geometryEl.GetChildOrThrow("point");

      PhysicalPosition point;
      GetDimensionalValue(pointEl, "m", point);
      return new extraction::SurfacePointSelector(point.as<float>());
    }

    extraction::OutputField SimConfig::DoIOForPropertyField(const io::xml::Element& fieldEl)
    {
      extraction::OutputField field;
      const std::string& type = fieldEl.GetAttributeOrThrow("type");
      // Default name is identical to type.
      field.name = fieldEl.GetAttributeMaybe("name").value_or(type);

      // Default offset is none
      field.noffsets = 0;
      field.offset = {};

      // The default type to be written is float
      field.typecode = float{0.0};

      // Check and assign the type.
      if (type == "pressure")
      {
        field.src = extraction::source::Pressure{};
	// Pressure has an offset of the reference pressure
	field.noffsets = 1;
	field.offset = {reference_pressure_mmHg};
      }
      else if (type == "velocity")
      {
        field.src = extraction::source::Velocity{};
      }
      else if (type == "vonmisesstress")
      {
        field.src = extraction::source::VonMisesStress{};
      }
      else if (type == "shearstress")
      {
        field.src = extraction::source::ShearStress{};
      }
      else if (type == "shearrate")
      {
        field.src = extraction::source::ShearRate{};
      }
      else if (type == "stresstensor")
      {
        field.src = extraction::source::StressTensor{};
      }
      else if (type == "traction")
      {
        field.src = extraction::source::Traction{};
      }
      else if (type == "tangentialprojectiontraction")
      {
        field.src = extraction::source::TangentialProjectionTraction{};
      }
      else if (type == "distributions")
      {
        field.src = extraction::source::Distributions{};
      }
      else if (type == "mpirank")
      {
        field.src = extraction::source::MpiRank{};
	field.typecode = int{0};
      }
      else
      {
        throw Exception() << "Unrecognised field type '" << type << "' in " << fieldEl.GetPath();
      }
      return field;
    }

    void SimConfig::DoIOForBaseInOutlet(const io::xml::Element& ioletEl,
                                        lb::iolets::InOutLet* value)
    {
      io::xml::Element positionEl = ioletEl.GetChildOrThrow("position");
      io::xml::Element normalEl = ioletEl.GetChildOrThrow("normal");

      PhysicalPosition pos;
      GetDimensionalValue(positionEl, "m", pos);
      value->SetPosition(unitConverter->ConvertPositionToLatticeUnits(pos));

      util::Vector3D<double> norm;
      GetDimensionalValue(normalEl, "dimensionless", norm);
      value->SetNormal(norm);
    }

    void SimConfig::DoIOForInitialConditions(io::xml::Element initialconditionsEl)
    {
      // The <time> element may be present - if so, it will set the
      // initial timestep value
      auto t0 = initialconditionsEl.GetChildOrNull("time").transform([](auto el) {
	  return GetDimensionalValue<LatticeTimeStep>(el, "lattice");
	});

      // Exactly one of {<pressure>, <checkpoint>} must be present
      // TODO: use something other than an if-tree
      auto pressureEl = initialconditionsEl.GetChildOrNull("pressure");
      auto checkpointEl = initialconditionsEl.GetChildOrNull("checkpoint");
      if (pressureEl) {
	if (checkpointEl) {
	  // Both are present - this is an error
	  throw Exception()
	    << "XML contains both <pressure> and <checkpoint> sub elements of <initialconditions>";
	} else {
	  // Only pressure
	  io::xml::Element uniformEl = pressureEl.GetChildOrThrow("uniform");
	  PhysicalPressure p0_mmHg;
	  GetDimensionalValue(uniformEl, "mmHg", p0_mmHg);
	  icConfig = EquilibriumIC(unitConverter, t0, p0_mmHg);
	}
      } else {
	if (checkpointEl) {
	  // Only checkpoint
	  icConfig = CheckpointIC(unitConverter, t0,
				  checkpointEl.GetAttributeOrThrow("file"),
				  checkpointEl.GetAttributeMaybe("offsets"));
	} else {
	  // No IC!
	  throw Exception() << "XML <initialconditions> element contains no known initial condition type";
	}
      }
    }

    auto SimConfig::DoIOForCosinePressureInOutlet(
        const io::xml::Element& ioletEl) -> IoletPtr
    {
      auto newIolet = util::make_clone_ptr<lb::iolets::InOutLetCosine>();
      DoIOForBaseInOutlet(ioletEl, newIolet.get());

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      // Amplitude is a pressure DIFFERENCE (no use of REFERENCE_PRESSURE)
      auto p_amp = GetDimensionalValue<PhysicalPressure>(conditionEl.GetChildOrThrow("amplitude"), "mmHg");
      newIolet->SetPressureAmp(unitConverter->ConvertPressureDifferenceToLatticeUnits(p_amp));

      // Mean is an absolute pressure
      auto p_mean = GetDimensionalValue<PhysicalPressure>(conditionEl.GetChildOrThrow("mean"), "mmHg");
      newIolet->SetPressureMean(unitConverter->ConvertPressureToLatticeUnits(p_mean));

      auto phase = GetDimensionalValueInLatticeUnits<Angle>(conditionEl.GetChildOrThrow("phase"), "rad");
      newIolet->SetPhase(phase);

      auto period = GetDimensionalValueInLatticeUnits<LatticeTime>(conditionEl.GetChildOrThrow("period"), "s");
      newIolet->SetPeriod(period);

      if (warmUpSteps != 0)
      {
        newIolet->SetWarmup(warmUpSteps);
      }
      return newIolet;
    }

    auto SimConfig::DoIOForFilePressureInOutlet(
        const io::xml::Element& ioletEl) -> IoletPtr
    {
      auto newIolet = util::make_clone_ptr<lb::iolets::InOutLetFile>();
      DoIOForBaseInOutlet(ioletEl, newIolet.get());

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const io::xml::Element pathEl = conditionEl.GetChildOrThrow("path");
      newIolet->SetFilePath(pathEl.GetAttributeOrThrow("value"));

      return newIolet;
    }

    auto SimConfig::DoIOForMultiscalePressureInOutlet(
        const io::xml::Element& ioletEl) -> IoletPtr
    {
      auto newIolet = util::make_clone_ptr<lb::iolets::InOutLetMultiscale>();
      DoIOForBaseInOutlet(ioletEl, newIolet.get());

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      const io::xml::Element pressureEl = conditionEl.GetChildOrThrow("pressure");
      GetDimensionalValue(pressureEl, "mmHg", newIolet->GetPressureReference());

      const io::xml::Element velocityEl = conditionEl.GetChildOrThrow("velocity");
      GetDimensionalValue(velocityEl, "m/s", newIolet->GetVelocityReference());

      newIolet->GetLabel() = conditionEl.GetChildOrThrow("label").GetAttributeOrThrow("value");
      return newIolet;
    }

    auto SimConfig::DoIOForParabolicVelocityInOutlet(
        const io::xml::Element& ioletEl) -> IoletPtr
    {
      auto newIolet = util::make_clone_ptr<lb::iolets::InOutLetParabolicVelocity>();
      DoIOForBaseInOutlet(ioletEl, newIolet.get());

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      const io::xml::Element radiusEl = conditionEl.GetChildOrThrow("radius");
      newIolet->SetRadius(GetDimensionalValueInLatticeUnits < LatticeDistance > (radiusEl, "m"));

      const io::xml::Element maximumEl = conditionEl.GetChildOrThrow("maximum");
      newIolet->SetMaxSpeed(GetDimensionalValueInLatticeUnits < PhysicalSpeed > (maximumEl, "m/s"));

      if (warmUpSteps != 0)
      {
        newIolet->SetWarmup(warmUpSteps);
      }

      return newIolet;
    }

    auto SimConfig::DoIOForWomersleyVelocityInOutlet(
        const io::xml::Element& ioletEl) -> IoletPtr
    {
      auto newIolet = util::make_clone_ptr<lb::iolets::InOutLetWomersleyVelocity>();
      DoIOForBaseInOutlet(ioletEl, newIolet.get());

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      const io::xml::Element radiusEl = conditionEl.GetChildOrThrow("radius");
      newIolet->SetRadius(GetDimensionalValueInLatticeUnits < LatticeDistance > (radiusEl, "m"));

      const io::xml::Element pgAmpEl = conditionEl.GetChildOrThrow("pressure_gradient_amplitude");
      newIolet->SetPressureGradientAmplitude(GetDimensionalValueInLatticeUnits
          < LatticePressureGradient > (pgAmpEl, "mmHg/m"));

      const io::xml::Element periodEl = conditionEl.GetChildOrThrow("period");
      newIolet->SetPeriod(GetDimensionalValueInLatticeUnits < LatticeTime > (periodEl, "s"));

      const io::xml::Element womNumEl = conditionEl.GetChildOrThrow("womersley_number");
      newIolet->SetWomersleyNumber(GetDimensionalValueInLatticeUnits < Dimensionless
          > (womNumEl, "dimensionless"));

      return newIolet;
    }

    auto SimConfig::DoIOForFileVelocityInOutlet(
        const io::xml::Element& ioletEl) -> IoletPtr
    {
      auto newIolet = util::make_clone_ptr<lb::iolets::InOutLetFileVelocity>();
      DoIOForBaseInOutlet(ioletEl, newIolet.get());

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      std::string velocityFilePath = RelPathToFullPath(conditionEl.GetChildOrThrow("path").GetAttributeOrThrow("value"));
      newIolet->SetFilePath(velocityFilePath);

      const io::xml::Element radiusEl = conditionEl.GetChildOrThrow("radius");
      newIolet->SetRadius(GetDimensionalValueInLatticeUnits < LatticeDistance > (radiusEl, "m"));

      return newIolet;
    }

    bool SimConfig::HasColloidSection() const
    {
      return hasColloidSection;
    }

    const util::UnitConverter& SimConfig::GetUnitConverter() const
    {
      if (unitConverter == nullptr)
        throw Exception() << "Invalid UnitConverter (nullptr)";

      return *unitConverter;
    }

    void SimConfig::DoIOForMonitoring(const io::xml::Element& monEl)
    {
      io::xml::Element convEl = monEl.GetChildOrNull("steady_flow_convergence");
      if (convEl != io::xml::Element::Missing())
      {
        DoIOForSteadyFlowConvergence(convEl);
      }

      monitoringConfig.doIncompressibilityCheck = (monEl.GetChildOrNull("incompressibility")
          != io::xml::Element::Missing());
    }

    void SimConfig::DoIOForSteadyFlowConvergence(const io::xml::Element& convEl)
    {
      monitoringConfig.doConvergenceCheck = true;
      monitoringConfig.convergenceRelativeTolerance = convEl.GetAttributeOrThrow<double>("tolerance");
      monitoringConfig.convergenceTerminate = (convEl.GetAttributeOrThrow("terminate") == "true");

      if (convEl.IterChildren("criterion").AtEnd())
      {
        throw Exception() << "At least one convergence criterion must be provided in "
            << convEl.GetPath();
      }

      for (io::xml::ChildIterator criteriaIt = convEl.IterChildren("criterion");
          !criteriaIt.AtEnd(); ++criteriaIt)
      {
        DoIOForConvergenceCriterion(*criteriaIt);
      }
    }

    void SimConfig::DoIOForConvergenceCriterion(const io::xml::Element& criterionEl)
    {
      const std::string& criterionType = criterionEl.GetAttributeOrThrow("type");

      // We only allow velocity-based convergence check for the time being
      if (criterionType != "velocity")
      {
        throw Exception() << "Invalid convergence criterion type " << criterionType << " in "
            << criterionEl.GetPath();
      }
      monitoringConfig.convergenceVariable = extraction::source::Velocity{};
      monitoringConfig.convergenceReferenceValue = GetDimensionalValueInLatticeUnits < LatticeSpeed
          > (criterionEl, "m/s");
    }

    const MonitoringConfig* SimConfig::GetMonitoringConfiguration() const
    {
      return &monitoringConfig;
    }
  }
}
