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
#include "util/fileutils.h"
#include "lb/InitialCondition.h"
#include "redblood/FlowExtension.h"
#include "redblood/RBCConfig.h"

namespace hemelb
{
  namespace configuration
  {
    // Base IC
    ICConfigBase::ICConfigBase(const util::UnitConverter* units, boost::optional<LatticeTimeStep> t) : unitConverter(units), t0(t) {
    }

    // Uniform equilibrium IC
    EquilibriumIC::EquilibriumIC(const util::UnitConverter* units, boost::optional<LatticeTimeStep> t, PhysicalPressure p) : ICConfigBase(units, t), p_mmHg(p), v_ms(0.0) {
    }

    EquilibriumIC::EquilibriumIC(const util::UnitConverter* units, boost::optional<LatticeTimeStep> t, PhysicalPressure p, const PhysicalVelocity& v) : ICConfigBase(units, t), p_mmHg(p), v_ms(v) {
    }

    // checkpoint IC
    CheckpointIC::CheckpointIC(const util::UnitConverter* units, boost::optional<LatticeTimeStep> t, const std::string& cp) : ICConfigBase(units, t), cpFile(cp) {
    }


    SimConfig* SimConfig::New(const std::string& path)
    {
      SimConfig* ans = new SimConfig(path);
      ans->Init();
      return ans;
    }


    SimConfig::SimConfig(const std::string& path) :
        xmlFilePath(path), rawXmlDoc(nullptr), hasColloidSection(false), warmUpSteps(0),
            unitConverter(nullptr)
    {
    }

    void SimConfig::Init()
    {
      if (!util::file_exists(xmlFilePath.c_str()))
      {
        throw Exception() << "Config file '" << xmlFilePath << "' does not exist";
      }
      rawXmlDoc = new io::xml::Document(xmlFilePath);
      colloidConfigPath = xmlFilePath;
      DoIO(rawXmlDoc->GetRoot());
    }

    SimConfig::~SimConfig()
    {
      for (unsigned outputNumber = 0; outputNumber < propertyOutputs.size(); ++outputNumber)
      {
        delete propertyOutputs[outputNumber];
      }

      delete rawXmlDoc;
      rawXmlDoc = nullptr;

      delete unitConverter;
      unitConverter = nullptr;

      delete rbcConf;
      rbcConf = nullptr;
    }

    // Turn an input XML-relative path into a full path
    std::string SimConfig::RelPathToFullPath(const std::string& path) const {
      return util::NormalizePathRelativeToPath(path, xmlFilePath);
    }

    void SimConfig::DoIO(io::xml::Element topNode)
    {
      // Top element must be:
      // <hemelbsettings version="4" />
      if (topNode.GetName() != "hemelbsettings")
        throw Exception() << "Invalid root element: " << topNode.GetPath();

      unsigned version;
      const std::string& versionStr = topNode.GetAttributeOrThrow("version", version);
      if (version != 4U)
        throw Exception() << "Unrecognised XML version. Expected 4, got " << versionStr;

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

      DoIOForVisualisation(topNode.GetChildOrThrow("visualisation"));

      // Optional element <properties>
      io::xml::Element propertiesEl = topNode.GetChildOrNull("properties");
      if (propertiesEl != io::xml::Element::Missing())
        DoIOForProperties(propertiesEl);

      // Optional element <monitoring>
      io::xml::Element monitoringEl = topNode.GetChildOrNull("monitoring");
      if (monitoringEl != io::xml::Element::Missing())
        DoIOForMonitoring(monitoringEl);

      // The RBC section must be parsed *after* the inlets and outlets have been
      // defined
      auto rbcEl = topNode.GetChildOrNull("redbloodcells");
      if (rbcEl != io::xml::Element::Missing()) {
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
      unsigned tmp;
      simEl.GetChildOrThrow("stresstype").GetAttributeOrThrow("value", tmp);
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
      const io::xml::Element wuEl = simEl.GetChildOrNull("extra_warmup_steps");
      if (wuEl != io::xml::Element::Missing())
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
      auto maybeDensityEl = simEl.GetChildOrNull("fluid_density");
      if (maybeDensityEl) {
	GetDimensionalValue(maybeDensityEl, "kg/m3", fluidDensityKgm3);
      } else {
	fluidDensityKgm3 = DEFAULT_FLUID_DENSITY_Kg_per_m3;
      }

      // Optional element
      // <fluid_viscosity value="float" units="Pa.s" />
      auto maybeViscosityEl = simEl.GetChildOrNull("fluid_viscosity");
      if (maybeViscosityEl) {
	GetDimensionalValue(maybeViscosityEl, "Pa.s", fluidViscosityPas);
      } else {
	fluidViscosityPas = DEFAULT_FLUID_VISCOSITY_Pas;
      }

      // Optional element (default = 0)
      // <reference_pressure value="float" units="mmHg" />
      auto maybeRefPresEl = simEl.GetChildOrNull("reference_pressure");
      if (maybeRefPresEl) {
	GetDimensionalValue(maybeRefPresEl, "mmHg", reference_pressure_mmHg);
      } else {
	reference_pressure_mmHg = 0;
      }
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

    std::vector<lb::iolets::InOutLet*> SimConfig::DoIOForInOutlets(const io::xml::Element ioletsEl)
    {
      const std::string& nodeName = ioletsEl.GetName();

      const std::string childNodeName = nodeName.substr(0, nodeName.size() - 1);
      std::vector<lb::iolets::InOutLet*> ioletList;
      for (io::xml::Element currentIoletNode = ioletsEl.GetChildOrNull(childNodeName);
          currentIoletNode != io::xml::Element::Missing();
          currentIoletNode = currentIoletNode.NextSiblingOrNull(childNodeName))
      {
        // Determine which InOutlet to create
        io::xml::Element conditionEl = currentIoletNode.GetChildOrThrow("condition");
        const std::string& conditionType = conditionEl.GetAttributeOrThrow("type");

        lb::iolets::InOutLet* newIolet = nullptr;

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
        DoIOForFlowExtension(newIolet, currentIoletNode);
        ioletList.push_back(newIolet);
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

    lb::iolets::InOutLet* SimConfig::DoIOForPressureInOutlet(const io::xml::Element& ioletEl)
    {
      CheckIoletMatchesCMake(ioletEl, "NASHZEROTHORDERPRESSUREIOLET");
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      lb::iolets::InOutLet* newIolet = nullptr;
      if (conditionSubtype == "cosine")
      {
        newIolet = DoIOForCosinePressureInOutlet(ioletEl);
      }
      else if (conditionSubtype == "file")
      {
        newIolet = DoIOForFilePressureInOutlet(ioletEl);
      }
      else if (conditionSubtype == "multiscale")
      {
        newIolet = DoIOForMultiscalePressureInOutlet(ioletEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << ioletEl.GetPath();
      }

      return newIolet;
    }

    lb::iolets::InOutLet* SimConfig::DoIOForVelocityInOutlet(const io::xml::Element& ioletEl)
    {
      CheckIoletMatchesCMake(ioletEl, "LADDIOLET");
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      lb::iolets::InOutLet* newIolet = nullptr;
      if (conditionSubtype == "parabolic")
      {
        newIolet = DoIOForParabolicVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "womersley")
      {
        newIolet = DoIOForWomersleyVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "file")
      {
        newIolet = DoIOForFileVelocityInOutlet(ioletEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << ioletEl.GetPath();
      }

      return newIolet;
    }
    void SimConfig::DoIOForVisualisation(const io::xml::Element& visEl)
    {
      GetDimensionalValue(visEl.GetChildOrThrow("centre"), "m", visualisationCentre);

      io::xml::Element orientationEl = visEl.GetChildOrThrow("orientation");
      GetDimensionalValue(orientationEl.GetChildOrThrow("longitude"),
                          "deg",
                          visualisationLongitude);
      GetDimensionalValue(orientationEl.GetChildOrThrow("latitude"), "deg", visualisationLatitude);

      io::xml::Element displayEl = visEl.GetChildOrThrow("display");
      displayEl.GetAttributeOrThrow("zoom", visualisationZoom);
      displayEl.GetAttributeOrThrow("brightness", visualisationBrightness);

      io::xml::Element rangeEl = visEl.GetChildOrThrow("range");
      GetDimensionalValue(rangeEl.GetChildOrThrow("maxvelocity"), "m/s", maxVelocity);
      GetDimensionalValue(rangeEl.GetChildOrThrow("maxstress"), "Pa", maxStress);
    }

    void SimConfig::DoIOForProperties(const io::xml::Element& propertiesEl)
    {
      for (io::xml::ChildIterator poPtr = propertiesEl.IterChildren("propertyoutput");
          !poPtr.AtEnd(); ++poPtr)
      {
        propertyOutputs.push_back(DoIOForPropertyOutputFile(*poPtr));
      }
    }

    extraction::PropertyOutputFile* SimConfig::DoIOForPropertyOutputFile(
        const io::xml::Element& propertyoutputEl)
    {
      extraction::PropertyOutputFile* file = new extraction::PropertyOutputFile();
      file->filename = propertyoutputEl.GetAttributeOrThrow("file");

      propertyoutputEl.GetAttributeOrThrow("period", file->frequency);

      io::xml::Element geometryEl = propertyoutputEl.GetChildOrThrow("geometry");
      const std::string& type = geometryEl.GetAttributeOrThrow("type");

      if (type == "plane")
      {
        file->geometry.reset(DoIOForPlaneGeometry(geometryEl));
      }
      else if (type == "line")
      {
        file->geometry.reset(DoIOForLineGeometry(geometryEl));
      }
      else if (type == "whole")
      {
        file->geometry.reset(new extraction::WholeGeometrySelector());
      }
      else if (type == "surface")
      {
        file->geometry.reset(new extraction::GeometrySurfaceSelector());
      }
      else if (type == "surfacepoint")
      {
        file->geometry.reset(DoIOForSurfacePoint(geometryEl));
      }
      else
      {
        throw Exception() << "Unrecognised property output geometry selector '" << type
            << "' in element " << geometryEl.GetPath();
      }

      for (io::xml::ChildIterator fieldPtr = propertyoutputEl.IterChildren("field");
          !fieldPtr.AtEnd(); ++fieldPtr)
        file->fields.push_back(DoIOForPropertyField(*fieldPtr));

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

      io::xml::Element radiusEl = geometryEl.GetChildOrNull("radius");

      if (radiusEl == io::xml::Element::Missing())
      {
        return new extraction::PlaneGeometrySelector(point.as<float>(), normal);
      }
      else
      {
        PhysicalDistance radius;
        GetDimensionalValue(radiusEl, "m", radius);
        return new extraction::PlaneGeometrySelector(point.as<float>(), normal, radius);
      }

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
      const std::string type = fieldEl.GetAttributeOrThrow("type");
      const std::string* name = fieldEl.GetAttributeOrNull("name");

      // Default name is identical to type.
      if (name == nullptr)
      {
        field.name = type;
      }
      else
      {
        field.name = *name;
      }

      // Default offset is zero
      field.offset = 0.0;

      // Check and assign the type.
      if (type == "pressure")
      {
        field.type = extraction::OutputField::Pressure;
	// Pressure uses the reference pressure
	field.offset = reference_pressure_mmHg;
      }
      else if (type == "velocity")
      {
        field.type = extraction::OutputField::Velocity;
      }
      else if (type == "vonmisesstress")
      {
        field.type = extraction::OutputField::VonMisesStress;
      }
      else if (type == "shearstress")
      {
        field.type = extraction::OutputField::ShearStress;
      }
      else if (type == "shearrate")
      {
        field.type = extraction::OutputField::ShearRate;
      }
      else if (type == "stresstensor")
      {
        field.type = extraction::OutputField::StressTensor;
      }
      else if (type == "traction")
      {
        field.type = extraction::OutputField::Traction;
      }
      else if (type == "tangentialprojectiontraction")
      {
        field.type = extraction::OutputField::TangentialProjectionTraction;
      }
      else if (type == "distributions")
      {
        field.type = extraction::OutputField::Distributions;
      }
      else if (type == "mpirank")
      {
        field.type = extraction::OutputField::MpiRank;
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
      boost::optional<LatticeTimeStep> t0;
      if (auto timeEl = initialconditionsEl.GetChildOrNull("time")) {
	LatticeTimeStep tmp;
	GetDimensionalValue(timeEl, "lattice", tmp);
	t0 = tmp;
      }

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
	  icConfig = CheckpointIC(unitConverter, t0, checkpointEl.GetAttributeOrThrow("file"));
	} else {
	  // No IC!
	  throw Exception() << "XML <initialconditions> element contains no known initial condition type";
	}
      }
    }

    lb::iolets::InOutLetCosine* SimConfig::DoIOForCosinePressureInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetCosine* newIolet = new lb::iolets::InOutLetCosine();
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      PhysicalPressure tempP;
      // Amplitude is a pressure DIFFERENCE (no use of REFERENCE_PRESSURE)
      GetDimensionalValue(conditionEl.GetChildOrThrow("amplitude"), "mmHg", tempP);
      newIolet->SetPressureAmp(unitConverter->ConvertPressureDifferenceToLatticeUnits(tempP));

      // Mean is an absolute pressure
      GetDimensionalValue(conditionEl.GetChildOrThrow("mean"), "mmHg", tempP);
      newIolet->SetPressureMean(unitConverter->ConvertPressureToLatticeUnits(tempP));

      newIolet->SetPhase(GetDimensionalValueInLatticeUnits < Angle
          > (conditionEl.GetChildOrThrow("phase"), "rad"));

      LatticeTime period;
      GetDimensionalValueInLatticeUnits(conditionEl.GetChildOrThrow("period"), "s", period);
      newIolet->SetPeriod(period);

      if (warmUpSteps != 0)
      {
        newIolet->SetWarmup(warmUpSteps);
      }
      return newIolet;
    }

    lb::iolets::InOutLetFile* SimConfig::DoIOForFilePressureInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetFile* newIolet = new lb::iolets::InOutLetFile();
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const io::xml::Element pathEl = conditionEl.GetChildOrThrow("path");
      newIolet->SetFilePath(pathEl.GetAttributeOrThrow("value"));

      return newIolet;
    }

    lb::iolets::InOutLetMultiscale* SimConfig::DoIOForMultiscalePressureInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetMultiscale* newIolet = new lb::iolets::InOutLetMultiscale();
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      const io::xml::Element pressureEl = conditionEl.GetChildOrThrow("pressure");
      GetDimensionalValue(pressureEl, "mmHg", newIolet->GetPressureReference());

      const io::xml::Element velocityEl = conditionEl.GetChildOrThrow("velocity");
      GetDimensionalValue(velocityEl, "m/s", newIolet->GetVelocityReference());

      newIolet->GetLabel() = conditionEl.GetChildOrThrow("label").GetAttributeOrThrow("value");
      return newIolet;
    }

    lb::iolets::InOutLetParabolicVelocity* SimConfig::DoIOForParabolicVelocityInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetParabolicVelocity* newIolet = new lb::iolets::InOutLetParabolicVelocity();
      DoIOForBaseInOutlet(ioletEl, newIolet);

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

    lb::iolets::InOutLetWomersleyVelocity* SimConfig::DoIOForWomersleyVelocityInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetWomersleyVelocity* newIolet = new lb::iolets::InOutLetWomersleyVelocity();
      DoIOForBaseInOutlet(ioletEl, newIolet);

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

    lb::iolets::InOutLetFileVelocity* SimConfig::DoIOForFileVelocityInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetFileVelocity* newIolet = new lb::iolets::InOutLetFileVelocity();
      DoIOForBaseInOutlet(ioletEl, newIolet);

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
      convEl.GetAttributeOrThrow("tolerance", monitoringConfig.convergenceRelativeTolerance);
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
      monitoringConfig.convergenceVariable = extraction::OutputField::Velocity;
      monitoringConfig.convergenceReferenceValue = GetDimensionalValueInLatticeUnits < LatticeSpeed
          > (criterionEl, "m/s");
    }

    const MonitoringConfig* SimConfig::GetMonitoringConfiguration() const
    {
      return &monitoringConfig;
    }
  }
}
