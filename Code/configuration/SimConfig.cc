
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
#include "log/Logger.h"
#include "util/fileutils.h"

namespace hemelb
{
  namespace configuration
  {
    SimConfig* SimConfig::New(const std::string& path)
    {
      SimConfig* ans = new SimConfig(path);
      ans->Init();
      return ans;
    }

    SimConfig::SimConfig(const std::string& path) :
        xmlFilePath(path), rawXmlDoc(NULL), hasColloidSection(false), warmUpSteps(0),
            unitConverter(NULL)
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
      dataFilePath = util::NormalizePathRelativeToPath(dataFilePath, xmlFilePath);
    }

    SimConfig::~SimConfig()
    {
      for (unsigned outputNumber = 0; outputNumber < propertyOutputs.size(); ++outputNumber)
      {
        delete propertyOutputs[outputNumber];
      }

      delete rawXmlDoc;
      rawXmlDoc = NULL;

      delete unitConverter;
      unitConverter = NULL;
    }

    void SimConfig::DoIO(io::xml::Element topNode)
    {
      // Top element must be:
      // <hemelbsettings version="2" />
      if (topNode.GetName() != "hemelbsettings")
        throw Exception() << "Invalid root element: " << topNode.GetPath();

      unsigned version;
      const std::string& versionStr = topNode.GetAttributeOrThrow("version", version);
      if (version != 3U)
        throw Exception() << "Unrecognised XML version. Expected 3, got " << versionStr;

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

      stents = DoIOForStents(topNode.GetChildOrThrow("stents"));

      DoIOForVisualisation(topNode.GetChildOrThrow("visualisation"));

      // Optional element <properties>
      io::xml::Element propertiesEl = topNode.GetChildOrNull("properties");
      if (propertiesEl != io::xml::Element::Missing())
        DoIOForProperties(propertiesEl);

      // Optional element <monitoring>
      io::xml::Element monitoringEl = topNode.GetChildOrNull("monitoring");
      if (monitoringEl != io::xml::Element::Missing())
        DoIOForMonitoring(monitoringEl);

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
    }

    void SimConfig::DoIOForGeometry(const io::xml::Element geometryEl)
    {
      // Required element
      // <geometry>
      //  <datafile path="relative path to GMY" />
      // </geometry>
      dataFilePath = geometryEl.GetChildOrThrow("datafile").GetAttributeOrThrow("path");
      // Convert to a full path
      dataFilePath = util::NormalizePathRelativeToPath(dataFilePath, xmlFilePath);

    }

    void SimConfig::CreateUnitConverter()
    {
      unitConverter = new util::UnitConverter(timeStepSeconds,
                                              voxelSizeMetres,
                                              geometryOriginMetres);
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

#define QUOTE_RAW(x) #x
#define QUOTE_CONTENTS(x) QUOTE_RAW(x)
      if (ioletTypeName == "inlet")
        hemeIoletBC = QUOTE_CONTENTS(HEMELB_INLET_BOUNDARY);
      else if (ioletTypeName == "outlet")
        hemeIoletBC = QUOTE_CONTENTS(HEMELB_OUTLET_BOUNDARY);
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

        lb::iolets::InOutLet* newIolet = NULL;

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
        ioletList.push_back(newIolet);
      }
      return ioletList;
    }

    std::vector<lb::stents::Stent*> SimConfig::DoIOForStents(const io::xml::Element stentsEl)
    {
      const std::string& nodeName = stentsEl.GetName();

      const std::string childNodeName = nodeName.substr(0, nodeName.size() - 1);
      std::vector<lb::stents::Stent*> stentList;
      for (io::xml::Element currentStentNode = stentsEl.GetChildOrNull(childNodeName);
          currentStentNode != io::xml::Element::Missing();
          currentStentNode = currentStentNode.NextSiblingOrNull(childNodeName))
      {
        // Determine which InOutlet to create
        io::xml::Element conditionEl = currentStentNode.GetChildOrThrow("condition");
        const std::string& conditionType = conditionEl.GetAttributeOrThrow("type");

        lb::stents::Stent* newStent = NULL;

        if (conditionType == "concentration")
        {
          newStent = DoIOForDensityStent(currentStentNode);
        }
        else if (conditionType == "flux")
        {
          newStent = DoIOForFluxStent(currentStentNode);
        }
        else
        {
          throw Exception() << "Invalid boundary condition type '" << conditionType << "' in "
              << conditionEl.GetPath();
        }
        stentList.push_back(newStent);
      }
      return stentList;
    }

    lb::iolets::InOutLet* SimConfig::DoIOForPressureInOutlet(const io::xml::Element& ioletEl)
    {
      CheckIoletMatchesCMake(ioletEl, "NASHZEROTHORDERPRESSUREIOLET");
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      lb::iolets::InOutLet* newIolet = NULL;
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

    lb::stents::Stent* SimConfig::DoIOForDensityStent(const io::xml::Element& stentEl)
    {
      io::xml::Element conditionEl = stentEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      lb::stents::Stent* newStent = NULL;
      if (conditionSubtype == "constant")
      {
        newStent = DoIOForConstantDensityStent(stentEl);
      }
      else if (conditionSubtype == "file")
      {
        newStent = DoIOForFileDensityStent(stentEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << stentEl.GetPath();
      }

      return newStent;
    }

    lb::iolets::InOutLet* SimConfig::DoIOForVelocityInOutlet(const io::xml::Element& ioletEl)
    {
      CheckIoletMatchesCMake(ioletEl, "LADDIOLET");
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      lb::iolets::InOutLet* newIolet = NULL;
      if (conditionSubtype == "parabolic")
      {
        newIolet = DoIOForParabolicVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "hollowparabolic")
      {
        newIolet = DoIOForHollowParabolicVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "womersley")
      {
        newIolet = DoIOForWomersleyVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "file")
      {
        newIolet = DoIOForFileVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "hollowfile")
      {
        newIolet = DoIOForHollowFileVelocityInOutlet(ioletEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << ioletEl.GetPath();
      }

      return newIolet;
    }

    lb::stents::Stent* SimConfig::DoIOForFluxStent(const io::xml::Element& stentEl)
    {
      io::xml::Element conditionEl = stentEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      lb::stents::Stent* newStent = NULL;
      if (conditionSubtype == "constant")
      {
        newStent = DoIOForConstantFluxStent(stentEl);
      }
      else if (conditionSubtype == "file")
      {
        newStent = DoIOForFileFluxStent(stentEl);
      }
      else
      {
        throw Exception() << "Invalid boundary condition subtype '" << conditionSubtype << "' in "
            << stentEl.GetPath();
      }

      return newStent;
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
        file->geometry = DoIOForPlaneGeometry(geometryEl);
      }
      else if (type == "line")
      {
        file->geometry = DoIOForLineGeometry(geometryEl);
      }
      else if (type == "whole")
      {
        file->geometry = new extraction::WholeGeometrySelector();
      }
      else if (type == "surface")
      {
        file->geometry = new extraction::GeometrySurfaceSelector();
      }
      else if (type == "surfacepoint")
      {
        file->geometry = DoIOForSurfacePoint(geometryEl);
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

      return new extraction::StraightLineGeometrySelector(point1, point2);
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
        return new extraction::PlaneGeometrySelector(point, normal);
      }
      else
      {
        PhysicalDistance radius;
        GetDimensionalValue(radiusEl, "m", radius);
        return new extraction::PlaneGeometrySelector(point, normal, radius);
      }

    }

    extraction::SurfacePointSelector* SimConfig::DoIOForSurfacePoint(
        const io::xml::Element& geometryEl)
    {
      io::xml::Element pointEl = geometryEl.GetChildOrThrow("point");

      PhysicalPosition point;
      GetDimensionalValue(pointEl, "m", point);
      return new extraction::SurfacePointSelector(point);
    }

    extraction::OutputField SimConfig::DoIOForPropertyField(const io::xml::Element& fieldEl)
    {
      extraction::OutputField field;
      const std::string type = fieldEl.GetAttributeOrThrow("type");
      const std::string* name = fieldEl.GetAttributeOrNull("name");

      // Default name is identical to type.
      if (name == NULL)
      {
        field.name = type;
      }
      else
      {
        field.name = *name;
      }

      // Check and assign the type.
      if (type == "pressure")
      {
        field.type = extraction::OutputField::Pressure;
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
      else if (type == "mpirank")
      {
        field.type = extraction::OutputField::MpiRank;
      }
      else if (type == "concentration") {
	field.type = extraction::OutputField::Concentration;
      }
      else if (type == "flux") {
        field.type = extraction::OutputField::Flux;
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
      //, isLoading, initialPressure
      io::xml::Element pressureEl = initialconditionsEl.GetChildOrThrow("pressure");
      io::xml::Element uniformEl = pressureEl.GetChildOrThrow("uniform");

      GetDimensionalValue(uniformEl, "mmHg", initialPressure_mmHg);

      //, isLoading, initialDensity
      io::xml::Element densityEl = initialconditionsEl.GetChildOrThrow("concentration");
      io::xml::Element uniformDensityEl = densityEl.GetChildOrThrow("uniform");

      GetDimensionalValue(uniformDensityEl, "kg/m3", initialDensity_Kg_per_m3);
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

      newIolet->SetPhase(GetDimensionalValueInLatticeUnits<Angle>(conditionEl.GetChildOrThrow("phase"),
                                                                  "rad"));

      LatticeTime period;
      GetDimensionalValueInLatticeUnits(conditionEl.GetChildOrThrow("period"), "s", period);
      newIolet->SetPeriod(period);

      if (warmUpSteps != 0)
      {
        newIolet->SetWarmup(warmUpSteps);
      }
      return newIolet;
    }

    lb::stents::StentConstant* SimConfig::DoIOForConstantDensityStent(
        const io::xml::Element& stentEl)
    {
      lb::stents::StentConstant* newStent = new lb::stents::StentConstant();

      const io::xml::Element conditionEl = stentEl.GetChildOrThrow("condition");

      PhysicalDensity tempD;
      // Amplitude is a pressure DIFFERENCE (no use of REFERENCE_PRESSURE)
      GetDimensionalValue(conditionEl.GetChildOrThrow("mean"), "kg/m3", tempD);
      newStent->SetDensity(unitConverter->ConvertDensityToLatticeUnits(tempD));

      return newStent;
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

    lb::stents::StentFile* SimConfig::DoIOForFileDensityStent(
        const io::xml::Element& stentEl)
    {
      lb::stents::StentFile* newStent = new lb::stents::StentFile();

      const io::xml::Element conditionEl = stentEl.GetChildOrThrow("condition");
      const io::xml::Element pathEl = conditionEl.GetChildOrThrow("path");
      newStent->SetFilePath(pathEl.GetAttributeOrThrow("value"));

      return newStent;
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
      newIolet->SetRadius(GetDimensionalValueInLatticeUnits<LatticeDistance>(radiusEl, "m"));

      const io::xml::Element maximumEl = conditionEl.GetChildOrThrow("maximum");
      newIolet->SetMaxSpeed(GetDimensionalValueInLatticeUnits<PhysicalSpeed>(maximumEl, "m/s"));

      const io::xml::Element concentrationEl = conditionEl.GetChildOrThrow("concentration");
      newIolet->SetConcentration(GetDimensionalValueInLatticeUnits<PhysicalDensity>(concentrationEl, "kg/m3"));

      if (warmUpSteps != 0)
      {
        newIolet->SetWarmup(warmUpSteps);
      }

      return newIolet;
    }

    lb::iolets::InOutLetHollowParabolicVelocity* SimConfig::DoIOForHollowParabolicVelocityInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetHollowParabolicVelocity* newIolet = new lb::iolets::InOutLetHollowParabolicVelocity();
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      const io::xml::Element radiusEl = conditionEl.GetChildOrThrow("radius");
      newIolet->SetRadius(GetDimensionalValueInLatticeUnits<LatticeDistance>(radiusEl, "m"));

      const io::xml::Element innerRadiusEl = conditionEl.GetChildOrThrow("innerradius");
      newIolet->SetInnerRadius(GetDimensionalValueInLatticeUnits<LatticeDistance>(innerRadiusEl, "m"));

      const io::xml::Element maximumEl = conditionEl.GetChildOrThrow("maximum");
      newIolet->SetMaxSpeed(GetDimensionalValueInLatticeUnits<PhysicalSpeed>(maximumEl, "m/s"));

      const io::xml::Element concentrationEl = conditionEl.GetChildOrThrow("concentration");
      newIolet->SetConcentration(GetDimensionalValueInLatticeUnits<PhysicalDensity>(concentrationEl, "kg/m3"));

      if (warmUpSteps != 0)
      {
        newIolet->SetWarmup(warmUpSteps);
      }

      return newIolet;
    }

    lb::stents::StentConstantFlux* SimConfig::DoIOForConstantFluxStent(
        const io::xml::Element& stentEl)
    {
      lb::stents::StentConstantFlux* newStent = new lb::stents::StentConstantFlux();

      const io::xml::Element conditionEl = stentEl.GetChildOrThrow("condition");

      const io::xml::Element meanEl = conditionEl.GetChildOrThrow("mean");
      newStent->SetFlux(GetDimensionalValueInLatticeUnits<PhysicalSpeed>(meanEl, "m/s"));

      return newStent;
    }

    lb::iolets::InOutLetWomersleyVelocity* SimConfig::DoIOForWomersleyVelocityInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetWomersleyVelocity* newIolet = new lb::iolets::InOutLetWomersleyVelocity();
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      const io::xml::Element radiusEl = conditionEl.GetChildOrThrow("radius");
      newIolet->SetRadius(GetDimensionalValueInLatticeUnits<LatticeDistance>(radiusEl, "m"));

      const io::xml::Element concentrationEl = conditionEl.GetChildOrThrow("concentration");
      newIolet->SetConcentration(GetDimensionalValueInLatticeUnits<PhysicalDensity>(concentrationEl, "kg/m3"));

      const io::xml::Element pgAmpEl = conditionEl.GetChildOrThrow("pressure_gradient_amplitude");
      newIolet->SetPressureGradientAmplitude(GetDimensionalValueInLatticeUnits<
          LatticePressureGradient>(pgAmpEl, "mmHg/m"));

      const io::xml::Element periodEl = conditionEl.GetChildOrThrow("period");
      newIolet->SetPeriod(GetDimensionalValueInLatticeUnits<LatticeTime>(periodEl, "s"));

      const io::xml::Element womNumEl = conditionEl.GetChildOrThrow("womersley_number");
      newIolet->SetWomersleyNumber(GetDimensionalValueInLatticeUnits<Dimensionless>(womNumEl,
                                                                                    "dimensionless"));

      return newIolet;
    }

    lb::iolets::InOutLetFileVelocity* SimConfig::DoIOForFileVelocityInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetFileVelocity* newIolet = new lb::iolets::InOutLetFileVelocity();
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      std::string velocityFilePath = conditionEl.GetChildOrThrow("path").GetAttributeOrThrow("value");

      velocityFilePath = util::NormalizePathRelativeToPath(velocityFilePath, xmlFilePath);
      newIolet->SetFilePath(velocityFilePath);

      const io::xml::Element radiusEl = conditionEl.GetChildOrThrow("radius");
      newIolet->SetRadius(GetDimensionalValueInLatticeUnits<LatticeDistance>(radiusEl, "m"));

      const io::xml::Element concentrationEl = conditionEl.GetChildOrThrow("concentration");
      newIolet->SetConcentration(GetDimensionalValueInLatticeUnits<PhysicalDensity>(concentrationEl, "kg/m3"));

      return newIolet;
    }

    lb::iolets::InOutLetHollowFileVelocity* SimConfig::DoIOForHollowFileVelocityInOutlet(
        const io::xml::Element& ioletEl)
    {
      lb::iolets::InOutLetHollowFileVelocity* newIolet = new lb::iolets::InOutLetHollowFileVelocity();
      DoIOForBaseInOutlet(ioletEl, newIolet);

      const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

      std::string velocityFilePath = conditionEl.GetChildOrThrow("path").GetAttributeOrThrow("value");

      velocityFilePath = util::NormalizePathRelativeToPath(velocityFilePath, xmlFilePath);
      newIolet->SetFilePath(velocityFilePath);

      const io::xml::Element radiusEl = conditionEl.GetChildOrThrow("radius");
      newIolet->SetRadius(GetDimensionalValueInLatticeUnits<LatticeDistance>(radiusEl, "m"));

      const io::xml::Element innerRadiusEl = conditionEl.GetChildOrThrow("innerradius");
      newIolet->SetInnerRadius(GetDimensionalValueInLatticeUnits<LatticeDistance>(innerRadiusEl, "m"));

      const io::xml::Element concentrationEl = conditionEl.GetChildOrThrow("concentration");
      newIolet->SetConcentration(GetDimensionalValueInLatticeUnits<PhysicalDensity>(concentrationEl, "kg/m3"));

      return newIolet;
    }

    lb::stents::StentFileFlux* SimConfig::DoIOForFileFluxStent(
        const io::xml::Element& stentEl)
    {
      lb::stents::StentFileFlux* newStent = new lb::stents::StentFileFlux();

      const io::xml::Element conditionEl = stentEl.GetChildOrThrow("condition");

      std::string fluxFilePath = conditionEl.GetChildOrThrow("path").GetAttributeOrThrow("value");

      fluxFilePath = util::NormalizePathRelativeToPath(fluxFilePath, xmlFilePath);
      newStent->SetFilePath(fluxFilePath);

      return newStent;
    }

    bool SimConfig::HasColloidSection() const
    {
      return hasColloidSection;
    }

    PhysicalPressure SimConfig::GetInitialPressure() const
    {
      return initialPressure_mmHg;
    }

    PhysicalDensity SimConfig::GetInitialDensity() const
    {
      return initialDensity_Kg_per_m3;
    }

    const util::UnitConverter& SimConfig::GetUnitConverter() const
    {
      if (unitConverter == NULL)
        throw Exception() << "Invalid UnitConverter (NULL)";

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
        throw Exception() << "At least one convergence criterion must be provided in " << convEl.GetPath();
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
      monitoringConfig.convergenceReferenceValue =
          GetDimensionalValueInLatticeUnits<LatticeSpeed>(criterionEl, "m/s");
    }

    const SimConfig::MonitoringConfig* SimConfig::GetMonitoringConfiguration() const
    {
      return &monitoringConfig;
    }
  }
}
