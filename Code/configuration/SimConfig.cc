// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
    SimConfig::SimConfig() :
      xmlFilePath(""), hasColloidSection(false), warmUpSteps(0), rawXmlDoc(NULL)
    {
    }

    SimConfig::SimConfig(const std::string& path) :
      xmlFilePath(path), hasColloidSection(false), warmUpSteps(0), rawXmlDoc(NULL)
    {
      if (!util::file_exists(path.c_str()))
      {
        throw Exception() << "Config file '" << path << "' does not exist";
      }
      rawXmlDoc = new io::xml::Document(path);
      colloidConfigPath = path;
      DoIO(rawXmlDoc->GetRoot());
      dataFilePath = util::NormalizePathRelativeToPath(dataFilePath, path);
    }

    SimConfig::~SimConfig()
    {
      for (unsigned outputNumber = 0; outputNumber < propertyOutputs.size(); ++outputNumber)
      {
        delete propertyOutputs[outputNumber];
      }

      delete rawXmlDoc;
      rawXmlDoc = NULL;
    }

    void SimConfig::DoIO(io::xml::Element topNode)
    {
      // Top element must be:
      // <hemelbsettings version="2" />
      if (topNode.GetName() != "hemelbsettings")
        throw Exception() << "Invalid root element: " << topNode.GetPath();

      unsigned version;
      const std::string& versionStr = topNode.GetAttributeOrThrow("version", version);
      if (version != 2U)
        throw Exception() << "Unrecognised XML version. Expected 2, got " << versionStr;

      DoIOForSimulation(topNode.GetChildOrThrow("simulation"));

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
    }

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
      GetDimensionalValue(tsEl, "s", timeStepLength);

      // Optional element
      // <extra_warmup_steps value="unsigned" units="lattice" />
      const io::xml::Element wuEl = simEl.GetChildOrNull("extra_warmup_steps");
      if (wuEl != io::xml::Element::Missing())
      {
        GetDimensionalValue(wuEl, "lattice", warmUpSteps);
        totalTimeSteps += warmUpSteps;
      }
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
      // TODO: peek into the GMY to get the space step.
    }

    std::vector<lb::iolets::InOutLet*> SimConfig::DoIOForInOutlets(const io::xml::Element ioletsEl)
    {
      const std::string& nodeName = ioletsEl.GetName();

      const std::string childNodeName = nodeName.substr(0, nodeName.size() - 1);
      std::vector<lb::iolets::InOutLet*> ioletList;
      for (io::xml::Element currentIoletNode = ioletsEl.GetChildOrNull(childNodeName); currentIoletNode
          != io::xml::Element::Missing(); currentIoletNode
          = currentIoletNode.NextSiblingOrNull(childNodeName))
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

    lb::iolets::InOutLet* SimConfig::DoIOForPressureInOutlet(const io::xml::Element& ioletEl)
    {
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

    lb::iolets::InOutLet* SimConfig::DoIOForVelocityInOutlet(const io::xml::Element& ioletEl)
    {
      io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");
      const std::string& conditionSubtype = conditionEl.GetAttributeOrThrow("subtype");

      lb::iolets::InOutLet* newIolet = NULL;
      if (conditionSubtype == "parabolic")
      {
        newIolet = DoIOForParabolicVelocityInOutlet(ioletEl);
      }
      else if (conditionSubtype == "womersley")
      {
        newIolet = DoIOForWomersleyVelocityInOutlet(ioletEl);
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
      GetDimensionalValue(orientationEl.GetChildOrThrow("longitude"), "deg", visualisationLongitude);
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
    for (io::xml::ChildIterator poPtr = propertiesEl.IterChildren("propertyoutput"); !poPtr.AtEnd(); ++poPtr)
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

    for (io::xml::ChildIterator fieldPtr = propertyoutputEl.IterChildren("field"); !fieldPtr.AtEnd(); ++fieldPtr)
    file->fields.push_back(DoIOForPropertyField(*fieldPtr));

    return file;
  }

  extraction::StraightLineGeometrySelector* SimConfig::DoIOForLineGeometry(
      const io::xml::Element& geometryEl)
  {
    io::xml::Element point1El = geometryEl.GetChildOrThrow("point");
    io::xml::Element point2El = point1El.NextSiblingOrThrow("point");

    util::Vector3D<float> point1;
    util::Vector3D<float> point2;

    GetDimensionalValue(point1El, "m", point1);
    GetDimensionalValue(point2El, "m", point2);

    return new extraction::StraightLineGeometrySelector(point1, point2);
  }

  extraction::PlaneGeometrySelector* SimConfig::DoIOForPlaneGeometry(
      const io::xml::Element& geometryEl)
  {
    io::xml::Element pointEl = geometryEl.GetChildOrThrow("point");
    io::xml::Element normalEl = geometryEl.GetChildOrThrow("normal");

    util::Vector3D<float> point;
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
      float radius;
      GetDimensionalValue(radiusEl, "m", radius);
      return new extraction::PlaneGeometrySelector(point, normal, radius);
    }

  }

  extraction::SurfacePointSelector* SimConfig::DoIOForSurfacePoint(
      const io::xml::Element& geometryEl)
  {
    io::xml::Element pointEl = geometryEl.GetChildOrThrow("point");

    util::Vector3D<float> point;
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

    util::Vector3D<double> temp;
    GetDimensionalValue(positionEl, "m", temp);
    value->SetPosition(temp);

    GetDimensionalValue(normalEl, "dimensionless", temp);
    value->SetNormal(temp);
  }

  void SimConfig::DoIOForInitialConditions(io::xml::Element initialconditionsEl)
  {
    //, isLoading, initialPressure
    io::xml::Element pressureEl = initialconditionsEl.GetChildOrThrow("pressure");
    io::xml::Element uniformEl = pressureEl.GetChildOrThrow("uniform");

    GetDimensionalValue(uniformEl, "mmHg", initialPressure);
  }

  lb::iolets::InOutLetCosine* SimConfig::DoIOForCosinePressureInOutlet(
      const io::xml::Element& ioletEl)
  {
    lb::iolets::InOutLetCosine* newIolet = new lb::iolets::InOutLetCosine();
    DoIOForBaseInOutlet(ioletEl, newIolet);

    const io::xml::Element conditionEl = ioletEl.GetChildOrThrow("condition");

    GetDimensionalValue(conditionEl.GetChildOrThrow("amplitude"),
        "mmHg",
        newIolet->GetPressureAmp());
    GetDimensionalValue(conditionEl.GetChildOrThrow("mean"), "mmHg", newIolet->GetPressureMean());
    GetDimensionalValue(conditionEl.GetChildOrThrow("phase"), "dimensionless", newIolet->GetPhase());
    GetDimensionalValue(conditionEl.GetChildOrThrow("period"), "s", newIolet->GetPeriod());
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
    newIolet->GetFilePath() = pathEl.GetAttributeOrThrow("value");

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
    GetDimensionalValue(radiusEl, "lattice", newIolet->GetRadius());

    const io::xml::Element maximumEl = conditionEl.GetChildOrThrow("maximum");
    GetDimensionalValue(maximumEl, "lattice", newIolet->GetMaxSpeed());

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
    GetDimensionalValue(radiusEl, "lattice", newIolet->GetRadius());

    const io::xml::Element pgAmpEl = conditionEl.GetChildOrThrow("pressure_gradient_amplitude");
    GetDimensionalValue(pgAmpEl, "lattice", newIolet->GetPressureGradientAmplitude());

    const io::xml::Element periodEl = conditionEl.GetChildOrThrow("period");
    GetDimensionalValue(periodEl, "lattice", newIolet->GetPeriod());

    const io::xml::Element womNumEl = conditionEl.GetChildOrThrow("womersley_number");
    GetDimensionalValue(womNumEl, "dimensionless", newIolet->GetWomersleyNumber());

    return newIolet;
  }

  bool SimConfig::HasColloidSection() const
  {
    return hasColloidSection;
  }

  PhysicalPressure SimConfig::GetInitialPressure() const
  {
    return initialPressure;
  }
}
}
