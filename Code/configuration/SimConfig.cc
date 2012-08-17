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
    // TODO really, the functions that take a TiXML object, a bool, a string and some basic type
    // should go into some Tiny XML abstraction layer. This would remove the need for the IOLet
    // types to know about the SimConfig object (and get rid of a circular dependency).

    SimConfig::SimConfig() :
        LEGACY_PULSATILE_PERIOD(60.0 / 70.0)
    {
      // This constructor only exists to prevent instantiation without
      // using the static load method.
    }

    SimConfig::~SimConfig()
    {
      for (unsigned outputNumber = 0; outputNumber < propertyOutputs.size(); ++outputNumber)
      {
        delete propertyOutputs[outputNumber];
      }
    }

    SimConfig *SimConfig::Load(const char *path)
    {
      util::check_file(path);
      TiXmlDocument configFile;
      configFile.LoadFile(path);

      SimConfig *result = new SimConfig();
      result->colloidConfigPath = path;
      result->DoIO(configFile.FirstChildElement(), true);
      result->dataFilePath = util::NormalizePathRelativeToPath(result->dataFilePath, path);

      return result;
    }

    void SimConfig::Save(std::string path)
    {
      TiXmlDocument configFile;
      TiXmlDeclaration * declaration = new TiXmlDeclaration("1.0", "", "");
      TiXmlElement *topElement = new TiXmlElement("hemelbsettings");

      configFile.LinkEndChild(declaration);
      configFile.LinkEndChild(topElement);

      DoIO(topElement, false);

      configFile.SaveFile(path);
    }

    void SimConfig::DoIO(TiXmlElement *topNode, bool isLoading)
    {
      TiXmlElement* simulationElement = GetChild(topNode, "simulation", isLoading);

      // Backwards compatibility with 0.2.0 input file.
      unsigned long numCycles;
      long stepsPerCycle;
      // short cut means DOIO will not be called for legacy values when saving.
      if (isLoading && DoIOForULong(simulationElement, "cycles", isLoading, numCycles)
          && DoIOForLong(simulationElement, "cyclesteps", isLoading, stepsPerCycle))
      {
        totalTimeSteps = numCycles * stepsPerCycle;
      }
      else
      {
        DoIOForULong(simulationElement, "steps", isLoading, totalTimeSteps);
      }

      if ( (!DoIOForDouble(simulationElement, "step_length", isLoading, timeStepLength)) && isLoading)
      {
        timeStepLength = LEGACY_PULSATILE_PERIOD / stepsPerCycle;
      }

      DoIOForStressType(simulationElement, "stresstype", isLoading, stressType);

      TiXmlElement* geometryElement = GetChild(topNode, "geometry", isLoading);
      if (geometryElement != NULL)
      {
        DoIOForString(GetChild(geometryElement, "datafile", isLoading), "path", isLoading, dataFilePath);
      }

      DoIOForInOutlets(GetChild(topNode, "inlets", isLoading), isLoading, inlets, "inlet");
      DoIOForInOutlets(GetChild(topNode, "outlets", isLoading), isLoading, outlets, "outlet");

      TiXmlElement* visualisationElement = GetChild(topNode, "visualisation", isLoading);
      DoIOForFloatVector(GetChild(visualisationElement, "centre", isLoading), isLoading, visualisationCentre);
      TiXmlElement *lOrientationElement = GetChild(visualisationElement, "orientation", isLoading);
      DoIOForFloat(lOrientationElement, "longitude", isLoading, visualisationLongitude);
      DoIOForFloat(lOrientationElement, "latitude", isLoading, visualisationLatitude);

      TiXmlElement *displayElement = GetChild(visualisationElement, "display", isLoading);

      DoIOForFloat(displayElement, "zoom", isLoading, visualisationZoom);
      DoIOForFloat(displayElement, "brightness", isLoading, visualisationBrightness);

      TiXmlElement *rangeElement = GetChild(visualisationElement, "range", isLoading);

      DoIOForFloat(rangeElement, "maxvelocity", isLoading, maxVelocity);
      DoIOForFloat(rangeElement, "maxstress", isLoading, maxStress);

      TiXmlElement* propertyExtractionElement = GetChild(topNode, "properties", isLoading);
      DoIOForProperties(propertyExtractionElement, isLoading);
    }

    void SimConfig::DoIOForFloat(TiXmlElement* parent, std::string attributeName, bool isLoading, float &value)
    {
      if (isLoading)
      {
        char *dummy;
        value = (float) std::strtod(parent->Attribute(attributeName)->c_str(), &dummy);
      }
      else
      {
        // This should be ample.
        std::stringstream output(std::stringstream::out);

        output.precision(6);

        // %g uses the shorter of decimal / mantissa-exponent notations.
        // 6 significant figures will be written.
        output << value;

        parent->SetAttribute(attributeName, output.str());
      }
    }

    bool SimConfig::DoIOForDouble(TiXmlElement* parent, std::string attributeName, bool isLoading, double &value)
    {
      if (isLoading)
      {
        const std::string* data = parent->Attribute(attributeName);

        if (data != NULL)
        {
          char *dummy;
          value = std::strtod(parent->Attribute(attributeName)->c_str(), &dummy);
          return true;
        }
        else
        {
          value = 0.0;
          return false;
        }
      }
      else
      {
        // This should be ample.
        std::stringstream output(std::stringstream::out);

        output.precision(6);

        // %g uses the shorter of decimal / mantissa-exponent notations.
        // 6 significant figures will be written.
        output << value;

        parent->SetAttribute(attributeName, output.str());
        return true;
      }
    }

    void SimConfig::DoIOForString(TiXmlElement* parent, std::string attributeName, bool isLoading, std::string &value)
    {
      if (isLoading)
      {
        // Compare to 0 not NULL, because that's what Attribute in TinyXml returns
        if (parent->Attribute(attributeName) == 0)
        {
          value = "";
        }
        else
        {
          value = std::string(parent->Attribute(attributeName)->c_str());
        }
      }
      else
      {
        if (value != "")
        {
          parent->SetAttribute(attributeName, value);
        }
      }
    }

    bool SimConfig::DoIOForLong(TiXmlElement* parent, std::string attributeName, bool isLoading, long &value)
    {
      if (isLoading)
      {
        const std::string *read_result = parent->Attribute(attributeName);
        if (read_result == NULL)
        {
          value = 0.0;
          return false;
        }
        char *dummy;
        // Read in, in base 10.
        value = std::strtol(read_result->c_str(), &dummy, 10);
        return true;
      }
      else
      {
        // This should be ample.
        std::stringstream output(std::stringstream::out);

        output << value;

        parent->SetAttribute(attributeName, output.str());
        return true;
      }
    }

    bool SimConfig::DoIOForULong(TiXmlElement* parent, std::string attributeName, bool isLoading, unsigned long &value)
    {
      if (isLoading)
      {
        const std::string *read_result = parent->Attribute(attributeName);
        if (read_result == NULL)
        {
          value = 0.0;
          return false;
        }
        char *dummy;
        // Read in, in base 10.
        value = std::strtoul(read_result->c_str(), &dummy, 10);
        return true;
      }
      else
      {
        // This should be ample.
        std::stringstream output(std::stringstream::out);

        output << value;

        parent->SetAttribute(attributeName, output.str());
        return true;
      }
    }

    void SimConfig::DoIOForStressType(TiXmlElement* xmlNode,
                                      std::string attributeName,
                                      bool isLoading,
                                      lb::StressTypes &value)
    {
      if (isLoading)
      {
        char *dummy;
        // Read in, in base 10.
        value = (lb::StressTypes) std::strtol(xmlNode->Attribute(attributeName)->c_str(), &dummy, 10);
      }
      else
      {
        std::stringstream output(std::stringstream::out);

        output << (int) value;

        xmlNode->SetAttribute(attributeName, output.str());
      }
    }

    void SimConfig::DoIOForInOutlets(TiXmlElement *parent,
                                     bool isLoading,
                                     std::vector<lb::boundaries::iolets::InOutLet*> &bResult,
                                     std::string childNodeName)
    {
      if (isLoading)
      {
        TiXmlElement *currentIoletNode = parent->FirstChildElement(childNodeName);

        while (currentIoletNode != NULL)
        {
          // Determine which InOutlet to create
          // This is done by checking if a path is specified
          std::string PFilePath;
          std::string MultiscaleLabel;
          DoIOForString(GetChild(GetChild(parent, childNodeName, isLoading), "pressure", isLoading),
                        "path",
                        isLoading,
                        PFilePath);
          DoIOForString(GetChild(GetChild(parent, childNodeName, isLoading), "pressure", isLoading),
                        "label",
                        isLoading,
                        MultiscaleLabel);
          lb::boundaries::iolets::InOutLet *newIolet;
          if (PFilePath != "")
          {
            // If there is a file specified we use it
            newIolet = new lb::boundaries::iolets::InOutLetFile();

          }
          else if (MultiscaleLabel != "")
          {
            newIolet = new lb::boundaries::iolets::InOutLetVelocityAware();
          }
          else
          {
            // If no file is specified we use a cosine trace
            newIolet = new lb::boundaries::iolets::InOutLetCosine();
          }

          newIolet->DoIO(GetChild(parent, childNodeName, isLoading), isLoading, this);
          bResult.push_back(newIolet);
          currentIoletNode = currentIoletNode->NextSiblingElement(childNodeName);
        }
      }
      else
      {
        for (unsigned int ii = 0; ii < bResult.size(); ii++)
        {
          // NB we're good up to 99 io-lets here.
          bResult[ii]->DoIO(GetChild(parent, childNodeName, isLoading), isLoading, this);
        }
      }
    }

    void SimConfig::DoIOForProperties(TiXmlElement *xmlNode, bool isLoading)
    {
      if (isLoading && xmlNode == NULL)
      {
        return;
      }

      TiXmlElement *currentPropertyNode = xmlNode->FirstChildElement();

      while (currentPropertyNode != NULL)
      {
        extraction::PropertyOutputFile* output = new extraction::PropertyOutputFile();

        DoIOForPropertyOutputFile(currentPropertyNode, isLoading, output);
        propertyOutputs.push_back(output);

        currentPropertyNode = currentPropertyNode->NextSiblingElement();
      }
    }

    void SimConfig::DoIOForPropertyOutputFile(TiXmlElement *xmlNode,
                                              bool isLoading,
                                              extraction::PropertyOutputFile* file)
    {
      if (isLoading)
      {
        DoIOForString(xmlNode, "file", isLoading, file->filename);

        char* dummy;
        file->frequency = strtoul(xmlNode->Attribute("frequency"), &dummy, 10);

        TiXmlElement* propertyElement = xmlNode->FirstChildElement();

        if (isLoading)
        {
          if (propertyElement->ValueStr().compare("planegeometry") == 0)
          {
            extraction::PlaneGeometrySelector* plane = NULL;
            DoIOForPlaneGeometry(propertyElement, isLoading, plane);
            file->geometry = plane;
          }
          else if (propertyElement->ValueStr().compare("linegeometry") == 0)
          {
            extraction::StraightLineGeometrySelector* line = NULL;
            DoIOForLineGeometry(propertyElement, isLoading, line);
            file->geometry = line;
          }
          else if (propertyElement->ValueStr().compare("wholegeometry") == 0)
          {
            extraction::WholeGeometrySelector* whole = new extraction::WholeGeometrySelector();
            file->geometry = whole;
          }
          else if (propertyElement->ValueStr().compare("geometrysurface") == 0)
          {
            extraction::GeometrySurfaceSelector* surface = new extraction::GeometrySurfaceSelector();
            file->geometry = surface;
          }
          else
          {
            log::Logger::Log<log::Info, log::OnePerCore>("Unrecognised geometry type: %s", xmlNode->Value());
            exit(1);
          }

          TiXmlElement *currentFieldNode = propertyElement->FirstChildElement("field");

          while (currentFieldNode != NULL)
          {
            extraction::OutputField outputField;

            DoIOForPropertyField(currentFieldNode, isLoading, outputField);
            file->fields.push_back(outputField);

            currentFieldNode = currentFieldNode->NextSiblingElement("field");
          }
        }
      }
    }

    void SimConfig::DoIOForLineGeometry(TiXmlElement *xmlNode,
                                        bool isLoading,
                                        extraction::StraightLineGeometrySelector*& line)
    {
      TiXmlElement* point1 = GetChild(xmlNode, "point", isLoading);
      TiXmlElement* point2 = isLoading ?
        point1->NextSiblingElement("point") :
        GetChild(xmlNode, "point", isLoading);

      util::Vector3D<float> mutableVector;
      util::Vector3D<float> mutableVector2;

      if (isLoading)
      {
        DoIOForFloatVector(point1, isLoading, mutableVector);
        DoIOForFloatVector(point2, isLoading, mutableVector2);

        line = new extraction::StraightLineGeometrySelector(mutableVector, mutableVector2);
      }
      else
      {
        mutableVector = line->GetEndpoint1();
        mutableVector2 = line->GetEndpoint2();

        DoIOForFloatVector(point1, isLoading, mutableVector);
        DoIOForFloatVector(point2, isLoading, mutableVector2);
      }
    }

    void SimConfig::DoIOForPlaneGeometry(TiXmlElement *xmlNode,
                                         bool isLoading,
                                         extraction::PlaneGeometrySelector*& plane)
    {
      TiXmlElement* point1 = GetChild(xmlNode, "point", isLoading);
      TiXmlElement* normal = GetChild(xmlNode, "normal", isLoading);

      util::Vector3D<float> mutableVector;
      util::Vector3D<float> mutableVector2;

      if (isLoading)
      {
        DoIOForFloatVector(point1, isLoading, mutableVector);
        DoIOForFloatVector(normal, isLoading, mutableVector2);

        const char* radiusAttribute = xmlNode->Attribute("radius");

        if (radiusAttribute == NULL)
        {
          plane = new extraction::PlaneGeometrySelector(mutableVector, mutableVector2);
        }
        else
        {
          char* dummy;
          float radius = (float) std::strtod(radiusAttribute, &dummy);
          plane = new extraction::PlaneGeometrySelector(mutableVector, mutableVector2, radius);
        }
      }
      else
      {
        mutableVector = plane->GetPoint();
        mutableVector = plane->GetNormal();
        float radius = plane->GetRadius();

        DoIOForFloatVector(point1, isLoading, mutableVector);
        DoIOForFloatVector(normal, isLoading, mutableVector2);
        xmlNode->SetAttribute("radius", radius);
      }
    }

    void SimConfig::DoIOForPropertyField(TiXmlElement *xmlNode, bool isLoading, extraction::OutputField& field)
    {
      std::string type;
      DoIOForString(xmlNode, "type", isLoading, type);
      DoIOForString(xmlNode, "name", isLoading, field.name);

      // Default name is identical to type.
      if (isLoading && field.name.length() == 0)
      {
        field.name = type;
      }

      if (type.compare("pressure") == 0)
      {
        field.type = extraction::OutputField::Pressure;
      }
      else if (type.compare("velocity") == 0)
      {
        field.type = extraction::OutputField::Velocity;
      }
      else if (type.compare("vonmisesstress") == 0)
      {
        field.type = extraction::OutputField::VonMisesStress;
      }
      else if (type.compare("shearstress") == 0)
      {
        field.type = extraction::OutputField::ShearStress;
      }
      else if (type.compare("shearrate") == 0)
      {
        field.type = extraction::OutputField::ShearRate;
      }
      else
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Unrecognised field type (%s) in xml file", type.c_str());
        exit(1);
      }
    }

    void SimConfig::DoIOForBaseInOutlet(TiXmlElement *parent,
                                        bool isLoading,
                                        lb::boundaries::iolets::InOutLet* const value)
    {
      TiXmlElement* lPositionElement = GetChild(parent, "position", isLoading);
      TiXmlElement* lNormalElement = GetChild(parent, "normal", isLoading);
      DoIOForFloatVector(lPositionElement, isLoading, value->GetPosition());
      DoIOForFloatVector(lNormalElement, isLoading, value->GetNormal());
    }

    void SimConfig::DoIOForCosineInOutlet(TiXmlElement *parent,
                                          bool isLoading,
                                          lb::boundaries::iolets::InOutLetCosine* const value)
    {
      DoIOForBaseInOutlet(parent,isLoading,value);
      TiXmlElement* lPressureElement = GetChild(parent, "pressure", isLoading);

      DoIOForDouble(lPressureElement, "mean", isLoading, value->GetPressureMean());
      DoIOForDouble(lPressureElement, "amplitude", isLoading, value->GetPressureAmp());
      DoIOForDouble(lPressureElement, "phase", isLoading, value->GetPhase());

      if (!DoIOForDouble(lPressureElement, "period", isLoading, value->GetPeriod()) && isLoading)
      {
        value->SetPeriod(LEGACY_PULSATILE_PERIOD);
      }
    }

    void SimConfig::DoIOForFileInOutlet(TiXmlElement *parent,
                                        bool isLoading,
                                        lb::boundaries::iolets::InOutLetFile* const value)
    {
      DoIOForBaseInOutlet(parent,isLoading,value);

      TiXmlElement* lPressureElement = GetChild(parent, "pressure", isLoading);

      DoIOForString(lPressureElement, "path", isLoading, value->GetFilePath());

    }

    void SimConfig::DoIOForMultiscaleInOutlet(TiXmlElement *parent,
                                              bool isLoading,
                                              lb::boundaries::iolets::InOutLetMultiscale* const value)
    {
      DoIOForBaseInOutlet(parent,isLoading,value);

      TiXmlElement* lPressureElement = GetChild(parent, "pressure", isLoading);
      DoIOForDouble(lPressureElement,"pressure",isLoading,value->GetPressureReference());
      DoIOForDouble(lPressureElement,"velocity",isLoading,value->GetVelocityReference());
      DoIOForString(lPressureElement, "label", isLoading, value->GetLabel());
    }

    void SimConfig::DoIOForFloatVector(TiXmlElement *parent, bool isLoading, util::Vector3D<float> &value)
    {
      DoIOForFloat(parent, "x", isLoading, value.x);
      DoIOForFloat(parent, "y", isLoading, value.y);
      DoIOForFloat(parent, "z", isLoading, value.z);
    }

    TiXmlElement *SimConfig::GetChild(TiXmlElement *parent, std::string childNodeName, bool isLoading)
    {
      if (isLoading)
      {
        return parent->FirstChildElement(childNodeName);
      }
      else
      {
        TiXmlElement* newChild = new TiXmlElement(childNodeName);
        parent->LinkEndChild(newChild);
        return newChild;
      }
    }
  }
}
