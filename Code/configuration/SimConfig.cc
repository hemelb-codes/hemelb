#include "configuration/SimConfig.h"

#include <string>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>

#include "util/fileutils.h"
#include "debug/Debugger.h"

namespace hemelb
{
  namespace configuration
  {
    // TODO really, the functions that take a TiXML object, a bool, a string and some basic type
    // should go into some Tiny XML abstraction layer. This would remove the need for the IOLet
    // types to know about the SimConfig object (and get rid of a circular dependency).

    SimConfig::SimConfig()
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

    SimConfig *SimConfig::Load(const char *iPath)
    {
      util::check_file(iPath);
      TiXmlDocument *lConfigFile = new TiXmlDocument();
      lConfigFile->LoadFile(iPath);

      SimConfig *lRet = new SimConfig();
      lRet->DoIO(lConfigFile->FirstChildElement(), true);
      lRet->DataFilePath = util::NormalizePathRelativeToPath(lRet->DataFilePath, iPath);
      delete lConfigFile;

      return lRet;
    }

    void SimConfig::Save(std::string iPath)
    {
      TiXmlDocument lConfigFile;
      TiXmlDeclaration * lDeclaration = new TiXmlDeclaration("1.0", "", "");
      TiXmlElement *lTopElement = new TiXmlElement("hemelbsettings");

      lConfigFile.LinkEndChild(lDeclaration);
      lConfigFile.LinkEndChild(lTopElement);

      DoIO(lTopElement, false);

      lConfigFile.SaveFile(iPath);
    }

    void SimConfig::DoIO(TiXmlElement *topNode, bool isLoading)
    {
      TiXmlElement* simulationElement = GetChild(topNode, "simulation", isLoading);

      // Backwards compatibility with 0.2.0 input file.
      unsigned long NumCycles;
      long StepsPerCycle;
      // short cut means DOIO will not be called for legacy values when saving.
      if (isLoading
          && DoIOForULong(simulationElement, "cycles", isLoading, NumCycles)
          && DoIOForLong(simulationElement, "cyclesteps", isLoading, StepsPerCycle))
      {
        TotalTimeSteps=NumCycles*StepsPerCycle;
      }
      else
      {
        DoIOForULong(simulationElement, "steps", isLoading, TotalTimeSteps);
      }

      if ((!DoIOForDouble(simulationElement, "step_length",isLoading,TimeStepLength)) && isLoading)
      {
        TimeStepLength=60.0/(70.0*StepsPerCycle);
      }

      DoIOForStressType(simulationElement, "stresstype", isLoading, StressType);

      TiXmlElement* geometryElement = GetChild(topNode, "geometry", isLoading);
      if (geometryElement != NULL)
      {
        DoIOForString(GetChild(geometryElement, "datafile", isLoading), "path", isLoading, DataFilePath);
      }

      DoIOForInOutlets(GetChild(topNode, "inlets", isLoading), isLoading, Inlets, "inlet");
      DoIOForInOutlets(GetChild(topNode, "outlets", isLoading), isLoading, Outlets, "outlet");

      TiXmlElement* visualisationElement = GetChild(topNode, "visualisation", isLoading);
      DoIOForFloatVector(GetChild(visualisationElement, "centre", isLoading), isLoading, VisCentre);
      TiXmlElement *lOrientationElement = GetChild(visualisationElement, "orientation", isLoading);
      DoIOForFloat(lOrientationElement, "longitude", isLoading, VisLongitude);
      DoIOForFloat(lOrientationElement, "latitude", isLoading, VisLatitude);

      TiXmlElement *displayElement = GetChild(visualisationElement, "display", isLoading);

      DoIOForFloat(displayElement, "zoom", isLoading, VisZoom);
      DoIOForFloat(displayElement, "brightness", isLoading, VisBrightness);

      TiXmlElement *rangeElement = GetChild(visualisationElement, "range", isLoading);

      DoIOForFloat(rangeElement, "maxvelocity", isLoading, MaxVelocity);
      DoIOForFloat(rangeElement, "maxstress", isLoading, MaxStress);

      TiXmlElement* propertyExtractionElement = GetChild(topNode, "properties", isLoading);
      DoIOForProperties(propertyExtractionElement, isLoading);
    }

    void SimConfig::DoIOForFloat(TiXmlElement* iParent, std::string iAttributeName, bool iIsLoading, float &value)
    {
      if (iIsLoading)
      {
        char *dummy;
        value = (float) std::strtod(iParent->Attribute(iAttributeName)->c_str(), &dummy);
      }
      else
      {
        // This should be ample.
        std::stringstream output(std::stringstream::out);

        output.precision(6);

        // %g uses the shorter of decimal / mantissa-exponent notations.
        // 6 significant figures will be written.
        output << value;

        iParent->SetAttribute(iAttributeName, output.str());
      }
    }

    bool SimConfig::DoIOForDouble(TiXmlElement* iParent, std::string iAttributeName, bool iIsLoading, double &value)
    {
      if (iIsLoading)
      {
        const std::string* data = iParent->Attribute(iAttributeName);

        if (data != NULL)
        {
          char *dummy;
          value = std::strtod(iParent->Attribute(iAttributeName)->c_str(), &dummy);
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

        iParent->SetAttribute(iAttributeName, output.str());
        return true;
      }
    }

    void SimConfig::DoIOForString(TiXmlElement* iParent,
                                  std::string iAttributeName,
                                  bool iIsLoading,
                                  std::string &iValue)
    {
      if (iIsLoading)
      {
        // Compare to 0 not NULL, because that's what Attribute in TinyXml returns
        if (iParent->Attribute(iAttributeName) == 0)
        {
          iValue = "";
        }
        else
        {
          iValue = std::string(iParent->Attribute(iAttributeName)->c_str());
        }
      }
      else
      {
        if (iValue != "")
        {
          iParent->SetAttribute(iAttributeName, iValue);
        }
      }
    }

    bool SimConfig::DoIOForLong(TiXmlElement* iParent, std::string iAttributeName, bool iIsLoading, long &bValue)
    {
      if (iIsLoading)
      {
        const std::string *read_result = iParent->Attribute(iAttributeName);
        if (read_result == NULL)
        {
          bValue = 0.0;
          return false;
        }
        char *dummy;
        // Read in, in base 10.
        bValue = std::strtol(read_result->c_str(), &dummy, 10);
        return true;
      }
      else
      {
        // This should be ample.
        std::stringstream output(std::stringstream::out);

        output << bValue;

        iParent->SetAttribute(iAttributeName, output.str());
        return true;
      }
    }

    bool SimConfig::DoIOForULong(TiXmlElement* iParent,
                                 std::string iAttributeName,
                                 bool iIsLoading,
                                 unsigned long &bValue)
    {
      if (iIsLoading)
      {
        const std::string *read_result = iParent->Attribute(iAttributeName);
        if (read_result == NULL)
        {
          bValue = 0.0;
          return false;
        }
        char *dummy;
        // Read in, in base 10.
        bValue = std::strtoul(read_result->c_str(), &dummy, 10);
        return true;
      }
      else
      {
        // This should be ample.
        std::stringstream output(std::stringstream::out);

        output << bValue;

        iParent->SetAttribute(iAttributeName, output.str());
        return true;
      }
    }

    void SimConfig::DoIOForStressType(TiXmlElement* iXmlNode,
                                      std::string iAttributeName,
                                      bool iIsLoading,
                                      lb::StressTypes &value)
    {
      if (iIsLoading)
      {
        char *dummy;
        // Read in, in base 10.
        value = (lb::StressTypes) std::strtol(iXmlNode->Attribute(iAttributeName)->c_str(), &dummy, 10);
      }
      else
      {
        std::stringstream output(std::stringstream::out);

        output << (int) value;

        iXmlNode->SetAttribute(iAttributeName, output.str());
      }
    }

    void SimConfig::DoIOForInOutlets(TiXmlElement *iParent,
                                     bool iIsLoading,
                                     std::vector<lb::boundaries::iolets::InOutLet*> &bResult,
                                     std::string iChildNodeName)
    {
      if (iIsLoading)
      {
        TiXmlElement *lCurrentLet = iParent->FirstChildElement(iChildNodeName);

        while (lCurrentLet != NULL)
        {
          // Determine which InOutlet to create
          // This is done by checking if a path is specified
          std::string PFilePath;
          DoIOForString(GetChild(GetChild(iParent, iChildNodeName, iIsLoading), "pressure", iIsLoading),
                        "path",
                        iIsLoading,
                        PFilePath);
          lb::boundaries::iolets::InOutLet *lNew;
          if (PFilePath == "")
          {
            // If no file is specified we use a cosine trace
            lNew = new lb::boundaries::iolets::InOutLetCosine();
          }
          else
          {
            // If there is a file specified we use it
            lNew = new lb::boundaries::iolets::InOutLetFile();
          }

          lNew->DoIO(GetChild(iParent, iChildNodeName, iIsLoading), iIsLoading, this);
          bResult.push_back(lNew);
          lCurrentLet = lCurrentLet->NextSiblingElement(iChildNodeName);
        }
      }
      else
      {
        for (unsigned int ii = 0; ii < bResult.size(); ii++)
        {
          // NB we're good up to 99 io-lets here.
          bResult[ii]->DoIO(GetChild(iParent, iChildNodeName, iIsLoading), iIsLoading, this);
        }
      }
    }

    void SimConfig::DoIOForProperties(TiXmlElement *iXmlNode, bool iIsLoading)
    {
      if (iIsLoading && iXmlNode == NULL)
      {
        return;
      }

      TiXmlElement *lCurrentLet = iXmlNode->FirstChildElement();

      while (lCurrentLet != NULL)
      {
        extraction::PropertyOutputFile* output = new extraction::PropertyOutputFile();

        DoIOForPropertyOutputFile(lCurrentLet, iIsLoading, output);
        propertyOutputs.push_back(output);

        lCurrentLet = lCurrentLet->NextSiblingElement();
      }
    }

    void SimConfig::DoIOForPropertyOutputFile(TiXmlElement *iXmlNode,
                                              bool iIsLoading,
                                              extraction::PropertyOutputFile* file)
    {
      if (iIsLoading)
      {
        DoIOForString(iXmlNode, "file", iIsLoading, file->filename);

        char* dummy;
        file->frequency = strtoul(iXmlNode->Attribute("frequency"), &dummy, 10);

        TiXmlElement* propertyElement = iXmlNode->FirstChildElement();

        if (iIsLoading)
        {
          if (propertyElement->ValueStr().compare("planegeometry") == 0)
          {
            extraction::PlaneGeometrySelector* plane = NULL;
            DoIOForPlaneGeometry(propertyElement, iIsLoading, plane);
            file->geometry = plane;
          }
          else if (propertyElement->ValueStr().compare("linegeometry") == 0)
          {
            extraction::StraightLineGeometrySelector* line = NULL;
            DoIOForLineGeometry(propertyElement, iIsLoading, line);
            file->geometry = line;
          }
          else if (propertyElement->ValueStr().compare("wholegeometry") == 0)
          {
            extraction::WholeGeometrySelector* whole = new extraction::WholeGeometrySelector();
            file->geometry = whole;
          }
          else
          {
            log::Logger::Log<log::Info, log::OnePerCore>("Unrecognised geometry type: %s", iXmlNode->Value());
            exit(1);
          }

          TiXmlElement *lCurrentLet = propertyElement->FirstChildElement("field");

          while (lCurrentLet != NULL)
          {
            extraction::OutputField outputField;

            DoIOForPropertyField(lCurrentLet, iIsLoading, outputField);
            file->fields.push_back(outputField);

            lCurrentLet = lCurrentLet->NextSiblingElement("field");
          }
        }
      }
    }

    void SimConfig::DoIOForLineGeometry(TiXmlElement *iXmlNode,
                                        bool iIsLoading,
                                        extraction::StraightLineGeometrySelector*& line)
    {
      TiXmlElement* point1 = GetChild(iXmlNode, "point", iIsLoading);
      TiXmlElement* point2 = iIsLoading ?
        point1->NextSiblingElement("point") :
        GetChild(iXmlNode, "point", iIsLoading);

      util::Vector3D<float> mutableVector;
      util::Vector3D<float> mutableVector2;

      if (iIsLoading)
      {
        DoIOForFloatVector(point1, iIsLoading, mutableVector);
        DoIOForFloatVector(point2, iIsLoading, mutableVector2);

        line = new extraction::StraightLineGeometrySelector(mutableVector, mutableVector2);
      }
      else
      {
        mutableVector = line->GetEndpoint1();
        mutableVector2 = line->GetEndpoint2();

        DoIOForFloatVector(point1, iIsLoading, mutableVector);
        DoIOForFloatVector(point2, iIsLoading, mutableVector2);
      }
    }

    void SimConfig::DoIOForPlaneGeometry(TiXmlElement *iXmlNode,
                                         bool iIsLoading,
                                         extraction::PlaneGeometrySelector*& plane)
    {
      TiXmlElement* point1 = GetChild(iXmlNode, "point", iIsLoading);
      TiXmlElement* normal = GetChild(iXmlNode, "normal", iIsLoading);

      util::Vector3D<float> mutableVector;
      util::Vector3D<float> mutableVector2;

      if (iIsLoading)
      {
        DoIOForFloatVector(point1, iIsLoading, mutableVector);
        DoIOForFloatVector(normal, iIsLoading, mutableVector2);

        const char* radiusAttribute = iXmlNode->Attribute("radius");

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

        DoIOForFloatVector(point1, iIsLoading, mutableVector);
        DoIOForFloatVector(normal, iIsLoading, mutableVector2);
        iXmlNode->SetAttribute("radius", radius);
      }
    }

    void SimConfig::DoIOForPropertyField(TiXmlElement *iXmlNode, bool iIsLoading, extraction::OutputField& field)
    {
      std::string type;
      DoIOForString(iXmlNode, "type", iIsLoading, type);
      DoIOForString(iXmlNode, "name", iIsLoading, field.name);

      // Default name is identical to type.
      if (iIsLoading && field.name.length() == 0)
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
      else
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Unrecognised field type (%s) in xml file", type.c_str());
        exit(1);
      }
    }

    void SimConfig::DoIOForCosineInOutlet(TiXmlElement *iParent,
                                          bool iIsLoading,
                                          lb::boundaries::iolets::InOutLetCosine* const value)
    {
      TiXmlElement* lPositionElement = GetChild(iParent, "position", iIsLoading);
      TiXmlElement* lNormalElement = GetChild(iParent, "normal", iIsLoading);
      TiXmlElement* lPressureElement = GetChild(iParent, "pressure", iIsLoading);

      DoIOForDouble(lPressureElement, "mean", iIsLoading, value->PressureMeanPhysical);
      DoIOForDouble(lPressureElement, "amplitude", iIsLoading, value->PressureAmpPhysical);
      DoIOForDouble(lPressureElement, "phase", iIsLoading, value->Phase);

      DoIOForFloatVector(lPositionElement, iIsLoading, value->Position);
      DoIOForFloatVector(lNormalElement, iIsLoading, value->Normal);
    }

    void SimConfig::DoIOForFileInOutlet(TiXmlElement *iParent,
                                        bool iIsLoading,
                                        lb::boundaries::iolets::InOutLetFile* const value)
    {
      TiXmlElement* lPositionElement = GetChild(iParent, "position", iIsLoading);
      TiXmlElement* lNormalElement = GetChild(iParent, "normal", iIsLoading);
      TiXmlElement* lPressureElement = GetChild(iParent, "pressure", iIsLoading);

      DoIOForString(lPressureElement, "path", iIsLoading, value->PressureFilePath);

      DoIOForFloatVector(lPositionElement, iIsLoading, value->Position);
      DoIOForFloatVector(lNormalElement, iIsLoading, value->Normal);
    }

    void SimConfig::DoIOForFloatVector(TiXmlElement *iParent, bool iIsLoading, util::Vector3D<float> &iValue)
    {
      DoIOForFloat(iParent, "x", iIsLoading, iValue.x);
      DoIOForFloat(iParent, "y", iIsLoading, iValue.y);
      DoIOForFloat(iParent, "z", iIsLoading, iValue.z);
    }

    TiXmlElement *SimConfig::GetChild(TiXmlElement *iParent, std::string iChildNodeName, bool iIsLoading)
    {
      if (iIsLoading)
      {
        return iParent->FirstChildElement(iChildNodeName);
      }
      else
      {
        TiXmlElement* lNewChild = new TiXmlElement(iChildNodeName);
        iParent->LinkEndChild(lNewChild);
        return lNewChild;
      }
    }
  }
}
