#include "SimConfig.h"

#include <string>
#include <iostream>
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

    void SimConfig::DoIO(TiXmlElement *iTopNode, bool iIsLoading)
    {
      TiXmlElement* lSimulationElement = GetChild(iTopNode, "simulation", iIsLoading);
      DoIO(lSimulationElement, "cycles", iIsLoading, NumCycles);
      DoIO(lSimulationElement, "cyclesteps", iIsLoading, StepsPerCycle);
      long dummyStress = -1;
      DoIO(lSimulationElement, "stresstype", iIsLoading, dummyStress);
      StressType = (lb::StressTypes) dummyStress;

      TiXmlElement* lGeometryElement = GetChild(iTopNode, "geometry", iIsLoading);

      if (lGeometryElement != NULL)
      {
        DoIO(GetChild(lGeometryElement, "datafile", iIsLoading), "path", iIsLoading, DataFilePath);
      }

      DoIO(GetChild(iTopNode, "inlets", iIsLoading), iIsLoading, Inlets, "inlet");

      DoIO(GetChild(iTopNode, "outlets", iIsLoading), iIsLoading, Outlets, "outlet");

      TiXmlElement* lVisualisationElement = GetChild(iTopNode, "visualisation", iIsLoading);
      DoIO(GetChild(lVisualisationElement, "centre", iIsLoading), iIsLoading, VisCentre);
      TiXmlElement *lOrientationElement = GetChild(lVisualisationElement,
                                                   "orientation",
                                                   iIsLoading);
      DoIO(lOrientationElement, "longitude", iIsLoading, VisLongitude);
      DoIO(lOrientationElement, "latitude", iIsLoading, VisLatitude);

      TiXmlElement *lDisplayElement = GetChild(lVisualisationElement, "display", iIsLoading);

      DoIO(lDisplayElement, "zoom", iIsLoading, VisZoom);
      DoIO(lDisplayElement, "brightness", iIsLoading, VisBrightness);

      TiXmlElement *lRangeElement = GetChild(lVisualisationElement, "range", iIsLoading);

      DoIO(lRangeElement, "maxvelocity", iIsLoading, MaxVelocity);
      DoIO(lRangeElement, "maxstress", iIsLoading, MaxStress);
    }

    void SimConfig::DoIO(TiXmlElement* iParent,
                         std::string iAttributeName,
                         bool iIsLoading,
                         float &value)
    {
      if (iIsLoading)
      {
        char *dummy;
        value = (float) std::strtod(iParent->Attribute(iAttributeName)->c_str(), &dummy);
      }
      else
      {
        // This should be ample.
        char lStringValue[20];

        // %g uses the shorter of decimal / mantissa-exponent notations.
        // 6 significant figures will be written.
        sprintf(lStringValue, "%.6g", value);

        iParent->SetAttribute(iAttributeName, lStringValue);
      }
    }

    void SimConfig::DoIO(TiXmlElement* iParent,
                         std::string iAttributeName,
                         bool iIsLoading,
                         double &value)
    {
      if (iIsLoading)
      {
        const std::string* data = iParent->Attribute(iAttributeName);

        if (data != NULL)
        {
          char *dummy;
          value = std::strtod(iParent->Attribute(iAttributeName)->c_str(), &dummy);
        }
        else
        {
          value = 0.0;
        }
      }
      else
      {
        // This should be ample.
        char lStringValue[20];

        // %g uses the shorter of decimal / mantissa-exponent notations.
        // 6 significant figures will be written.
        sprintf(lStringValue, "%.6g", value);

        iParent->SetAttribute(iAttributeName, lStringValue);
      }
    }

    void SimConfig::DoIO(TiXmlElement* iParent,
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

    void SimConfig::DoIO(TiXmlElement* iParent,
                         std::string iAttributeName,
                         bool iIsLoading,
                         long &bValue)
    {
      if (iIsLoading)
      {
        char *dummy;
        // Read in, in base 10.
        bValue = std::strtol(iParent->Attribute(iAttributeName)->c_str(), &dummy, 10);
      }
      else
      {
        // This should be ample.
        char lStringValue[20];

        // %ld specifies long integer style.
        sprintf(lStringValue, "%ld", bValue);

        iParent->SetAttribute(iAttributeName, lStringValue);
      }
    }

    void SimConfig::DoIO(TiXmlElement* iParent,
                         std::string iAttributeName,
                         bool iIsLoading,
                         unsigned long &bValue)
    {
      if (iIsLoading)
      {
        char *dummy;
        // Read in, in base 10.
        bValue = std::strtoul(iParent->Attribute(iAttributeName)->c_str(), &dummy, 10);
      }
      else
      {
        // This should be ample.
        char lStringValue[20];

        // %ld specifies long integer style.
        sprintf(lStringValue, "%ld", bValue);

        iParent->SetAttribute(iAttributeName, lStringValue);
      }
    }

    void SimConfig::DoIO(TiXmlElement *iParent,
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
          DoIO(GetChild(GetChild(iParent, iChildNodeName, iIsLoading), "pressure", iIsLoading),
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

    void SimConfig::DoIO(TiXmlElement *iParent,
                         bool iIsLoading,
                         lb::boundaries::iolets::InOutLetCosine* const value)
    {
      TiXmlElement* lPositionElement = GetChild(iParent, "position", iIsLoading);
      TiXmlElement* lNormalElement = GetChild(iParent, "normal", iIsLoading);
      TiXmlElement* lPressureElement = GetChild(iParent, "pressure", iIsLoading);

      DoIO(lPressureElement, "mean", iIsLoading, value->PressureMeanPhysical);
      DoIO(lPressureElement, "amplitude", iIsLoading, value->PressureAmpPhysical);
      DoIO(lPressureElement, "phase", iIsLoading, value->Phase);
      value->PressureMinPhysical = value->PressureMeanPhysical - value->PressureAmpPhysical;
      value->PressureMaxPhysical = value->PressureMeanPhysical + value->PressureAmpPhysical;

      DoIO(lPositionElement, iIsLoading, value->Position);
      DoIO(lNormalElement, iIsLoading, value->Normal);
    }

    void SimConfig::DoIO(TiXmlElement *iParent,
                         bool iIsLoading,
                         lb::boundaries::iolets::InOutLetFile* const value)
    {
      TiXmlElement* lPositionElement = GetChild(iParent, "position", iIsLoading);
      TiXmlElement* lNormalElement = GetChild(iParent, "normal", iIsLoading);
      TiXmlElement* lPressureElement = GetChild(iParent, "pressure", iIsLoading);

      DoIO(lPressureElement, "path", iIsLoading, value->PressureFilePath);

      DoIO(lPositionElement, iIsLoading, value->Position);
      DoIO(lNormalElement, iIsLoading, value->Normal);
    }

    void SimConfig::DoIO(TiXmlElement *iParent, bool iIsLoading, util::Vector3D<float> &iValue)
    {
      DoIO(iParent, "x", iIsLoading, iValue.x);
      DoIO(iParent, "y", iIsLoading, iValue.y);
      DoIO(iParent, "z", iIsLoading, iValue.z);
    }

    TiXmlElement *SimConfig::GetChild(TiXmlElement *iParent,
                                      std::string iChildNodeName,
                                      bool iIsLoading)
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
