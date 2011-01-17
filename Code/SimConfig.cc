#include "SimConfig.h"

#include <string>
#include <iostream>
#include <unistd.h>
#include <cstdlib>

#include "util/fileutils.h"
#include "debug/Debugger.h"

namespace hemelb
{

  SimConfig::SimConfig()
  {
    // This constructor only exists to prevent instantiation without
    // using the static load method.
  }

  SimConfig::~SimConfig()
  {
    for (unsigned int ii = 0; ii < Inlets.size(); ii++)
    {
      delete Inlets[ii];
    }
    for (unsigned int jj = 0; jj < Outlets.size(); jj++)
    {
      delete Outlets[jj];
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

  void SimConfig::DoIO(TiXmlElement *iTopNode, bool iIsLoading)
  {
    TiXmlElement* lSimulationElement = GetChild(iTopNode, "simulation",
                                                iIsLoading);
    DoIO(lSimulationElement, "cycles", iIsLoading, NumCycles);
    DoIO(lSimulationElement, "cyclesteps", iIsLoading, StepsPerCycle);

    TiXmlElement* lGeometryElement = GetChild(iTopNode, "geometry", iIsLoading);
    DoIO(lGeometryElement, "voxelsize", iIsLoading, VoxelSize);
    DoIO(GetChild(lGeometryElement, "datafile", iIsLoading), "path",
         iIsLoading, DataFilePath);

    DoIO(GetChild(iTopNode, "inlets", iIsLoading), iIsLoading, Inlets, "inlet");

    DoIO(GetChild(iTopNode, "outlets", iIsLoading), iIsLoading, Outlets,
         "outlet");

    TiXmlElement* lVisualisationElement = GetChild(iTopNode, "visualisation",
                                                   iIsLoading);
    DoIO(GetChild(lVisualisationElement, "centre", iIsLoading), iIsLoading,
         VisCentre);
    TiXmlElement *lOrientationElement = GetChild(lVisualisationElement,
                                                 "orientation", iIsLoading);
    DoIO(lOrientationElement, "longitude", iIsLoading, VisLongitude);
    DoIO(lOrientationElement, "latitude", iIsLoading, VisLatitude);

    TiXmlElement *lDisplayElement = GetChild(lVisualisationElement, "display",
                                             iIsLoading);

    DoIO(lDisplayElement, "zoom", iIsLoading, VisZoom);
    DoIO(lDisplayElement, "brightness", iIsLoading, VisBrightness);

    TiXmlElement *lRangeElement = GetChild(lVisualisationElement, "range",
                                           iIsLoading);

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
      value = std::strtod(iParent->Attribute(iAttributeName)->c_str(), &dummy);
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
      char *dummy;
      value = std::strtod(iParent->Attribute(iAttributeName)->c_str(), &dummy);
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
      iValue = std::string(iParent->Attribute(iAttributeName)->c_str());
    }
    else
    {
      iParent->SetAttribute(iAttributeName, iValue);
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

  void SimConfig::DoIO(TiXmlElement *iParent, bool iIsLoading, std::vector<
      InOutLet*> &bResult, std::string iChildNodeName)
  {
    if (iIsLoading)
    {
      TiXmlElement *lCurrentLet = iParent->FirstChildElement(iChildNodeName);

      while (lCurrentLet != NULL)
      {
        InOutLet *lNew = new InOutLet();
        bResult.push_back(lNew);
        DoIO(lCurrentLet, iIsLoading, lNew);
        lCurrentLet = lCurrentLet->NextSiblingElement(iChildNodeName);
      }
    }
    else
    {
      for (unsigned int ii = 0; ii < bResult.size(); ii++)
      {
        // NB we're good up to 99 inlets here.
        DoIO(GetChild(iParent, iChildNodeName, iIsLoading), iIsLoading,
             bResult[ii]);
      }
    }
  }

  void SimConfig::DoIO(TiXmlElement *iParent, bool iIsLoading, InOutLet *value)
  {
    TiXmlElement* lPositionElement = GetChild(iParent, "position", iIsLoading);
    TiXmlElement* lNormalElement = GetChild(iParent, "normal", iIsLoading);
    TiXmlElement* lPressureElement = GetChild(iParent, "pressure", iIsLoading);

    DoIO(lPressureElement, "mean", iIsLoading, value->PMean);
    DoIO(lPressureElement, "amplitude", iIsLoading, value->PAmp);
    DoIO(lPressureElement, "phase", iIsLoading, value->PPhase);

    DoIO(lPositionElement, iIsLoading, value->Position);
    DoIO(lNormalElement, iIsLoading, value->Normal);
  }

  void SimConfig::DoIO(TiXmlElement *iParent, bool iIsLoading, Vector &iValue)
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
