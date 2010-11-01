#include "SimConfig.h"
#include "xml/xpath_processor.h"
#include "xml/tinyxpath_conf.h"
#include "xml/ticpp.h"
#include "xml/xpath_static.h"

#include <string>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include "fileutils.h"
#include <cstdlib>
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
    for (int ii = 0; ii < Inlets.size(); ii++)
    {
      delete Inlets[ii];
    }
    for (int jj = 0; jj < Outlets.size(); jj++)
    {
      delete Outlets[jj];
    }
  }

  SimConfig *SimConfig::Load(char *iPath)
  {
    util::check_file(iPath);
    TiXmlDocument *lConfigFile = new TiXmlDocument();
    lConfigFile->LoadFile(iPath);

    SimConfig *lRet = new SimConfig();

    lRet->DoIO(lConfigFile->LastChild(), true);

    delete lConfigFile;

    return lRet;
  }

  void SimConfig::Save(char *iPath)
  {
    TiXmlDocument *lConfigFile = new TiXmlDocument(iPath);

    DoIO(lConfigFile, false);

    delete lConfigFile;
  }

  void SimConfig::DoIO(TiXmlNode *iTopNode, bool iIsLoading)
  {
    DoIO(TinyXPath::XAp_xpath_attribute(iTopNode, "//geometry/@voxelsize"),
         iIsLoading, VoxelSize);
    DoIO(TinyXPath::XAp_xpath_attribute(iTopNode, "//geometry/datafile/@path"),
         iIsLoading, DataFilePath);

    DoIO(TinyXPath::XNp_xpath_node(iTopNode, "//inlets"), iIsLoading, Inlets);
    DoIO(TinyXPath::XNp_xpath_node(iTopNode, "//outlets"), iIsLoading, Outlets);
    DoIO(TinyXPath::XNp_xpath_node(iTopNode, "//visualisation/centre"),
         iIsLoading, VisCentre);
    DoIO(
         TinyXPath::XAp_xpath_attribute(iTopNode,
                                        "//visualisation/orientation/@longitude"),
         iIsLoading, VisLongitude);
    DoIO(
         TinyXPath::XAp_xpath_attribute(iTopNode,
                                        "//visualisation/orientation/@latitude"),
         iIsLoading, VisLatitude);
    DoIO(TinyXPath::XAp_xpath_attribute(iTopNode,
                                        "//visualisation/display/@zoom"),
         iIsLoading, VisZoom);
    DoIO(TinyXPath::XAp_xpath_attribute(iTopNode,
                                        "//visualisation/display/@brightness"),
         iIsLoading, VisBrightness);
    DoIO(TinyXPath::XAp_xpath_attribute(iTopNode,
                                        "//visualisation/range/@maxvelocity"),
         iIsLoading, MaxVelocity);
    DoIO(TinyXPath::XAp_xpath_attribute(iTopNode,
                                        "//visualisation/range/@maxstress"),
         iIsLoading, MaxStress);
  }

  void SimConfig::DoIO(TiXmlAttribute *iXmlDoc, bool iIsLoading, float &value)
  {
    if (iIsLoading)
    {
      char *dummy;
      value = strtof(iXmlDoc->Value(), &dummy);
    }
    else
    {
      std::stringstream lTemp(std::stringstream::in);
      lTemp << value;
      iXmlDoc->SetValue(lTemp.str());
    }
  }

  void SimConfig::DoIO(TiXmlAttribute *iXmlDoc,
                       bool iIsLoading,
                       std::string &value)
  {
    if (iIsLoading)
    {
      value = std::string(iXmlDoc->Value());
    }
    else
    {
      iXmlDoc->SetValue(value);
    }
  }

  void SimConfig::DoIO(TiXmlNode *iXmlDoc, bool iIsLoading, std::vector<
      InOutLet*> &bResult)
  {
    if (iIsLoading)
    {
      const TiXmlNode *lCurrentLet = iXmlDoc->FirstChild();

      while (lCurrentLet != NULL)
      {
        InOutLet *lNew = new InOutLet();
        bResult.push_back(lNew);
        DoIO((TiXmlNode *) lCurrentLet, iIsLoading, lNew);
        lCurrentLet = lCurrentLet->NextSibling();
      }
    }
    else
    {
      for (int ii = 0; ii < bResult.size(); ii++)
      {
        // NB we're good up to 99 inlets here.
        char lTempString[11];
        sprintf(lTempString, "/inlet[%i]", ii + 1);
        TiXmlNode *lNew = TinyXPath::XNp_xpath_node(iXmlDoc, lTempString);
        DoIO(lNew, iIsLoading, bResult[ii]);
      }
    }
  }

  void SimConfig::DoIO(TiXmlNode *iXmlDoc, bool iIsLoading, InOutLet *value)
  {
    TiXmlNode * lPositionElement = TinyXPath::XNp_xpath_node(iXmlDoc,
                                                             "//position");
    TiXmlNode *lNormalElement = TinyXPath::XNp_xpath_node(iXmlDoc, "//normal");

    DoIO(TinyXPath::XAp_xpath_attribute(iXmlDoc, "//pressure/@mean"),
         iIsLoading, value->PMean);
    DoIO(TinyXPath::XAp_xpath_attribute(iXmlDoc, "//pressure/@amplitude"),
         iIsLoading, value->PAmp);
    DoIO(TinyXPath::XAp_xpath_attribute(iXmlDoc, "//pressure/@phase"),
         iIsLoading, value->PPhase);

    DoIO(lPositionElement, iIsLoading, value->Position);
    DoIO(lNormalElement, iIsLoading, value->Normal);
  }

  void SimConfig::DoIO(TiXmlNode *iXmlDoc, bool iIsLoading, Vector &iValue)
  {
    DoIO(TinyXPath::XAp_xpath_attribute(iXmlDoc, "//@x"), iIsLoading, iValue.x);
    DoIO(TinyXPath::XAp_xpath_attribute(iXmlDoc, "//@y"), iIsLoading, iValue.y);
    DoIO(TinyXPath::XAp_xpath_attribute(iXmlDoc, "//@z"), iIsLoading, iValue.z);
  }

}
