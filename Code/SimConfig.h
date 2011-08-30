#ifndef HEMELB_SIMCONFIG_H
#define HEMELB_SIMCONFIG_H

#include <string.h>
#include <vector>
#include "xml/tinyxml.h"

namespace hemelb
{
  class SimConfig
  {
    public:

      struct Vector
      {
          float x;
          float y;
          float z;
      };

      struct InOutLet
      {
          double PMean;
          double PAmp;
          double PPhase;
          double PMin;
          double PMax;
          Vector Position;
          Vector Normal;
          std::string PFilePath;
      };

      static SimConfig *Load(const char *iPath);
      ~SimConfig();

      void Save(std::string iPath);

      std::string DataFilePath;
      std::vector<InOutLet> Inlets;
      std::vector<InOutLet> Outlets;
      Vector VisCentre;
      float VisLongitude;
      float VisLatitude;
      float VisZoom;
      float VisBrightness;
      float MaxVelocity;
      float MaxStress;
      unsigned long NumCycles;
      long StepsPerCycle;

    private:
      SimConfig();
      void DoIO(TiXmlElement *iXmlNode, bool iIsLoading);
      void DoIO(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, long &value);
      void DoIO(TiXmlElement* iXmlNode,
                std::string iAttributeName,
                bool iIsLoading,
                unsigned long &value);
      void DoIO(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, float &value);
      void DoIO(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, double &value);
      void DoIO(TiXmlElement *iXmlNode,
                std::string iAttributeName,
                bool iIsLoading,
                std::string &iValue);
      void DoIO(TiXmlElement *iXmlNode,
                bool iIsLoading,
                std::vector<InOutLet> &value,
                std::string iChildNodeName);
      void DoIO(TiXmlElement *iXmlNode, bool iIsLoading, InOutLet &value);
      void DoIO(TiXmlElement *iXmlNode, bool iIsLoading, Vector &iValue);
      TiXmlElement* GetChild(TiXmlElement *iParent, std::string iChildNodeName, bool iIsLoading);
  };
}

#endif /* HEMELB_SIMCONFIG_H */
