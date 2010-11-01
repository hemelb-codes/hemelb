#ifndef HEMELB_SIMCONFIG_H
#define HEMELB_SIMCONFIG_H

#include <string.h>
#include <vector>
#include "xml/ticpp.h"
#include "xml/xpath_static.h"

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
          float PMean;
          float PAmp;
          float PPhase;
          Vector Position;
          Vector Normal;
      };

      static SimConfig *Load(char *iPath);
      ~SimConfig();

      void Save(char *iPath);

      float VoxelSize;
      std::string DataFilePath;
      std::vector<InOutLet*> Inlets;
      std::vector<InOutLet*> Outlets;
      Vector VisCentre;
      float VisLongitude;
      float VisLatitude;
      float VisZoom;
      float VisBrightness;
      float MaxVelocity;
      float MaxStress;

    private:
      SimConfig();
      void DoIO(TiXmlNode *iXmlNode, bool iIsLoading);
      void DoIO(TiXmlAttribute *iXmlNode, bool iIsLoading, float &value);
      void DoIO(TiXmlAttribute *iXmlNode, bool iIsLoading, std::string &value);
      void DoIO(TiXmlNode *iXmlNode,
                bool iIsLoading,
                std::vector<InOutLet*> &value);
      void DoIO(TiXmlNode *iXmlNode, bool iIsLoading, InOutLet *value);
      void DoIO(TiXmlNode *iXmlNode, bool iIsLoading, Vector &iValue);
  };
}

#endif /* HEMELB_SIMCONFIG_H */
