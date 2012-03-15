#ifndef HEMELB_CONFIGURATION_SIMCONFIG_H
#define HEMELB_CONFIGURATION_SIMCONFIG_H

#include <vector>
#include "tinyxml.h"
#include "util/Vector3D.h"
#include "lb/boundaries/iolets/InOutLets.h"
#include "extraction/PropertyOutputFile.h"
#include "extraction/StraightLineGeometrySelector.h"
#include "extraction/PlaneGeometrySelector.h"
#include "extraction/WholeGeometrySelector.h"

namespace hemelb
{
  namespace configuration
  {
    class SimConfig
    {
      public:
        static SimConfig *Load(const char *iPath);
        ~SimConfig();

        void Save(std::string iPath); // TODO this method should be able to be CONST
        // but because it uses DoIo, which uses one function signature for both reading and writing, it cannot be.

        std::string DataFilePath;
        // These have to contain pointers because there are multiple derived types that might be
        // instantiated.
        std::vector<lb::boundaries::iolets::InOutLet*> Inlets;
        std::vector<lb::boundaries::iolets::InOutLet*> Outlets;
        util::Vector3D<float> VisCentre;
        float VisLongitude;
        float VisLatitude;
        float VisZoom;
        float VisBrightness;
        float MaxVelocity;
        float MaxStress;
        double TimeStepLength;
        unsigned long TotalTimeSteps;
        lb::StressTypes StressType;
        std::vector<extraction::PropertyOutputFile*> propertyOutputs;

        void DoIOForCosineInOutlet(TiXmlElement *iXmlNode,
                                   bool iIsLoading,
                                   lb::boundaries::iolets::InOutLetCosine* value);
        void DoIOForFileInOutlet(TiXmlElement *iXmlNode, bool iIsLoading, lb::boundaries::iolets::InOutLetFile* value);
      protected:
        SimConfig();

      private:
        void DoIO(TiXmlElement *iXmlNode, bool iIsLoading);
        bool DoIOForLong(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, long &value);
        void DoIOForStressType(TiXmlElement* iXmlNode,
                               std::string iAttributeName,
                               bool iIsLoading,
                               lb::StressTypes &value);
        bool DoIOForULong(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, unsigned long &value);
        void DoIOForFloat(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, float &value);
        bool DoIOForDouble(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, double &value);
        void DoIOForString(TiXmlElement* iXmlNode, std::string iAttributeName, bool iIsLoading, std::string &iValue);
        void DoIOForInOutlets(TiXmlElement *iXmlNode,
                              bool iIsLoading,
                              std::vector<lb::boundaries::iolets::InOutLet*> &value,
                              std::string iChildNodeName);
        void DoIOForProperties(TiXmlElement *iXmlNode, bool iIsLoading);
        void DoIOForProperty(TiXmlElement *iXmlNode, bool iIsLoading);
        void DoIOForPropertyField(TiXmlElement *iXmlNode, bool iIsLoading, extraction::OutputField& field);
        void DoIOForPropertyOutputFile(TiXmlElement *iXmlNode, bool iIsLoading, extraction::PropertyOutputFile* file);
        void DoIOForLineGeometry(TiXmlElement *iXmlNode, bool iIsLoading, extraction::StraightLineGeometrySelector*& line);
        void DoIOForPlaneGeometry(TiXmlElement *iXmlNode, bool iIsLoading, extraction::PlaneGeometrySelector*& plane);

        void DoIOForFloatVector(TiXmlElement *iXmlNode, bool iIsLoading, util::Vector3D<float> &iValue);
        TiXmlElement* GetChild(TiXmlElement *iParent, std::string iChildNodeName, bool iIsLoading);
        const double LEGACY_PULSATILE_PERIOD;
    };
  }
}

#endif /* HEMELB_CONFIGURATION_SIMCONFIG_H */
