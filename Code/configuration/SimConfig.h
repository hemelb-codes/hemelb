// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_CONFIGURATION_SIMCONFIG_H
#define HEMELB_CONFIGURATION_SIMCONFIG_H

#include <vector>
#include "tinyxml.h"
#include "util/Vector3D.h"
#include "lb/boundaries/iolets/InOutLets.h"
#include "extraction/PropertyOutputFile.h"
#include "extraction/GeometrySelectors.h"

namespace hemelb
{
  namespace configuration
  {
    class SimConfig
    {
      public:
        static SimConfig *Load(const char *path);
        ~SimConfig();

        void Save(std::string path); // TODO this method should be able to be CONST
        // but because it uses DoIo, which uses one function signature for both reading and writing, it cannot be.

        void
        DoIOForCosineInOutlet(TiXmlElement *xmlNode, bool isLoading, lb::boundaries::iolets::InOutLetCosine* value);
        void DoIOForFileInOutlet(TiXmlElement *xmlNode, bool isLoading, lb::boundaries::iolets::InOutLetFile* value);
        void DoIOForMultiscaleInOutlet(TiXmlElement *xmlNode,
                                       bool isLoading,
                                       lb::boundaries::iolets::InOutLetMultiscale* value);
        const util::Vector3D<float> & GetVisualisationCentre() const
        {
          return visualisationCentre;
        }
        const std::vector<lb::boundaries::iolets::InOutLet*> & GetInlets() const
        {
          return inlets;
        }
        const std::vector<lb::boundaries::iolets::InOutLet*> & GetOutlets() const
        {
          return outlets;
        }
        lb::StressTypes GetStressType() const
        {
          return stressType;
        }
        float GetVisualisationLongitude() const
        {
          return visualisationLongitude;
        }
        float GetVisualisationLatitude() const
        {
          return visualisationLatitude;
        }
        float GetVisualisationZoom() const
        {
          return visualisationZoom;
        }
        float GetVisualisationBrightness() const
        {
          return visualisationBrightness;
        }
        float GetMaximumVelocity() const
        {
          return maxVelocity;
        }
        float GetMaximumStress() const
        {
          return maxStress;
        }
        const std::string & GetDataFilePath() const
        {
          return dataFilePath;
        }
        LatticeTime GetTotalTimeSteps() const
        {
          return totalTimeSteps;
        }
        LatticeTime GetWarmUpSteps() const
        {
          return warmUpSteps;
        }
        PhysicalTime GetTimeStepLength() const
        {
          return timeStepLength;
        }
        unsigned int PropertyOutputCount() const
        {
          return propertyOutputs.size();
        }
        extraction::PropertyOutputFile * GetPropertyOutput(unsigned int index) const
        {
          return propertyOutputs[index];
        }
        std::vector<extraction::PropertyOutputFile*> const GetPropertyOutputs() const
        {
          return propertyOutputs;
        }
        const std::string GetColloidConfigPath() const
        {
          return colloidConfigPath;
        }
        /**
         * True if the XML file has a section specifying colloids.
         * @return
         */
        bool HasColloidSection() const;

      protected:
        SimConfig();

      private:
        void DoIO(TiXmlElement *xmlNode, bool isLoading);
        bool DoIOForLong(TiXmlElement* xmlNode, std::string attributeName, bool isLoading, long &value);
        void DoIOForStressType(TiXmlElement* xmlNode,
                               std::string attributeName,
                               bool isLoading,
                               lb::StressTypes &value);
        bool DoIOForULong(TiXmlElement* xmlNode, std::string attributeName, bool isLoading, unsigned long &value);
        void DoIOForFloat(TiXmlElement* xmlNode, std::string attributeName, bool isLoading, float &value);
        bool DoIOForDouble(TiXmlElement* xmlNode, std::string attributeName, bool isLoading, double &value);
        void DoIOForString(TiXmlElement* xmlNode, std::string attributeName, bool isLoading, std::string &value);
        void DoIOForInOutlets(TiXmlElement *xmlNode,
                              bool isLoading,
                              std::vector<lb::boundaries::iolets::InOutLet*> &value,
                              std::string childNodeName);
        void DoIOForProperties(TiXmlElement *xmlNode, bool isLoading);
        void DoIOForProperty(TiXmlElement *xmlNode, bool isLoading);
        void DoIOForPropertyField(TiXmlElement *xmlNode, bool isLoading, extraction::OutputField& field);
        void DoIOForPropertyOutputFile(TiXmlElement *xmlNode, bool isLoading, extraction::PropertyOutputFile* file);
        void DoIOForLineGeometry(TiXmlElement *xmlNode,
                                 bool isLoading,
                                 extraction::StraightLineGeometrySelector*& line);
        void DoIOForPlaneGeometry(TiXmlElement *xmlNode, bool isLoading, extraction::PlaneGeometrySelector*& plane);
        void DoIOForSurfacePoint(TiXmlElement *xmlNode, bool isLoading, extraction::SurfacePointSelector*& plane);

        void DoIOForFloatVector(TiXmlElement *xmlNode, bool isLoading, util::Vector3D<float> &value);
        void DoIOForBaseInOutlet(TiXmlElement *parent, bool isLoading, lb::boundaries::iolets::InOutLet* const value);
        TiXmlElement* GetChild(TiXmlElement *parent, std::string childNodeName, bool isLoading);

        const double LEGACY_PULSATILE_PERIOD;
        std::string dataFilePath;

        util::Vector3D<float> visualisationCentre;
        float visualisationLongitude;
        float visualisationLatitude;
        float visualisationZoom;
        float visualisationBrightness;
        float maxVelocity;
        float maxStress;
        lb::StressTypes stressType;
        std::vector<extraction::PropertyOutputFile*> propertyOutputs;
        std::string colloidConfigPath;
        /**
         * True if the file has a colloids section.
         */
        bool hasColloidSection;

      protected:
        // These have to contain pointers because there are multiple derived types that might be
        // instantiated.
        std::vector<lb::boundaries::iolets::InOutLet*> inlets;
        std::vector<lb::boundaries::iolets::InOutLet*> outlets;
        double timeStepLength;
        unsigned long totalTimeSteps;
        unsigned long warmUpSteps;

    };
  }
}

#endif /* HEMELB_CONFIGURATION_SIMCONFIG_H */
