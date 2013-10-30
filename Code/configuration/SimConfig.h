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
#include "util/Vector3D.h"
#include "lb/LbmParameters.h"
#include "lb/iolets/InOutLets.h"
#include "extraction/PropertyOutputFile.h"
#include "extraction/GeometrySelectors.h"
#include "io/xml/XmlAbstractionLayer.h"

namespace hemelb
{
  namespace configuration
  {
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

    class SimConfig
    {
      public:
        SimConfig(const std::string& path);
        ~SimConfig();

        void Save(std::string path); // TODO this method should be able to be CONST
        // but because it uses DoIo, which uses one function signature for both reading and writing, it cannot be.

        const util::Vector3D<float> & GetVisualisationCentre() const
        {
          return visualisationCentre;
        }
        const std::vector<lb::iolets::InOutLet*> & GetInlets() const
        {
          return inlets;
        }
        const std::vector<lb::iolets::InOutLet*> & GetOutlets() const
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
        PhysicalDistance GetVoxelSize() const
        {
          return voxelSizeMetres;
        }
        PhysicalPosition GetGeometryOrigin() const
        {
          return geometryOriginMetres;
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

        /**
         * Returns the pressure to be used to initialise all the fluid sites in the domain
         * @return initial pressure
         */
        LatticeDensity GetInitialPressure() const;

      protected:
        /**
         * Protected default ctor to allow derived test fixture classes to create mocks.
         */
        SimConfig();

      private:
        void DoIO(const io::xml::Element xmlNode);
        void DoIOForSimulation(const io::xml::Element simEl);
        void DoIOForGeometry(const io::xml::Element geometryEl);

        std::vector<lb::iolets::InOutLet*> DoIOForInOutlets(const io::xml::Element xmlNode);

        void DoIOForBaseInOutlet(const io::xml::Element& ioletEl, lb::iolets::InOutLet* value);

        lb::iolets::InOutLet* DoIOForPressureInOutlet(const io::xml::Element& ioletEl);
        lb::iolets::InOutLetCosine* DoIOForCosinePressureInOutlet(const io::xml::Element& ioletEl);
        lb::iolets::InOutLetFile* DoIOForFilePressureInOutlet(const io::xml::Element& ioletEl);
        lb::iolets::InOutLetMultiscale
            * DoIOForMultiscalePressureInOutlet(const io::xml::Element& ioletEl);

        lb::iolets::InOutLet* DoIOForVelocityInOutlet(const io::xml::Element& ioletEl);
        lb::iolets::InOutLetParabolicVelocity
            * DoIOForParabolicVelocityInOutlet(const io::xml::Element& ioletEl);
        /**
         * Reads/writes Womersley velocity inlet from/to XML file.
         *
         * @param parent parent XML element
         * @param isLoading whether the method is reading or writing
         * @param value womersley iolet instance to be configured
         */
        lb::iolets::InOutLetWomersleyVelocity
            * DoIOForWomersleyVelocityInOutlet(const io::xml::Element& ioletEl);

        void DoIOForProperties(const io::xml::Element& xmlNode);
        void DoIOForProperty(io::xml::Element xmlNode, bool isLoading);
        extraction::OutputField DoIOForPropertyField(const io::xml::Element& xmlNode);
        extraction::PropertyOutputFile
            * DoIOForPropertyOutputFile(const io::xml::Element& propertyoutputEl);
        extraction::StraightLineGeometrySelector
            * DoIOForLineGeometry(const io::xml::Element& xmlNode);
        extraction::PlaneGeometrySelector* DoIOForPlaneGeometry(const io::xml::Element&);
        extraction::SurfacePointSelector* DoIOForSurfacePoint(const io::xml::Element&);

        void DoIOForInitialConditions(io::xml::Element parent);
        void DoIOForVisualisation(const io::xml::Element& visEl);

        const std::string& xmlFilePath;
        io::xml::Document* rawXmlDoc;
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
        PhysicalPressure initialPressure; ///< Pressure used to initialise the domain

      protected:
        // These have to contain pointers because there are multiple derived types that might be
        // instantiated.
        std::vector<lb::iolets::InOutLet*> inlets;
        std::vector<lb::iolets::InOutLet*> outlets;
        PhysicalTime timeStepLength;
        unsigned long totalTimeSteps;
        unsigned long warmUpSteps;
        PhysicalDistance voxelSizeMetres;
        PhysicalPosition geometryOriginMetres;
    };
  }
}

#endif /* HEMELB_CONFIGURATION_SIMCONFIG_H */
