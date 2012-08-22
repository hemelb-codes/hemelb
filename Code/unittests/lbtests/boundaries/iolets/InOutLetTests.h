// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_LBTESTS_BOUNDARIES_IOLETS_INOUTLETTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_BOUNDARIES_IOLETS_INOUTLETTESTS_H
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "unittests/helpers/FolderTestFixture.h"
#include "lb/boundaries/iolets/InOutLets.h"
#include "resources/Resource.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      namespace boundaries
      {
        namespace iolets
        {
          using namespace lb::boundaries::iolets;
          using namespace resources;
          /***
           * Class asserting behaviour of in- and out- lets.
           */
          class InOutLetTests : public helpers::FolderTestFixture
          {
            CPPUNIT_TEST_SUITE(InOutLetTests);
                CPPUNIT_TEST(TestCosineConstruct);
                CPPUNIT_TEST(TestFileConstruct);
              CPPUNIT_TEST_SUITE_END();
            public:
              void setUp()
              {

              }
              void tearDown()
              {

              }
            private:
              void TestCosineConstruct()
              {

                // Bootstrap ourselves a in inoutlet, by loading config.xml.
                configuration::SimConfig *config =
                    configuration::SimConfig::Load(Resource("config.xml").Path().c_str());
                cosine = static_cast<InOutLetCosine*>(config->GetInlets()[0]);

                // Bootstrap ourselves a unit converter, which the cosine needs in initialisation
                lb::SimulationState state = lb::SimulationState(config->GetTimeStepLength(),
                                                                config->GetTotalTimeSteps());
                double voxelSize = 0.0001;
                lb::LbmParameters lbmParams = lb::LbmParameters(state.GetTimeStepLength(), voxelSize);
                util::UnitConverter converter = util::UnitConverter(&lbmParams, &state, voxelSize, PhysicalPosition());
                // at this stage, Initialise() has not been called, so the unit converter will be invalid, so we will not be able to convert to physical units.
                cosine->Initialise(&converter);
                cosine->Reset(state);

                // Check the cosine IOLET contains the values expected given the file.
                /*
                 * <inlet>
                 <pressure amplitude="0.0" mean="80.1" phase="0.0" period="0.6"/>
                 <normal x="0.0" y="0.0" z="1.0" />
                 <position x="-1.66017717834e-05" y="-4.58437586355e-05" z="-0.05" />
                 </inlet>
                 */
                CPPUNIT_ASSERT_EQUAL(80.1, cosine->GetPressureMean());
                CPPUNIT_ASSERT_EQUAL(0.0, cosine->GetPressureAmp());
                CPPUNIT_ASSERT_EQUAL(0.0, cosine->GetPhase());
                CPPUNIT_ASSERT_EQUAL(0.6, cosine->GetPeriod());
                CPPUNIT_ASSERT_EQUAL(util::Vector3D<float>(-1.66017717834e-05,-4.58437586355e-05,-0.05),
                                     cosine->GetPosition());
                CPPUNIT_ASSERT_EQUAL(util::Vector3D<float>(0.0,0.0,1.0), cosine->GetNormal());

                // Set an approriate target value for the density, the maximum.
                double temp = state.GetTimeStepLength() / voxelSize;
                double targetMeanDensity = 1
                    + (80.1 - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * temp * temp / (Cs2 * BLOOD_DENSITY_Kg_per_m3);

                // Check that the cosine formula correctly produces mean value
                CPPUNIT_ASSERT_EQUAL(targetMeanDensity, cosine->GetDensityMean());
                CPPUNIT_ASSERT_EQUAL(targetMeanDensity, cosine->GetDensity(0));
              }
              void TestFileConstruct()
              {

                // Bootstrap ourselves a file inlet, by loading an appropriate config file.
                // We have to move to a tempdir, as the path from the inlet to the iolet.txt file is a relative path

                FolderTestFixture::setUp();
                CopyResourceToTempdir("iolet.txt");
                MoveToTempdir();

                configuration::SimConfig *config =
                    configuration::SimConfig::Load(Resource("config_file_inlet.xml").Path().c_str());
                lb::SimulationState state = lb::SimulationState(config->GetTimeStepLength(),
                                                                config->GetTotalTimeSteps());
                double voxelSize = 0.0001;
                lb::LbmParameters lbmParams = lb::LbmParameters(state.GetTimeStepLength(), voxelSize);
                util::UnitConverter converter = util::UnitConverter(&lbmParams, &state, voxelSize, PhysicalPosition());
                file = static_cast<InOutLetFile*>(config->GetInlets()[0]);
                // at this stage, Initialise() has not been called, so the unit converter will be invalid, so we will not be able to convert to physical units.
                file->Initialise(&converter);
                file->Reset(state);


                // Ok, now we have an inlet, check the values are right.
                CPPUNIT_ASSERT_EQUAL(std::string("./iolet.txt"),file->GetFilePath());
                CPPUNIT_ASSERT_EQUAL(78.0, file->GetPressureMin());
                CPPUNIT_ASSERT_EQUAL(82.0, file->GetPressureMax());
                CPPUNIT_ASSERT_EQUAL(util::Vector3D<float>(-1.66017717834e-05,-4.58437586355e-05,-0.05),
                                     file->GetPosition());
                CPPUNIT_ASSERT_EQUAL(util::Vector3D<float>(0.0,0.0,1.0), file->GetNormal());

                // Set some target values for the density at various times.
                double temp = state.GetTimeStepLength() / voxelSize;
                double targetStartDensity = 1
                    + (78.0 - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * temp * temp / (Cs2 * BLOOD_DENSITY_Kg_per_m3);
                double targetMidDensity = 1
                    + (82.0 - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * temp * temp / (Cs2 * BLOOD_DENSITY_Kg_per_m3);

                CPPUNIT_ASSERT_DOUBLES_EQUAL(targetStartDensity, file->GetDensityMin(), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(targetStartDensity, file->GetDensity(0), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(targetMidDensity, file->GetDensity(state.GetTotalTimeSteps()/2), 1e-6);
                FolderTestFixture::tearDown();
              }
              InOutLetCosine *cosine;
              InOutLetFile *file;
          };
          CPPUNIT_TEST_SUITE_REGISTRATION(InOutLetTests);
        }
      }
    }
  }
}

#endif // HEMELB_UNITTESTS_LBTESTS_BOUNDARIES_IOLETS_INOUTLETTESTS_H
