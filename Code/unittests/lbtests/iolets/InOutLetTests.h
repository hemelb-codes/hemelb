
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_LBTESTS_IOLETS_INOUTLETTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_IOLETS_INOUTLETTESTS_H
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "unittests/helpers/FolderTestFixture.h"
#include "lb/iolets/InOutLets.h"
#include "resources/Resource.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      namespace iolets
      {
        using namespace lb::iolets;
        using namespace resources;

        class UncheckedSimConfig : public configuration::SimConfig
        {
          public:
            UncheckedSimConfig(const std::string& path) :
                configuration::SimConfig(path)
            {
              Init();
            }
          protected:
            virtual void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                                const std::string& requiredBC)
            {

            }
        };

        /***
         * Class asserting behaviour of in- and out- lets.
         */
        class InOutLetTests : public helpers::FolderTestFixture
        {
            CPPUNIT_TEST_SUITE(InOutLetTests);
            CPPUNIT_TEST(TestCosineConstruct);
            CPPUNIT_TEST(TestFileConstruct);
            CPPUNIT_TEST(TestIoletCoordinates);
            CPPUNIT_TEST(TestParabolicVelocityConstruct);
            CPPUNIT_TEST(TestWomersleyVelocityConstruct);
            CPPUNIT_TEST(TestFileVelocityConstruct);
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
              UncheckedSimConfig config(Resource("config.xml").Path());
              cosine = static_cast<InOutLetCosine*>(config.GetInlets()[0]);

              // Bootstrap ourselves a unit converter, which the cosine needs in initialisation
              lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                              config.GetTotalTimeSteps());
              double voxelSize = config.GetVoxelSize();
              const util::UnitConverter& converter = config.GetUnitConverter();
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
              CPPUNIT_ASSERT_DOUBLES_EQUAL(80.1,
                                           converter.ConvertPressureToPhysicalUnits(cosine->GetPressureMean()),
                                           1e-6);
              CPPUNIT_ASSERT_EQUAL(0.0, cosine->GetPressureAmp());
              CPPUNIT_ASSERT_EQUAL(0.0, cosine->GetPhase());
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.6,
                                           converter.ConvertTimeToPhysicalUnits(cosine->GetPeriod()),
                                           1e-6);

              PhysicalPosition expected(-1.66017717834e-05, -4.58437586355e-05, -0.05);
              PhysicalPosition actual =
                  converter.ConvertPositionToPhysicalUnits(cosine->GetPosition());
              for (unsigned i = 0; i < 3; ++i)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i], actual[i], 1e-9);

              CPPUNIT_ASSERT_EQUAL(util::Vector3D<Dimensionless>(0.0, 0.0, 1.0),
                                   cosine->GetNormal());

              // Set an approriate target value for the density, the maximum.
              double temp = state.GetTimeStepLength() / voxelSize;
              LatticeDensity targetMeanDensity = 1
                  + (80.1 - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * temp * temp
                      / (Cs2 * BLOOD_DENSITY_Kg_per_m3);
              // Check that the cosine formula correctly produces mean value
              CPPUNIT_ASSERT_DOUBLES_EQUAL(targetMeanDensity, cosine->GetDensityMean(), 1e-9);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(targetMeanDensity, cosine->GetDensity(0), 1e-9);
            }

            void TestFileConstruct()
            {

              // Bootstrap ourselves a file inlet, by loading an appropriate config file.
              // We have to move to a tempdir, as the path from the inlet to the iolet.txt file is a relative path

              FolderTestFixture::setUp();
              CopyResourceToTempdir("iolet.txt");
              MoveToTempdir();

              UncheckedSimConfig config(Resource("config_file_inlet.xml").Path());
              lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                              config.GetTotalTimeSteps());
              double voxelSize = config.GetVoxelSize();
              const util::UnitConverter& converter = config.GetUnitConverter();
              file = static_cast<InOutLetFile*>(config.GetInlets()[0]);
              // at this stage, Initialise() has not been called, so the unit converter will be invalid, so we will not be able to convert to physical units.
              file->Initialise(&converter);
              file->Reset(state);

              // Ok, now we have an inlet, check the values are right.
              CPPUNIT_ASSERT_EQUAL(std::string("./iolet.txt"), file->GetFilePath());
              CPPUNIT_ASSERT_DOUBLES_EQUAL(78.0,
                                           converter.ConvertPressureToPhysicalUnits(file->GetPressureMin()),
                                           1e-6);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(82.0,
                                           converter.ConvertPressureToPhysicalUnits(file->GetPressureMax()),
                                           1e-6);
              PhysicalPosition expected(-1.66017717834e-05, -4.58437586355e-05, -0.05);
              PhysicalPosition actual =
                  converter.ConvertPositionToPhysicalUnits(file->GetPosition());
              CPPUNIT_ASSERT_EQUAL(expected, actual);
              CPPUNIT_ASSERT_EQUAL(util::Vector3D<Dimensionless>(0.0, 0.0, 1.0), file->GetNormal());

              // Set some target values for the density at various times.
              double temp = state.GetTimeStepLength() / voxelSize;
              double targetStartDensity = 1
                  + (78.0 - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * temp * temp
                      / (Cs2 * BLOOD_DENSITY_Kg_per_m3);
              double targetMidDensity = 1
                  + (82.0 - REFERENCE_PRESSURE_mmHg) * mmHg_TO_PASCAL * temp * temp
                      / (Cs2 * BLOOD_DENSITY_Kg_per_m3);

              CPPUNIT_ASSERT_DOUBLES_EQUAL(targetStartDensity, file->GetDensityMin(), 1e-6);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(targetStartDensity, file->GetDensity(0), 1e-6);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(targetMidDensity,
                                           file->GetDensity(state.GetTotalTimeSteps() / 2),
                                           1e-6);
              FolderTestFixture::tearDown();
            }

            void TestParabolicVelocityConstruct()
            {

              // Bootstrap ourselves a in inoutlet, by loading config.xml.
              UncheckedSimConfig config(Resource("config-velocity-iolet.xml").Path());
              p_vel = dynamic_cast<InOutLetParabolicVelocity*>(config.GetInlets()[0]);
              CPPUNIT_ASSERT(p_vel != NULL);

              // Bootstrap ourselves a unit converter, which the cosine needs in initialisation
              lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                              config.GetTotalTimeSteps());
              const util::UnitConverter& converter = config.GetUnitConverter();
              // at this stage, Initialise() has not been called, so the unit converter will be invalid, so we will not be able to convert to physical units.
              p_vel->Initialise(&converter);
              p_vel->Reset(state);

              // Check the IOLET contains the values expected given the file.
              /*
               * <inlet>
               <velocity radius="0.005" maximum="0.10"/>
               <normal x="0.0" y="0.0" z="1.0" />
               <position x="-1.66017717834e-05" y="-4.58437586355e-05" z="-0.05" />
               </inlet>
               */
              CPPUNIT_ASSERT_EQUAL(5.0, p_vel->GetRadius());
              CPPUNIT_ASSERT_EQUAL(0.01, p_vel->GetMaxSpeed());
              CPPUNIT_ASSERT_EQUAL(PhysicalPosition(-1.66017717834e-05, -4.58437586355e-05, -0.05),
                                   converter.ConvertPositionToPhysicalUnits(p_vel->GetPosition()));
              CPPUNIT_ASSERT_EQUAL(util::Vector3D<Dimensionless>(0.0, 0.0, 1.0),
                                   p_vel->GetNormal());
            }

            void TestWomersleyVelocityConstruct()
            {

              // Bootstrap ourselves a in inoutlet, by loading config.xml.
              UncheckedSimConfig config(Resource("config_new_velocity_inlets.xml").Path());
              womersVel = static_cast<InOutLetWomersleyVelocity*>(config.GetInlets()[0]);

              // Bootstrap ourselves a unit converter, which the cosine needs in initialisation
              lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                              config.GetTotalTimeSteps());
              double voxelSize = config.GetVoxelSize();
              const util::UnitConverter& converter = config.GetUnitConverter();
              // at this stage, Initialise() has not been called, so the unit converter will be invalid, so we will not be able to convert to physical units.
              womersVel->Initialise(&converter);
              womersVel->Reset(state);

              // Check the IOLET contains the values expected given the file.
              CPPUNIT_ASSERT_EQUAL(10.0, womersVel->GetRadius());
              CPPUNIT_ASSERT_EQUAL(mmHg_TO_PASCAL * 1e-6,
                                   womersVel->GetPressureGradientAmplitude());
              CPPUNIT_ASSERT_EQUAL(5.0, womersVel->GetPeriod());
              CPPUNIT_ASSERT_EQUAL(2.0, womersVel->GetWomersleyNumber());
              CPPUNIT_ASSERT_EQUAL(PhysicalPosition(0, 0, -0.05),
                                   converter.ConvertPositionToPhysicalUnits(womersVel->GetPosition()));
              CPPUNIT_ASSERT_EQUAL(util::Vector3D<Dimensionless>(0.0, 0.0, 1.0),
                                   womersVel->GetNormal());

              /*
               *  Test that the analytical solution at r=R is 0
               */
              PhysicalPosition pointAtCylinderWall(womersVel->GetRadius() * voxelSize, 0, -0.05);
              LatticePosition pointAtCylinderWallLatticeUnits(converter.ConvertPositionToLatticeUnits(pointAtCylinderWall));
              LatticeVelocity zeroVelAtWall(womersVel->GetVelocity(pointAtCylinderWallLatticeUnits,
                                                                   0));
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, zeroVelAtWall[0], 1e-9);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, zeroVelAtWall[1], 1e-9);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, zeroVelAtWall[2], 1e-9);

              /*
               * With a small enough Womersley number, the solution should match the Poiseuille solution for the same pressure
               * difference after pi/2 radians and have changed direction after 3*pi/2 radians.
               */
              Dimensionless eta(1), nu(1);

              Dimensionless alpha = 1e-4;
              womersVel->SetWomersleyNumber(alpha);
              womersVel->SetPeriod(2 * PI * pow(womersVel->GetRadius(), 2) / (alpha * alpha * nu));

              LatticeSpeed poiseuilleSolution = womersVel->GetPressureGradientAmplitude()
                  * pow(womersVel->GetRadius(), 2) / (4 * eta);

              LatticePosition pointAtCentrelineLatticeUnits(womersVel->GetPosition());

              {
                LatticeVelocity poiseuilleVelAtCentreLine(womersVel->GetVelocity(pointAtCentrelineLatticeUnits,
                                                                                 0.25
                                                                                     * womersVel->GetPeriod()));
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, poiseuilleVelAtCentreLine[0], 1e-9);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, poiseuilleVelAtCentreLine[1], 1e-9);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(poiseuilleSolution,
                                             poiseuilleVelAtCentreLine[2],
                                             1e-9);
              }

              {
                LatticeVelocity poiseuilleVelAtCentreLine(womersVel->GetVelocity(pointAtCentrelineLatticeUnits,
                                                                                 0.75
                                                                                     * womersVel->GetPeriod()));
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, poiseuilleVelAtCentreLine[0], 1e-9);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, poiseuilleVelAtCentreLine[1], 1e-9);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(-poiseuilleSolution,
                                             poiseuilleVelAtCentreLine[2],
                                             1e-9);
              }

            }

            void TestFileVelocityConstruct()
            {
              // We have to move to a tempdir, as the path specified in the xml file is a relative path
              FolderTestFixture::setUp();
              CopyResourceToTempdir("velocity_inlet.txt");
              MoveToTempdir();

              // Bootstrap ourselves a file velocity inlet, by loading an appropriate config file.
              UncheckedSimConfig config(Resource("config_file_velocity_inlet.xml").Path());
              lb::SimulationState state = lb::SimulationState(config.GetTimeStepLength(),
                                                              config.GetTotalTimeSteps());
              double voxelSize = config.GetVoxelSize();
              const util::UnitConverter& converter = config.GetUnitConverter();
              fileVel = static_cast<InOutLetFileVelocity*>(config.GetInlets()[0]);
              // at this stage, Initialise() has not been called, so the unit converter will be invalid, so we will not be able to convert to physical units.
              fileVel->Initialise(&converter);
              fileVel->Reset(state);

              int FilePath_length = fileVel->GetFilePath().length();
              CPPUNIT_ASSERT_EQUAL(fileVel->GetFilePath().substr(FilePath_length-19), std::string("/velocity_inlet.txt"));
              CPPUNIT_ASSERT_EQUAL(fileVel->GetRadius(), 20.0);
              PhysicalPosition expected(0, 0, -0.05);
              PhysicalPosition actual =
                  converter.ConvertPositionToPhysicalUnits(fileVel->GetPosition());
              CPPUNIT_ASSERT_EQUAL(expected, actual);
              CPPUNIT_ASSERT_EQUAL(util::Vector3D<Dimensionless>(0.0, 0.0, 1.0),
                                   fileVel->GetNormal());

              LatticePosition pointAtCentrelineLatticeUnits(fileVel->GetPosition());

              {
                LatticeVelocity velAtCentreLine(fileVel->GetVelocity(pointAtCentrelineLatticeUnits,
                                                                     converter.ConvertTimeToLatticeUnits(3.0)));
                PhysicalVelocity physVelCentreLine =
                    converter.ConvertVelocityToPhysicalUnits(velAtCentreLine);

                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, physVelCentreLine[0], 1e-9);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, physVelCentreLine[1], 1e-9);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.01, physVelCentreLine[2], 1e-9);
              }

              // Point equidistant to the centreline and wall
              LatticePosition pointEquidistant = fileVel->GetPosition();
              pointEquidistant[0] -= fileVel->GetRadius() / 2.0;

              {
                LatticeVelocity velAtPointEquidistant(fileVel->GetVelocity(pointEquidistant,
                                                                           converter.ConvertTimeToLatticeUnits(3.0)));
                PhysicalVelocity physVelPointEqui =
                    converter.ConvertVelocityToPhysicalUnits(velAtPointEquidistant);

                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, physVelPointEqui[0], 1e-9);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, physVelPointEqui[1], 1e-9);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0075, physVelPointEqui[2], 1e-9);
              }

              FolderTestFixture::tearDown();
            }

            class ConcreteIolet : public InOutLet
            {
                virtual InOutLet* Clone() const
                {
                  ConcreteIolet* copy = new ConcreteIolet(*this);
                  return copy;
                }
                virtual LatticeDensity GetDensityMin() const
                {
                  return 1.0;
                }
                virtual LatticeDensity GetDensityMax() const
                {
                  return 1.0;
                }
                virtual LatticeDensity GetDensity(hemelb::LatticeTimeStep) const
                {
                  return 1.0;
                }
                virtual void Reset(hemelb::lb::SimulationState&)
                {
                }
            };

            void TestIoletCoordinates()
            {
              // unit converter - make physical and lattice units the same
              hemelb::util::UnitConverter units(1, 1, PhysicalPosition::Zero());

              ConcreteIolet iolet;
              // normal
              util::Vector3D<Dimensionless> n(5, 7, -4);
              n.Normalise();
              iolet.SetNormal(n);
              // position
              PhysicalPosition c(7.77438796, 9.21293516, 9.87122463);
              iolet.SetPosition(units.ConvertPositionToLatticeUnits(c));
              iolet.Initialise(&units);
              IoletExtraData extra(iolet);
              iolet.SetExtraData(&extra);

              // Convert the centre to iolet coords
              LatticePosition tmp = extra.WorldToIolet(c);
              // This should be zero
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmp.x, 1e-9);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmp.y, 1e-9);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmp.z, 1e-9);

              // Make a point 3 lattice units along the normal.
              LatticePosition zEqThree = c + n * 3.0;
              tmp = extra.WorldToIolet(zEqThree);
              // This should be zero in x & y but 3 in z
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmp.x, 1e-9);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, tmp.y, 1e-9);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, tmp.z, 1e-9);

            }

            InOutLetCosine *cosine;
            InOutLetFile *file;
            InOutLetParabolicVelocity* p_vel;
            InOutLetWomersleyVelocity* womersVel;
            InOutLetFileVelocity* fileVel;
        };
        CPPUNIT_TEST_SUITE_REGISTRATION(InOutLetTests);
      }
    }
  }
}

#endif // HEMELB_UNITTESTS_LBTESTS_IOLETS_INOUTLETTESTS_H
