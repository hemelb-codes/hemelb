#ifndef HEMELB_UNITTEST_TEST_INTERCOMMUNICATOR_EXAMPLE_H
#define HEMELB_UNITTEST_TEST_INTERCOMMUNICATOR_EXAMPLE_H
#include "unittests/multiscale/MockIntercommunicand.h"
#include "resources/Resource.h"
#include "unittests/multiscale/MockIntercommunicator.h"

using namespace hemelb::unittest::multiscale;
#include <iostream>

namespace hemelb
{
  namespace unittests
  {
    namespace multiscale
    {
      template<class IntercommuniatorImplementation> class MockHemeLB
      {
        public:
          MockHemeLB(double aspace_resolution, double atime_resolution, std::map<std::string, double> & buffer) :
              inlet(), outlet(), inoutlettype("inoutlet"), intercomms(buffer), time_resolution(atime_resolution), space_resolution(aspace_resolution), current_time(0)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inoutlettype.RegisterSharedValue<double>("pressure");
            inoutlettype.RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inoutlettype, inlet, "boundary1");
            intercomms.RegisterIntercommunicand(inoutlettype, outlet, "boundary2");
          }
          MockIntercommunicand inlet;
          MockIntercommunicand outlet;
          MockIntercommunicator::IntercommunicandTypeT inoutlettype;
          MockIntercommunicator intercomms;
          double time_resolution;
          double space_resolution;
          double current_time;

          void DoLB()
          {
            double resistance = 10.0;
            double velocity = (inlet.GetDensity() - outlet.GetDensity()) / resistance;
            inlet.SetVelocity(velocity);
            outlet.SetVelocity(velocity);
          }

          void Simulate()
          {
            if (!intercomms.ShouldAdvance())
              return;
            intercomms.GetFromMultiscale();
            DoLB();
            current_time += time_resolution;
            intercomms.AdvanceTime(current_time);
            intercomms.SendToMultiscale();
          }

      };
      template<class IntercommuniatorImplementation> class Mock0DModel
      {
        public:
          Mock0DModel(double aspace_resolution, double atime_resolution, std::map<std::string, double> & buffer) :
              inlet(), outlet(), inoutlettype("inoutlet"), intercomms(buffer), time_resolution(atime_resolution), space_resolution(aspace_resolution), current_time(0)
          {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inoutlettype.RegisterSharedValue<double>("pressure");
            inoutlettype.RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inoutlettype, inlet, "boundary2");
            intercomms.RegisterIntercommunicand(inoutlettype, outlet, "boundary1");
          }
          MockIntercommunicand inlet;
          MockIntercommunicand outlet;
          MockIntercommunicator::IntercommunicandTypeT inoutlettype;
          MockIntercommunicator intercomms;
          double time_resolution;
          double space_resolution;
          double current_time;

          void Do1D()
          {
            double capacitance = 10.0;
            double deltap = -1.0 * inlet.GetVelocity() / capacitance;

            outlet.SetDensity(outlet.GetDensity() + deltap);
          }

          void Simulate()
          {
            if (!intercomms.ShouldAdvance())
              return;
            intercomms.GetFromMultiscale();
            Do1D();
            current_time += time_resolution;
            intercomms.AdvanceTime(current_time);
            intercomms.SendToMultiscale();
          }

      };

      class MockIntercommunicatorTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(MockIntercommunicatorTests);
          CPPUNIT_TEST(TestCRRun);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            pbuffer = new std::map<std::string, double>();
            std::map<std::string, double> &buffer = *pbuffer;
            buffer["boundary1_pressure"] = 120.0;
            buffer["boundary2_pressure"] = 80.0;
            buffer["boundary1_velocity"] = 0.0;
            buffer["boundary2_velocity"] = 0.0;
            buffer["shared_time"] = 0.0;
            mockheme = new MockHemeLB<MockIntercommunicator>(25.0, 25.0, buffer);
            zerod = new Mock0DModel<MockIntercommunicator>(10.0, 10.0, buffer);

          }
          void tearDown()
          {

          }
        private:
          void TestCRRun()
          {
            // In this example, we imagine we are coupling inlets in two different instances of HemeLB.
            // Each has it's own instance of an intercommunicator - these would be running in different processes.
            // Our example intercommunicator has been hard coded with how to share the values.
            // In principle, though, this would be defined by whatever configuration file format the implementing intercommunicator wants to use.
            for (int steps = 0; steps < 100; steps++)
            {
              mockheme->Simulate();
              zerod->Simulate();
              //std::cout << "Heme: I: p:" << heme->inlet.GetDensity() << " v:" << heme->inlet.GetVelocity() << " O: p: "
              //    << heme->outlet.GetDensity() << " v:" << heme->outlet.GetVelocity() << std::endl;
              //std::cout << "Zerod: I: p:" << zerod->inlet.GetDensity() << " v:" << zerod->inlet.GetVelocity()
              //    << " O: p: " << zerod->outlet.GetDensity() << " v:" << zerod->outlet.GetVelocity() << std::endl;
            }

            CPPUNIT_ASSERT_EQUAL(mockheme->current_time, zerod->current_time);
          }
          MockHemeLB<MockIntercommunicator> *mockheme;
          Mock0DModel<MockIntercommunicator> *zerod;
          std::map<std::string, double> *pbuffer;
      };
      //class
      CPPUNIT_TEST_SUITE_REGISTRATION(MockIntercommunicatorTests);
    } //multiscale
  } //unittests
} //hemelb
#endif // HEMELB_UNITTEST_TEST_INTERCOMMUNICATOR_EXAMPLE_H
